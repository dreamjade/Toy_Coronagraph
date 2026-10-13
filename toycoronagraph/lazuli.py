"""
The Lazuli space telescope (UASAL; its configuration packages call it "Ultramarine", UM), computed with prysm.

    [[coronagraph]]
    type = "lazuli"
    charge = 0                  # 0: the telescope alone; even >= 2: a vortex mask and a Lyot stop
    wavelengths = ["632 nm"]    # "auto": [observation] wavelength
    throughput = "auto"         # weights from the mirror, filter and detector curves of the configs

The telescope and its wide-field coronagraph camera (WCC) come from the UASAL configuration
packages, read once (lazuli_parameters): config_um (focal length f_eff, mirror coatings) and
config_um_wcc (pupil sampling, beam size on M1, detector, filter). The pupil (pupil()) is a circle of
diameter D = the beam size on M1 (3.0 m) with an anti-aliased (gray) edge, filling pupil_samples x
pupil_samples samples (x along the columns, y along the rows, centred on the array).
OPD: an optional wavefront error map (opd) adds Lazuli's aberrations. Without one (the default) the
optics are ideal: a perfect telescope and, with a charge, an ideal coronagraph. The mirror OPD maps
of the configs ($UASAL_ARCHIVE) are not bundled.

A point source at angle (theta_x, theta_y), one wavelength at a time (LazuliBackend.source):
    pupil field   A exp(2 pi i (OPD + theta_x x + theta_y y) / lambda);
    coronagraph   charge >= 2: an ideal vortex (prysm's multi-resolution focal-plane propagation)
                  then a circular Lyot stop of lyot_mask_size D; charge 0: nothing;
    image         prysm's matrix DFT (MDFT) onto the image pixels; a pixel is pixel_scale scene
                  lambda/D, i.e. pixel_scale (wavelength0 / diameter0) f_eff in the focal plane
                  (only the angle matters).
The image is the weighted sum over the wavelengths. The image size must be odd: the star then sits
on the central pixel, which is prysm's centre n//2 as well as toycoronagraph's (n-1)/2, so the images
compare with prysm's pixel by pixel. Intensities follow toycoronagraph.optics: each pixel holds
|field|^2 / sum(A^2), the fraction of the source's light entering the aperture (prysm's executor
is unitary), so an unocculted source sums to 1 over the whole focal plane at every wavelength.
The focal field of P pupil samples repeats every P lambda/D (Lazuli's, at each wavelength): every
source is repeated P lambda/D away along x and y. optics_from_table refuses a sampling whose period
at the shortest wavelength does not clear the PSF frames (pupil samples are cheap; aliasing is not).

The Lazuli team's get_psf is here too (get_psf, on prysm 0.22, its inputs from tele_data), with
the wavelengths in um, prysm's unit (their version multiplies them by 1e-3, so it needs nm). It is
charge 0 of the backend:
    tele_data['Aperture'], ['dx'], ['fnum']   pupil(), its sample D / pupil_samples, f_eff / D
    tele_data['OPD']                          opd (a map in nm)
    tele_data['Pixsize'] / oversample         the image pixel: [grid] pixel_scale (detector_pixel_scale()
                                              gives the WCC detector's)
    tele_data['Wavelengths']                  wavelengths: um there; m here, or text ("0.632 um", "632 nm")
    tele_data['Throughput']                   throughput ("auto": spectral_weights)
    shift                                     moves get_psf's focal grid by +shift; here the source moves
                                              instead, to -shift
get_psf normalises each wavelength's PSF and the weighted sum to 1 over its own image; the backend
normalises each wavelength to the light entering the aperture, which keeps the light falling outside
the image (and the coronagraph's throughput).

This is the only module importing prysm and the config packages; toycoronagraph imports it when a
"lazuli" coronagraph is used (Optics.lazuli names it as its module).
"""
import os
import warnings
from functools import lru_cache
from types import SimpleNamespace
import numpy as np
from astropy import units as u
from astropy.io import fits
from prysm.geometry import antialias, circle_sdf
from prysm.propagation import (Wavefront, prepare_executor, prepare_multiresolution, unit_cell_focal_grid,
                               to_fpm_and_back_multiresolution, focus_dft, vortex_phase_mask)
from toycoronagraph.params import AUTO, RAD_TO_ARCSEC, ParamError
from toycoronagraph.optics import Optics, register_backend, MIN_SAMPLING, SamplingWarning
from toycoronagraph.psf import extended_frame_size

# Vortex propagation (prysm.propagation.prepare_multiresolution): the coarsest level spans the whole
# focal plane at Nyquist (unit_cell_focal_grid), each finer level zooms 4x into the singularity. Against
# finer settings (coarse level at 4 samples per lambda/D, 5 levels of 512 samples) the PSFs differ by
# < 5e-6 of a source's light (charges 2 and 6, pupil_samples 128 and 256), and the on-axis leak (3e-5 at
# 128 samples, 2e-6 at 512: the gray pupil edge) by ~5% of itself.
_VORTEX_LEVELS = 4
_VORTEX_FINE_SAMPLES = 256


def _quantity(value, unit):
    """
    A "parsed" config value in a unit.

    Args:
        value (dict): {"value": number, "unit": astropy unit text}, as config_um's load_config_values("parsed")
            gives it.
        unit (str): the wanted unit, e.g. "m".

    Returns:
        float: the value in unit.
    """
    return float(u.Quantity(value["value"], value["unit"]).to_value(unit))


def _csv_curve(path, header_lines):
    """
    A transmission curve from a CSV file of wavelengths in nm and values.

    Args:
        path (str): the file, e.g. config_um's coatings/NIST_1st_surface_Al.csv ("Wavelength nm,Specular Reflectance").
        header_lines (int): lines to skip before the data.

    Returns:
        tuple: wavelengths ((n,) float64, m) and values ((n,) float64, 0-1).
    """
    data = np.loadtxt(path, delimiter=",", skiprows=header_lines)
    return data[:, 0]*1e-9, data[:, 1]


def _fits_curve(path, unit):
    """
    A transmission curve from a FITS table with WAVELENGTH and THROUGHPUT columns (first extension).

    Args:
        path (str): the file.
        unit (str): the wavelength unit when the WAVELENGTH column declares none (TUNIT).

    Returns:
        tuple: wavelengths ((n,) float64, m) and values ((n,) float64).
    """
    with fits.open(path) as h:
        table = h[1]
        lam = table.data["WAVELENGTH"]*u.Unit(table.columns["WAVELENGTH"].unit or unit)
        return lam.to_value(u.m).astype(float), np.asarray(table.data["THROUGHPUT"], dtype=float)


@lru_cache(maxsize=None)
def lazuli_parameters():
    """
    The Lazuli telescope and WCC parameters, read once from config_um and config_um_wcc.

    Returns:
        SimpleNamespace:
            diameter (float): pupil diameter D in m: the beam size on M1 (config_um_wcc
                E2ES_default.optics.m1.beam_size; M1's clear aperture is larger, 3.065 m).
            focal_length (float): effective focal length in m (config_um telescope.general.f_eff, at the
                centre of the field; config_um's f_number = 15 is a provisional design requirement).
            pupil_samples (int): pupil samples across D (config_um_wcc E2ES_default.simulation.beam_sampling).
            pixel_size (float): WCC detector pixel size in m (config_um_wcc common_params.sensor.pixel_size).
            curves (dict): {name: (wavelengths (n,) float64 in m, values (n,) float64)} for "m1" (aluminium
                reflectance), "m2", "m3", "m4" (silver reflectance), "filter" (SDSS r') and "qe" (detector).
    """
    import config_um
    import config_um_wcc
    with warnings.catch_warnings():  # the mirror OPD map paths need $UASAL_ARCHIVE, not used here
        warnings.filterwarnings("ignore", message=".*UASAL_ARCHIVE")
        telescope = config_um.load_config_values("parsed")["telescope"]
        wcc = config_um_wcc.load_config_values("parsed")
    um_data, wcc_data = config_um.get_data_path(), config_um_wcc.get_data_path()
    mirrors = telescope["optics"]
    curves = {"m1": _csv_curve(os.path.join(um_data, mirrors["m1"]["coating_refl"]), 1)}
    for m in ("m2", "m3", "m4"):  # the silver table gives no wavelength unit; its 300-2500 range is nm
        curves[m] = _fits_curve(os.path.join(um_data, mirrors[m]["coating_refl"]), "nm")
    curves["filter"] = _fits_curve(os.path.join(wcc_data, wcc["WCC_ETC"]["filters"]["filter"]), "nm")
    # the config names ..._QE_curve.fits, but the package ships the same curve as .csv (nm, QE; no header)
    qe = os.path.splitext(os.path.join(wcc_data, wcc["common_params"]["sensor"]["qe"]))[0] + ".csv"
    curves["qe"] = _csv_curve(qe, 0)
    return SimpleNamespace(diameter=_quantity(wcc["E2ES_default"]["optics"]["m1"]["beam_size"], "m"),
                           focal_length=_quantity(telescope["general"]["f_eff"], "m"),
                           pupil_samples=int(wcc["E2ES_default"]["simulation"]["beam_sampling"]),
                           pixel_size=_quantity(wcc["common_params"]["sensor"]["pixel_size"], "m"), curves=curves)


def detector_pixel_scale():
    """
    The angular size of a WCC detector pixel, pixel_size / focal_length (config_um's f_eff = 42.047 m:
    0.018445 arcsec): e.g. for [grid] pixel_scale. Inside the backend the focal length cancels (only
    angles matter).

    Returns:
        float: arcsec per detector pixel.
    """
    t = lazuli_parameters()
    return t.pixel_size/t.focal_length*RAD_TO_ARCSEC


def spectral_weights(wavelengths):
    """
    The relative weight of each wavelength through Lazuli and the WCC: M1 aluminium x M2-M4 silver
    reflectance x SDSS r' filter x detector QE, interpolated at the wavelengths (0 outside a curve),
    times the width of the band each wavelength stands for (halfway to its neighbours, as far again
    beyond the first and last: np.gradient of the sorted wavelengths). The weighted sum is then the
    band integral however the band is sampled; evenly spaced wavelengths get the curves' product alone.

    Args:
        wavelengths (sequence of float): wavelengths in m, in any order (one wavelength: weight 1).

    Returns:
        np.ndarray: (n,) float64 weights, summing to 1.

    Raises:
        ValueError: no light at any of the wavelengths (all outside the r' band, about 540-730 nm).
    """
    lam = np.asarray(wavelengths, dtype=float)
    weights = np.prod([np.interp(lam, x, y, left=0, right=0) for x, y in lazuli_parameters().curves.values()],
                      axis=0)
    if lam.size > 1:
        order = np.argsort(lam)
        width = np.empty_like(lam)
        width[order] = np.gradient(lam[order])
        weights = weights*width
    if not weights.sum() > 0:
        raise ValueError("Lazuli: no light through the WCC's r' filter at %s nm: give wavelengths in the band "
                         "(about 540-730 nm), or the throughput explicitly" % ", ".join("%.4g" % (x*1e9) for x in lam))
    return weights/weights.sum()


def pupil(samples):
    """
    The Lazuli pupil: a circle of diameter D (the beam size on M1) with an anti-aliased (gray) edge,
    centred on the array (so the circle touching its edges is symmetric; prysm's grids put 0 at n//2).

    Args:
        samples (int): Pupil samples across D.

    Returns:
        tuple: x, y ((P, P) float64 pupil coordinates in mm; x along the columns, y along the rows) and
            the aperture ((P, P) float64 amplitude, 0 to 1).
    """
    D = lazuli_parameters().diameter*1e3
    line = (np.arange(samples)-(samples-1)/2)*D/samples
    x, y = np.meshgrid(line, line)
    return x, y, antialias(circle_sdf(D/2, np.hypot(x, y)), D/samples)


def read_opd(path):
    """
    An OPD (wavefront error) map over the pupil diameter.

    Args:
        path (str): a FITS file (primary HDU, in its BUNIT, converted with astropy; nm without BUNIT) or a
            .npy file (nm). The map is square, its rows along y and its columns along x.

    Returns:
        np.ndarray: (n, n) float64 OPD in nm (NaN, e.g. outside the pupil, set to 0).

    Raises:
        ValueError: the map is not square and 2-D, or its BUNIT is not a length.
    """
    if path.endswith(".npy"):
        opd = np.load(path).astype(float)
    else:
        with fits.open(path) as h:
            opd = h[0].data.astype(float)*u.Unit(h[0].header.get("BUNIT", "nm")).to(u.nm)
    if opd.ndim != 2 or opd.shape[0] != opd.shape[1]:
        raise ValueError("%s: an OPD map must be a square 2-D array (got shape %s)" % (path, opd.shape))
    return np.nan_to_num(opd)


def pupil_samples(opd, requested=AUTO):
    """
    The number of pupil samples across the aperture.

    Args:
        opd (str): the OPD map file; "" for none.
        requested (int or str): [[coronagraph]] pupil_samples; "auto": the OPD map's size, else the config's.

    Returns:
        int: the OPD map's size if there is one, else requested, else config_um_wcc's beam_sampling.

    Raises:
        ValueError: requested differs from the OPD map's size, or the map is not valid (see read_opd).
        OSError: the OPD map cannot be read.
    """
    if opd:
        n = read_opd(opd).shape[0]
        if requested not in (AUTO, n):
            raise ValueError("pupil_samples = %d but the OPD map %s is %d x %d: set pupil_samples = \"auto\""
                             % (requested, opd, n, n))
        return n
    return lazuli_parameters().pupil_samples if requested == AUTO else requested


def optics_from_table(values, img_pixel, psf_range, wavelength, diameter, size_fix):
    """
    The Optics of a "lazuli" [[coronagraph]] table, its "auto" settings resolved.

    Args:
        values (dict): the table's values, every key of params.CORONAGRAPHS["lazuli"] present (defaults
            filled in): charge, lyot_mask_size, wavelengths (list in m or "auto"), throughput (list or "auto"),
            opd (str), pupil_samples (int or "auto").
        img_pixel (int): Number of pixels on one edge of the image (odd).
        psf_range (float): Half-width of the image, in scene lambda/D.
        wavelength (float): The scene's wavelength in m ([observation] wavelength; "auto" wavelengths use it).
        diameter (float): The scene's telescope diameter in m ([observation] diameter).
        size_fix (str): How the user makes the image size odd, with two %d for the odd sizes around
            img_pixel, e.g. "set [grid] size = %d or %d (the image is resampled onto it)".

    Returns:
        Optics: see Optics.lazuli.

    Raises:
        ParamError: img_pixel is even, or the focal field's period (pupil samples x Lazuli's lambda/D at the
            shortest wavelength) does not clear the PSF frames: the farthest pixel of a star-centred radial
            frame (psf.extended_frame_size) from a source at the image edge, (N-1)/2 + (M-1)/2 pixels.

    Warns:
        SamplingWarning: fewer than MIN_SAMPLING x 8/3 pupil samples per lambda/D of the image's half-width
            (the PSFs' far wings come out too low near the image edge; see optics.MIN_SAMPLING).
        ValueError: "auto" throughput with no light at the wavelengths (spectral_weights), an OPD map that
            is not square and 2-D or whose BUNIT is not a length (read_opd), or pupil_samples differing
            from the OPD map's size.
        OSError: the OPD map cannot be read (e.g. FileNotFoundError).
    """
    if img_pixel % 2 == 0:
        raise ParamError("[[coronagraph]] lazuli needs an odd image size (got %d) so the star sits on a pixel as "
                         "in prysm: %s" % (img_pixel, size_fix % (img_pixel-1, img_pixel+1)))
    wavelengths = [wavelength] if values["wavelengths"] == AUTO else values["wavelengths"]
    weights = spectral_weights(wavelengths) if values["throughput"] == AUTO else values["throughput"]
    samples = pupil_samples(values["opd"], values["pupil_samples"])
    period = samples*min(wavelengths)/lazuli_parameters().diameter/(wavelength/diameter)  # scene lambda/D
    reach = ((img_pixel-1)/2 + (extended_frame_size(img_pixel)-1)/2)*2*psf_range/img_pixel  # scene lambda/D
    if period <= reach:
        needed = int(samples*reach/period)+1
        fix = ("use an OPD map of at least %d x %d samples" % (needed, needed) if values["opd"] else
               "set pupil_samples >= %d" % needed)
        raise ParamError("[[coronagraph]] lazuli: with %d pupil samples the focal field repeats every %.4g lambda/D "
                         "(at %.4g nm), within the %.4g lambda/D between a source and the farthest pixel of the PSF "
                         "frames, so every source would also appear that far away: %s, or a smaller image"
                         % (samples, period, min(wavelengths)*1e9, reach, fix))
    per_lambda_d = period/psf_range  # pupil samples per lambda/D (Lazuli's, shortest wavelength) of the half-width
    if per_lambda_d < MIN_SAMPLING*8/3:
        warnings.warn(SamplingWarning(
            "[[coronagraph]] lazuli: %d pupil samples give %.3g samples per lambda/D of the image's %.4g lambda/D "
            "half-width (at %.4g nm): the PSFs' far wings come out a few %% too low near the image edge. For about "
            "1%%, set pupil_samples >= %d, or use a smaller image"
            % (samples, per_lambda_d, psf_range, min(wavelengths)*1e9, int(np.ceil(samples*MIN_SAMPLING*8/3/per_lambda_d)))),
            stacklevel=2)
    return Optics.lazuli(values["charge"], img_pixel, psf_range, wavelengths, weights, wavelength, diameter, samples,
                         values["lyot_mask_size"], values["opd"])


@register_backend('lazuli')
class LazuliBackend:
    """
    The Lazuli telescope, alone (charge 0) or with an ideal vortex coronagraph (see the module docstring).

    A source enters only through the tilt across the pupil, whose width is D: for a fixed image pixel
    the field is a Fourier integral over the pupil, so the intensity is band-limited in the source
    angle at D / lambda cycles per radian, i.e. (D / diameter0) (wavelength0 / lambda) cycles per scene
    lambda/D. Without an OPD map, the pupil and the stops are circles: rotationally symmetric PSFs.

    Attributes:
        optics (Optics): kind 'lazuli' (see Optics.lazuli).
        rotationally_symmetric (bool): True without an OPD map.
        bandlimit (float): (D / diameter0) (wavelength0 / shortest wavelength), cycles per scene lambda/D.
        lambda_d (float): the scene's lambda/D, wavelength0 / diameter0, in radians.
        pupil_dx (float): pupil sample, in mm.
        efl (float): focal length (config_um's f_eff), in mm.
        um_per_lambda_d (float): focal-plane microns per scene lambda/D.
        x, y (np.ndarray): (P, P) float64 pupil coordinates in mm (x along the columns, y along the rows).
        aperture (np.ndarray): (P, P) float64 circle of diameter D (anti-aliased edge).
        opd (np.ndarray): (P, P) float64 OPD in nm (zeros without a map).
        lyot_stop (np.ndarray or None): (P, P) float64 circle of diameter lyot_mask_size D; None for charge 0.
        entering (float): sum(aperture^2), the light entering the aperture (the intensity unit).
        lyot_light (float): sum((aperture lyot_stop)^2), a source's light after the Lyot stop without the vortex
            (the throughput's reference; for charge 0, entering).
        vortex (callable or None): prysm's vortex_phase_mask(charge); None for charge 0.
        channels (list): per wavelength, (wavelength in um, weight, multi-resolution vortex executor or None).
        frame_executors (list): per wavelength, the MDFT onto the whole image.
    """

    def __init__(self, optics):
        """
        Args:
            optics (Optics): An Optics of kind 'lazuli' (see Optics.lazuli).
        """
        self.optics = optics
        p = dict(optics.params)
        telescope = lazuli_parameters()
        n, D = p['pupil_samples'], telescope.diameter*1e3  # samples, mm
        self.lambda_d = p['scene_wavelength']/p['scene_diameter']
        self.pupil_dx = D/n
        self.efl = telescope.focal_length*1e3
        self.um_per_lambda_d = self.lambda_d*self.efl*1e3
        self.rotationally_symmetric = not p['opd']
        self.bandlimit = telescope.diameter/min(p['wavelengths'])*self.lambda_d
        self.x, self.y, self.aperture = pupil(n)
        r = np.hypot(self.x, self.y)
        self.opd = read_opd(p['opd']) if p['opd'] else np.zeros((n, n))
        self.entering = np.sum(self.aperture**2)
        self.lyot_stop, self.lyot_light, self.vortex = None, self.entering, None
        if p['charge']:
            self.lyot_stop = antialias(circle_sdf(p['lyot_mask_size']*D/2, r), self.pupil_dx)
            self.lyot_light = np.sum((self.aperture*self.lyot_stop)**2)
            self.vortex = vortex_phase_mask(p['charge'])
        self.channels = [(lam*1e6, weight, self._multiresolution(lam*1e6) if self.vortex else None)
                         for lam, weight in zip(p['wavelengths'], p['weights'])]
        self.frame_executors = [self.executor(lam, (0.0, 0.0), (optics.img_pixel, optics.img_pixel))
                                for lam, _, _ in self.channels]

    def _multiresolution(self, wavelength):
        """
        prysm's multi-resolution executor for the vortex mask at one wavelength.

        Args:
            wavelength (float): in um.

        Returns:
            prysm.propagation.MultiResolutionExecutor: _VORTEX_LEVELS levels, the coarsest spanning the whole
                focal plane at Nyquist, the finer ones _VORTEX_FINE_SAMPLES across.
        """
        n = self.aperture.shape[0]
        focal_dx, focal_samples = unit_cell_focal_grid(self.pupil_dx, n*self.pupil_dx, wavelength, self.efl)
        return prepare_multiresolution(self.pupil_dx, n, focal_dx, focal_samples, wavelength, self.efl,
                                       _VORTEX_LEVELS, fine_samples=_VORTEX_FINE_SAMPLES)

    def executor(self, wavelength, center, shape):
        """
        prysm's MDFT from the pupil onto a block of image pixels.

        Args:
            wavelength (float): in um.
            center (tuple): (x, y) of the block's centre, in scene lambda/D from the star.
            shape (tuple): (ny, nx) size of the block, in pixels.

        Returns:
            prysm.fttools.MDFT: its output pixel (row, col) is at x = center_x + (col - (nx-1)/2) pixel_scale,
                y = center_y + (row - (ny-1)/2) pixel_scale.
        """
        ny, nx = shape
        s = self.optics.pixel_scale
        # prysm puts output sample n//2 at the shift, toycoronagraph the block centre at (n-1)/2 (half a pixel
        # lower for even n): shift to where sample n//2 is
        shift = tuple((c+(n//2-(n-1)/2)*s)*self.um_per_lambda_d for c, n in zip(center, (nx, ny)))
        return prepare_executor(self.pupil_dx, self.aperture.shape, s*self.um_per_lambda_d, shape, wavelength,
                                self.efl, focal_shift=shift)

    def source(self, x, y):
        """
        A point source at (x, y) scene lambda/D, through the telescope and the coronagraph.

        Args:
            x, y (float): Source position, in scene lambda/D from the star (x along columns, y along rows).

        Returns:
            LazuliSource: Its exit fields (after the Lyot stop with a charge), their total light and the
                throughput.
        """
        tilt = (x*self.x+y*self.y)*self.lambda_d*1e6  # nm (pupil coordinates in mm, angles in rad)
        fields = []
        for wavelength, _, stack in self.channels:
            field = Wavefront.from_amp_and_phase(self.aperture, self.opd+tilt, wavelength, self.pupil_dx).data
            if stack is not None:
                field = self.lyot_stop*to_fpm_and_back_multiresolution(field, self.vortex, stack)
            fields.append(field)
        return LazuliSource(self, fields)

    def unocculted(self):
        """
        The star without the vortex (the Lyot stop kept when there is a charge): the reference of the
        contrast_peak unit.

        Returns:
            LazuliSource: The star through the telescope (and the Lyot stop); throughput 1.
        """
        stop = 1.0 if self.lyot_stop is None else self.lyot_stop
        return LazuliSource(self, [stop*Wavefront.from_amp_and_phase(self.aperture, self.opd, wavelength, self.pupil_dx).data
                                   for wavelength, _, _ in self.channels])


class LazuliSource:
    """
    A point source through the Lazuli telescope: its fields leaving the last pupil, one per wavelength,
    focused onto any block of the image pixel grid.

    Attributes:
        total (float): Total light over the whole focal plane (Parseval: the exit fields' light, weighted
            over the wavelengths), in the package's units.
        throughput (float): The source's light leaving the system relative to the same source without
            the vortex (1 for charge 0).
    """

    def __init__(self, backend, fields):
        """
        Args:
            backend (LazuliBackend): The backend (channels, entering, lyot_light, executor()).
            fields (list): Per wavelength, the (P, P) complex128 field leaving the last pupil.
        """
        self._backend = backend
        self._fields = fields
        self.total = sum(w*np.sum(np.abs(f)**2) for (_, w, _), f in zip(backend.channels, fields))/backend.entering
        self.throughput = self.total*backend.entering/backend.lyot_light

    def _image(self, executors):
        """
        The weighted sum over the wavelengths of |MDFT(field)|^2, in the package's units.

        Args:
            executors (list): Per wavelength, the MDFT onto the wanted pixels.

        Returns:
            np.ndarray: (ny, nx) float64 intensity (rows = y, columns = x).
        """
        b = self._backend
        return sum(w*np.abs(focus_dft(f, ex))**2 for (_, w, _), f, ex in zip(b.channels, self._fields, executors))/b.entering

    def block(self, center, shape):
        """
        Intensity on a block of the image pixel grid (pixel size optics.pixel_scale).

        Args:
            center (tuple): (x, y) of the block's centre, in scene lambda/D from the star.
            shape (tuple): (ny, nx) size of the block, in pixels.

        Returns:
            np.ndarray: (ny, nx) float64 intensity (rows = y, columns = x), in the package's units.
        """
        return self._image([self._backend.executor(lam, center, shape) for lam, _, _ in self._backend.channels])

    def frame(self):
        """np.ndarray: (img_pixel, img_pixel) float64 intensity on the whole image."""
        return self._image(self._backend.frame_executors)


def tele_data(wavelengths, opd="", samples=AUTO):
    """
    The inputs of get_psf for Lazuli, from the same configuration values as the backend.

    Args:
        wavelengths (sequence of float): Wavelengths in um (prysm's unit), e.g. [0.632].
        opd (str): An OPD map file (see read_opd); "" for none (ideal optics).
        samples (int or str): Pupil samples across D; "auto": the OPD map's size, else the configuration's.

    Returns:
        dict: Aperture ((P, P) float64 amplitude, pupil()), OPD ((P, P) float64, nm), Wavelengths (list of
            float, um), dx (float, pupil sample in mm), fnum (float, f_eff / D), Pixsize (float, WCC detector
            pixel in um), Throughput ((n,) float64 weights, spectral_weights).
    """
    t = lazuli_parameters()
    n = pupil_samples(opd, samples)
    return {"Aperture": pupil(n)[2], "OPD": read_opd(opd) if opd else np.zeros((n, n)),
            "Wavelengths": list(wavelengths), "dx": t.diameter*1e3/n, "fnum": t.focal_length/t.diameter,
            "Pixsize": t.pixel_size*1e6, "Throughput": spectral_weights(np.asarray(wavelengths)*1e-6)}


def get_psf(tele_data, oversample=1, shift=(0, 0), samples=256):
    """
    The Lazuli team's get_psf on prysm 0.22 (prepare_executor + focus_dft replace focus_fixed_sampling;
    numpy instead of cupy), with tele_data['Wavelengths'] in um (their version multiplies the
    wavelengths by 1e-3, so it needs them in nm).

    Args:
        tele_data (dict): Aperture ((P, P) amplitude), OPD ((P, P) nm), Wavelengths (list, um), dx (pupil
            sample, mm), fnum (focal length / (dx P)), Pixsize (focal-plane pixel, um), Throughput (a
            float, or one weight per wavelength); see tele_data().
        oversample (int): Pixels per Pixsize.
        shift (tuple): (x, y) shift of the focal grid, um.
        samples (int): Pixels on one edge of the image.

    Returns:
        np.ndarray: (samples, samples) float64 PSF (rows y, columns x; the star, for shift (0, 0), at
            pixel samples // 2): each wavelength's PSF normalised to 1 over the image, weighted by the
            throughput, and the sum normalised to 1.
    """
    D = tele_data['dx']*tele_data['Aperture'].shape[0]
    weights = np.broadcast_to(tele_data['Throughput'], len(tele_data['Wavelengths']))
    psf = 0
    for wavelength, weight in zip(tele_data['Wavelengths'], weights):
        wavefront = Wavefront.from_amp_and_phase(tele_data['Aperture'], tele_data['OPD'], wavelength, tele_data['dx'])
        executor = wavefront.prepare_executor(tele_data['fnum']*D, tele_data['Pixsize']/oversample, samples, shift=shift)
        flux = np.abs(wavefront.focus_dft(executor).data)**2
        psf = psf + weight*flux/flux.sum()
    return psf/psf.sum()
