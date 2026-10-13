"""
The CME module: a coronal mass ejection seen in starlight scattered by its electrons.

The CME's shape is the model of astro_toolbox's cme_imaging (by iveydavis,
https://github.com/iveydavis/astro_toolbox; ported here with the author's agreement):
cme_mass_fraction() builds the share of the CME's mass in each pixel of a dim x dim box (dim
odd, the star at pixel dim // 2), a crescent-shaped shell pointing along +x. Its electrons
(mass / (mu_e x atomic mass unit)) scatter the starlight they receive at their projected
distance r from the star, so each pixel's contrast (a fraction of the star's flux) is

    N_e sigma_T / (4 pi r^2),

zero for pixels inside the star. Thomson scattering is grey: the contrast holds in any band.

With dim and extent "auto", the box is the scene's own grid (no resampling); otherwise its
dim pixels span extent au and the image is resampled onto the scene grid.

Newer cme_imaging versions (e.g. the halo/toroidal model) can be used directly instead: set
toolbox_path (the model is then built by cme_imaging with cme_args), or pass an existing
cme_imaging.Star with Cme.from_object(). Its contrast is read as grid_cme_photons * pix_res^2
over the star's photon output, taken from the model's own stellar flux grid (flux(r) x 4 pi r^2).
(cme_imaging's luminosity is 4 pi R^2 x radiance, without the factor pi of the surface flux, so
its luminosity / photon_energy is pi too small; a note says so when the two differ.)
"""
import io
import os
import sys
import importlib
import contextlib
import warnings
import numpy as np
from scipy.ndimage import gaussian_filter
from astropy import units as u, constants as const
from toycoronagraph.params import Param, ParamError, AUTO, _QUANTITY
from toycoronagraph.sources import Source, register_source
from toycoronagraph.units import SourceImage


# ---------------------------------------------------------------- the geometry (ported from cme_imaging)

def cme_mass_fraction(dim, dR_factor=0.2, vertical_height_factor=1.0, horizontal_height_factor=0.0, smear=False,
                      smear_sig=3.0):
    """
    The share of the CME's mass in each pixel: cme_imaging's CME.build_cme_normalized_grid,
    ported line by line (the same results, including its integer rounding).

    The CME is a shell between the radii Rmin and Rmax (in pixels; they grow with dR_factor and
    shift outwards with horizontal_height_factor), spanning +-45 degrees around +x, squeezed
    vertically by vertical_height_factor; its upper half is mirrored to the lower half.

    Method (on the upper-right quarter, rc = dim // 2 + 1 pixels on a side, dr = rc dR_factor):
    for each layer j < int(dr / 2), the edges r(phi) = R cos(2 phi), phi in [0, 45 deg], with
    R = Rmin_j = dr (1 + horizontal_height_factor) + j and Rmax_j = dr (2 + horizontal_height_factor) - j,
    are rounded to integer pixels; each column x = i is filled between the inner and outer edge
    (from y = 0 beyond the inner edge's reach; a column without an edge point takes the mean of
    the last edge sample found and the next one), each pixel holding the column's height, and the
    layer is added with the weight sqrt(|Rmin_j - Rmax_j| xmax) j (xmax: the outer edge's x
    samples, broadcast along the columns, as in cme_imaging; layer 0 has weight 0). Columns
    whose edges are not defined yet (or beyond the box) are skipped; y edges left from the
    previous column (or layer) are reused, as in cme_imaging.

    Do not restructure: tests/scene_stage4_cme.py checks it against cme_imaging bit for bit.

    Args:
        dim (int): Pixels on a side; an even dim is made odd (+1), as in cme_imaging. The box must
            hold some mass: at least 2 ceil(4 / dR_factor) - 1 px (39 at dR_factor = 0.2), so that
            int(dr / 2) >= 2 and layer 1 exists; a thin shell (small dR_factor) can need more,
            when all its filled columns get the weight 0 (xmax = 0).
        dR_factor (float): Shell thickness as a fraction of the half box (0 < dR_factor < 1).
        vertical_height_factor (float): Vertical stretch of the shell.
        horizontal_height_factor (float): Outward shift of the shell, in units of its thickness.
        smear (bool): Smooth the result with a Gaussian of smear_sig pixels.
        smear_sig (float): The Gaussian's sigma, in pixels.

    Returns:
        np.ndarray: (dim, dim) float64 array summing to 1; rows = y, columns = x, the star at
            pixel (dim // 2, dim // 2).

    Raises:
        ParamError: dR_factor is not between 0 and 1, or the box holds no mass (dim too small for
            dR_factor; cme_imaging then gives NaN everywhere): the message names a lower bound
            for dim (a thin shell can need more).
    """
    dim += dim % 2 == 0
    if not 0 < dR_factor < 1:
        raise ParamError("cme: dR_factor = %g must be between 0 and 1" % dR_factor)
    rc = int(dim/2) + 1                 # the quarter box, star at its pixel 0
    dr = rc*dR_factor                   # shell thickness, pixels
    grid = np.zeros((rc, rc))           # the quarter (rows y, columns x after the transpose below)
    phi = np.linspace(0, np.pi/4, rc)
    a = np.pi/2/(np.pi/4)               # r = R cos(a phi): 0 at phi = 45 deg
    for j in range(int(dr/2)):
        rmin_j = dr + j + dr*horizontal_height_factor
        rmax_j = dr*2 - j + dr*horizontal_height_factor
        rmin, rmax = rmin_j*np.cos(a*phi), rmax_j*np.cos(a*phi)
        # integer pixel coordinates of the inner and outer edges, from the tip (phi = 45 deg) inwards
        xmin = (rmin*np.cos(phi)).astype(int)[::-1]
        ymin = (vertical_height_factor*rmin*np.sin(phi)).astype(int)[::-1]
        xmax = (rmax*np.cos(phi)).astype(int)[::-1]
        ymax = (vertical_height_factor*rmax*np.sin(phi)).astype(int)[::-1]
        subset = np.zeros(grid.shape)   # indexed [x, y] (transposed when added)
        for i in range(int(rmax_j)):  # each column x = i: the shell's y range [yf_min, yf_max)
            idxmax, idxmin = np.where(xmax == i)[0], np.where(xmin == i)[0]
            if len(idxmax) != 0:  # the outer edge reaches this column
                idxmaxsave = idxmax[0]
                yf_max = ymax[idxmax[0]]
                if i > xmin.max():
                    yf_min = 0
                if len(idxmin) != 0:
                    yf_min = ymin[idxmin[0]]
                    idxminsave = idxmin[0]
                if len(idxmin) == 0 and i < xmin.max():
                    yf_min = int((ymin[idxminsave] + ymin[idxminsave + 1])/2)
            else:  # no outer edge in this column: interpolate between the last ones found
                if len(idxmin) != 0:
                    yf_max = int((ymax[idxmaxsave] + ymax[idxmaxsave + 1])/2)
                    yf_min = ymin[idxmin[0]]
                    idxminsave = idxmin[0]
                if len(idxmin) == 0 and i < xmin.max():
                    yf_max = int((ymax[idxmaxsave] + ymax[idxmaxsave + 1])/2)
                    yf_min = int((ymin[idxminsave] + ymin[idxminsave + 1])/2)
                if len(idxmin) == 0 and i >= xmin.max():
                    yf_max = int((ymax[idxmaxsave] + ymax[idxmaxsave + 1])/2)
                    yf_min = 0
            try:  # (cme_imaging skips columns whose edges are not defined)
                subset[i, yf_min:yf_max] = np.abs(yf_max - yf_min)
                if yf_max == yf_min:
                    subset[i, yf_max] = 1
            except Exception:
                pass
        grid += subset.T*np.abs((rmin_j - rmax_j)*xmax)**0.5*j   # (xmax broadcast along the columns)
    if not grid.any():  # no layer of non-zero weight: cme_imaging divides 0 by 0 below
        rc_min = int(np.ceil(4/dR_factor))      # the smallest rc with int(dr/2) >= 2, so that layer 1 exists
        rc_min += int(rc_min*dR_factor/2) < 2   # (4/dR_factor may be rounded down)
        need = ("at least %d px (a thin shell can need more)" % (2*rc_min - 1) if rc < rc_min
                else "more than %d px (a thin shell)" % dim)
        raise ParamError('cme: a model box of %d px (dim; "auto": the [grid] size) holds no mass of the shell at '
                         "dR_factor = %g: it needs %s; use a larger grid (or set dim and extent) or a larger dR_factor"
                         % (dim, dR_factor, need))
    # the quarter is the upper right (rows >= rc-1); mirrored about the row of the star for the lower right
    full = np.zeros((dim, dim))
    full[rc-1:, rc-1:] = grid
    full[0:rc, rc-1:] = np.flip(grid, axis=0)
    if smear:
        full = gaussian_filter(full, sigma=smear_sig)
    return full/np.nansum(full)


def cme_contrast(fraction, mass_g, pixel_au, radius_au, mu_e=1.0):
    """
    Each pixel's contrast: N_e sigma_T / (4 pi r^2), r the projected distance from the star.

    N_e = fraction x mass_g / (mu_e m_u) electrons per pixel, all at the distance r of the pixel
    centre (the electrons are taken in the plane of the sky).

    Args:
        fraction (np.ndarray): (n, n) share of the mass per pixel (cme_mass_fraction), n odd, the
            star at pixel (n // 2, n // 2).
        mass_g (float): The CME's mass in g.
        pixel_au (float): Pixel size in au.
        radius_au (float): The star's radius in au; pixels whose centre is inside it are 0.
        mu_e (float): Mass per electron in atomic mass units (cme_imaging: 1).

    Returns:
        np.ndarray: (n, n) float64 contrast per pixel (a fraction of the star's flux).
    """
    n = fraction.shape[0]
    offsets = (np.arange(n) - n//2)*pixel_au
    r_cm = (np.hypot(*np.meshgrid(offsets, offsets))*u.au).to_value(u.cm)
    electrons = fraction*mass_g/(mu_e*const.u.cgs.value)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(r_cm > radius_au*u.au.to(u.cm), electrons*const.sigma_T.cgs.value/(4*np.pi*r_cm**2), 0.0)


# ---------------------------------------------------------------- cme_imaging itself (newer versions)

def load_cme_imaging(toolbox_path=""):
    """
    The cme_imaging module of astro_toolbox.

    Args:
        toolbox_path (str): Folder holding cme_imaging.py (~ expanded; added to the front of
            sys.path); "" = importable as it is.

    Returns:
        module: cme_imaging (imported with its warnings silenced).

    Raises:
        FileNotFoundError: no cme_imaging.py in toolbox_path.
        ImportError: cme_imaging cannot be imported.
    """
    if toolbox_path:
        path = os.path.abspath(os.path.expanduser(toolbox_path))
        if not os.path.exists(os.path.join(path, "cme_imaging.py")):
            raise FileNotFoundError("cme: no cme_imaging.py in toolbox_path = %r" % toolbox_path)
        if path not in sys.path:
            sys.path.insert(0, path)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return importlib.import_module("cme_imaging")
    except ImportError as err:
        raise ImportError("cme: cannot import cme_imaging (%s); set toolbox_path to the astro_toolbox folder" % err)


def model_arguments(args):
    """
    cme_args for cme_imaging's make_cme.

    Args:
        args (dict): {name: value}; text with a unit ("1e19 g", matched by params._QUANTITY)
            becomes an astropy Quantity; other values are passed as they are.

    Returns:
        dict: the arguments, ready to pass as make_cme(**arguments).
    """
    return {k: u.Quantity(v) if isinstance(v, str) and _QUANTITY.match(v) else v for k, v in args.items()}


def star_photons(model):
    """
    The star's photon output (photons/s) of a cme_imaging.Star, from its stellar flux grid:
    flux(r) x 4 pi r^2, the same at every r (the median of the finite values is taken).

    Args:
        model: a cme_imaging.Star after particle_to_photon_flux().

    Returns:
        float: photons per second.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (model.grid_stellar_photon_fluxes*4*np.pi*model.grid_distances**2).to_value(1/u.s)
    return float(np.nanmedian(np.where(np.isfinite(out), out, np.nan)))


def contrast_of(model):
    """
    The CME's contrast per pixel of a cme_imaging.Star.

    Args:
        model: a cme_imaging.Star after particle_to_photon_flux().

    Returns:
        np.ndarray: (dim, dim) float64 contrast per pixel (grid_cme_photons x pix_res^2 / the star's
            photon output); NaN where cme_imaging has no stellar flux (inside the star).
    """
    photons = (model.grid_cme_photons*model.pix_res**2).to_value(1/u.s)
    return np.asarray(photons, dtype=float)/star_photons(model)


# ---------------------------------------------------------------- the module

# parameters of the ported geometry, not used when cme_imaging itself builds the CME (toolbox_path)
GEOMETRY = ("mass", "dR_factor", "vertical_height_factor", "horizontal_height_factor", "smear", "smear_sig", "mu_e")


@register_source
class Cme(Source):
    """
    A CME in Thomson-scattered starlight, as a contrast (family "relative", valid in any band).

    Three ways to build it (see the module docstring): the ported geometry (default; uses
    the GEOMETRY parameters, extent and dim), cme_imaging itself with cme_args (toolbox_path
    set; then also requires the star's temperature), or an existing cme_imaging.Star
    (from_object(); requires nothing, its distance, radius, temperature and grid are provided).
    """
    kind = "cme"
    doc = "a coronal mass ejection (the cme_imaging model of astro_toolbox)"
    family = "relative"
    requires = ("star.distance", "star.radius")
    params = (
        Param("mass", 1e19, unit="g", doc="CME mass", min=0, exclusive_min=True),
        Param("dR_factor", 0.2, doc="shell thickness, as a fraction of the half box", min=0, max=1, exclusive_min=True),
        Param("vertical_height_factor", 1.0, doc="vertical stretch of the shell", min=0, exclusive_min=True),
        Param("horizontal_height_factor", 0.0, doc="outward shift of the shell, in shell thicknesses", min=0),
        Param("smear", False, doc="smooth the CME with a Gaussian", type="bool", advanced=True),
        Param("smear_sig", 3.0, doc="sigma of the smoothing, in model pixels", min=0, exclusive_min=True, advanced=True),
        Param("mu_e", 1.0, doc="mass per electron, in atomic mass units (cme_imaging: 1)", min=0, exclusive_min=True,
              advanced=True),
        Param("extent", AUTO, unit="AU", doc='width of the model box; "auto" = the scene field of view',
              auto=True, min=0, exclusive_min=True, advanced=True),
        Param("dim", AUTO, doc='pixels of the model box (odd); "auto" = the scene grid size', type="int",
              auto=True, min=3, advanced=True),
        Param("toolbox_path", "", doc='use cme_imaging from this folder instead (newer versions)', type="str",
              advanced=True),
        Param("cme_args", {}, doc='with toolbox_path: arguments of make_cme(), e.g. {mass = "1e19 g"}', type="table",
              advanced=True,
              help="Passed to cme_imaging's Star.make_cme(**cme_args), so they follow the cme_imaging version you "
                   "have (e.g. toroidal_height and halo in newer versions). Values with units are written as text, "
                   "e.g. mass = \"1e19 g\"; numbers and true/false are passed as they are."),
    )

    def __init__(self, **settings):
        """
        Check the parameters; with toolbox_path set, also require the star's temperature.

        Args:
            **settings: the module's parameters (see `params`).
        """
        super().__init__(**settings)
        self._object = None   # a cme_imaging.Star given from Python
        self._built = {}      # models built with cme_imaging, by (dim, extent, star)
        if self.settings["toolbox_path"]:
            self.requires = Cme.requires + ("star.temperature",)  # cme_imaging's Star is a blackbody

    @classmethod
    def from_object(cls, star):
        """
        The module for an existing cme_imaging.Star (any version).

        Args:
            star (cme_imaging.Star): with its CME made (particle_to_photon_flux() is called if it
                was not, its output silenced).

        Returns:
            Cme: the module, with default settings and no requirements.
        """
        source = cls()
        if getattr(star, "grid_cme_photons", None) is None:
            with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
                warnings.simplefilter("ignore")
                star.particle_to_photon_flux()
        source._object = star
        source.requires = ()  # everything comes from the object
        return source

    def inspect(self):
        """
        With a cme_imaging.Star (from_object): its distance, radius, temperature and grid.

        Returns:
            list of tuple: [(key, value, origin)] for star.distance (pc), grid.pixel_scale (arcsec),
                star.radius (R_sun) and star.temperature (K) when they are Quantities, grid.size;
                [] without a cme_imaging.Star.
        """
        m = self._object
        if m is None:
            return []
        distance = m.distance.to_value(u.pc)
        found = [("star.distance", distance), ("grid.pixel_scale", m.pix_res.to_value(u.au)/distance)]
        for key, attr, unit in [("star.radius", "radius", u.R_sun), ("star.temperature", "temp", u.K)]:
            if isinstance(getattr(m, attr, None), u.Quantity):
                found.append((key, getattr(m, attr).to_value(unit)))
        found.append(("grid.size", int(contrast_of(m).shape[0])))
        return [(key, value, "the cme_imaging Star") for key, value in found]

    def _box(self, grid, star):
        """
        The model box: "auto" = the scene grid (made odd, same pixel size).

        Args:
            grid (toycoronagraph.scene.Grid): the scene grid (size, pixel_scale arcsec).
            star (toycoronagraph.scene.Star): its distance (pc).

        Returns:
            tuple: (dim (int, odd), extent (float, au): the box's width).
        """
        s = self.settings
        dim = grid.size if s["dim"] == AUTO else s["dim"]
        dim += dim % 2 == 0
        extent = grid.pixel_scale*star.distance*dim if s["extent"] == AUTO else s["extent"]
        return dim, extent

    def _external(self, grid, star):
        """
        A cme_imaging.Star built with cme_args (toolbox_path set), cached per box and star.

        Args:
            grid (toycoronagraph.scene.Grid): the scene grid.
            star (toycoronagraph.scene.Star): temperature (K), radius (R_sun), distance (pc).

        Returns:
            cme_imaging.Star: after make_cme(**cme_args) and particle_to_photon_flux().

        Raises:
            ParamError: make_cme rejects cme_args.
        """
        s = self.settings
        dim, extent = self._box(grid, star)
        key = (dim, extent, star)
        if key not in self._built:
            ci = load_cme_imaging(s["toolbox_path"])
            args = dict(temp=star.temperature*u.K, radius=star.radius*u.R_sun, distance=star.distance*u.pc, dim=dim,
                        linear_extent=extent*u.au)
            with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    model = ci.Star(**args)                     # the published version: keyword arguments
                except TypeError:
                    model = ci.Star(*args.values())             # newer versions: positional arguments
                try:
                    model.make_cme(**model_arguments(s["cme_args"]))
                except (u.UnitsError, TypeError) as err:
                    raise ParamError("cme: make_cme(%s) failed (%s); give values with units as text, e.g. "
                                     'mass = "1e19 g"' % (s["cme_args"], err))
                model.particle_to_photon_flux()
            self._built[key] = model
        return self._built[key]

    def render(self, grid, star):
        """
        The CME's contrast per pixel.

        Args:
            grid (toycoronagraph.scene.Grid): The scene grid (size, pixel_scale in arcsec).
            star (toycoronagraph.scene.Star): The star (distance pc, radius R_sun; with toolbox_path
                also the temperature K); not used with from_object().

        Returns:
            SourceImage: (dim, dim) contrast_pixel on the model box (pixel_scale set in arcsec;
                resampled by the scene when it differs from the grid), with notes on unused
                parameters, cme_imaging's luminosity and NaN pixels.

        Raises:
            ParamError: A model box too small to hold the shell's mass (dim, or the grid size with
                dim "auto": < 39 px at dR_factor = 0.2), with the ported geometry or cme_imaging's;
                with the ported geometry, dR_factor = 1; with toolbox_path, make_cme rejects cme_args.
        """
        s, notes = self.settings, []
        if self._object is None and not s["toolbox_path"]:  # the ported geometry
            dim, extent = self._box(grid, star)
            fraction = cme_mass_fraction(dim, s["dR_factor"], s["vertical_height_factor"], s["horizontal_height_factor"],
                                         s["smear"], s["smear_sig"])
            radius_au = (star.radius*u.R_sun).to_value(u.au)
            contrast = cme_contrast(fraction, s["mass"], extent/dim, radius_au, s["mu_e"])
            pixel_scale = extent/dim/star.distance
            origin = "CME (%d px, %.3g au, %.3g g)" % (dim, extent, s["mass"])
        else:  # cme_imaging itself
            changed = [p.name for p in self.params if p.name in GEOMETRY and s[p.name] != p.default]
            if s["toolbox_path"] and changed:
                notes.append("toolbox_path is set: the CME comes from cme_imaging with cme_args; %s not used"
                             % ", ".join(changed))
            model = self._object if self._object is not None else self._external(grid, star)
            contrast = contrast_of(model)
            ratio = star_photons(model)/(model.luminosity/model.photon_energy).to_value(1/u.s)
            if abs(ratio-1) > 0.01:
                notes.append("cme_imaging's luminosity / photon_energy is %.4g times smaller than the star's photon output "
                             "from its flux grid (it lacks the factor pi of the surface flux); the contrast uses the "
                             "latter, so it is %.4g times below cme_imaging's contrast cell" % (ratio, ratio))
            bad = ~np.isfinite(contrast)
            if bad.all():  # cme_imaging's 0/0: its box holds no mass of the shell
                raise ParamError("cme: cme_imaging gave NaN everywhere: its model box of %d px holds no mass of the "
                                 "shell; use a larger grid (or a larger dim, or dR_factor, in cme_args)"
                                 % contrast.shape[0])
            if bad.any():
                notes.append("%d pixel(s) inside the star (NaN in cme_imaging) set to 0" % bad.sum())
                contrast = np.where(bad, 0.0, contrast)
            pixel_scale = model.pix_res.to_value(u.au)/model.distance.to_value(u.pc)
            origin = "cme_imaging CME (%d px, %.3g au)" % (contrast.shape[0], model.linear_extent.to_value(u.au))
        if contrast.shape[0] != grid.size and grid.size % 2 == 0:
            origin += " [the CME model has an odd number of pixels: an odd grid size avoids resampling]"
        return SourceImage(contrast, "contrast_pixel", pixel_scale=pixel_scale, origin=origin, notes=tuple(notes))
