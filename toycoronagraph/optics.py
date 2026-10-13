"""
Optical systems behind one interface.

The spatial PSF machinery (toycoronagraph.psf, toycoronagraph.spatial) only asks an optical
system for the image of a point source: its throughput, its total light, and its intensity
on blocks of the image pixel grid. So a different coronagraph, or a starshade, can be added
as a new backend without touching that machinery:

    @register_backend('my_system')
    class MySystemBackend:
        rotationally_symmetric = False  # True allows the fast 'rotate' mode
        bandlimit = None                # cycles per lambda/D of the intensity vs. source
                                        # position, if known (also needed by 'rotate')
        def __init__(self, optics): ...
        def source(self, x, y): ...     # returns an object with .throughput, .total,
                                        # .block(center, shape) and .frame()
        def unocculted(self): ...       # optional: the star without the starlight-suppression
                                        # element (same kind of object; the contrast_peak unit)

and then Optics('my_system', img_pixel, psf_range, (('name', value), ...), module='my.package').
The module is imported when the kind is not registered yet (e.g. in parallel worker processes).
Two systems come with the package: Optics.vortex (hcipy; VortexBackend below) and Optics.lazuli
(the Lazuli space telescope computed with prysm; toycoronagraph.lazuli, imported only when used).

Conventions shared by every backend:
    * positions are in lambda/D, with the star at the origin; image pixel (row, col) is at
      x = (col - c) * pixel_scale, y = (row - c) * pixel_scale with c = (img_pixel - 1) / 2;
    * intensities are in the package's units (hcipy intensity / normal_factor, as in
      radial_psf_series): each pixel holds the fraction of the source's light entering the
      aperture, so kernels, frames and totals are directly comparable;
    * throughput is the source's total light relative to the same source without the
      starlight-suppression element.
"""
import os
import hashlib
import numpy as np
from dataclasses import dataclass
from functools import lru_cache
from hcipy import (make_pupil_grid, make_uniform_grid, make_circular_aperture, evaluate_supersampled,
                   FraunhoferPropagator, VortexCoronagraph, Apodizer, Wavefront)


@dataclass(frozen=True)
class Optics:
    """
    Description of an optical system and of the image sampling (hashable, picklable).

    Attributes:
        kind (str): Name of a registered backend, e.g. 'vortex'.
        img_pixel (int): Number of pixels on one edge of the image.
        psf_range (float): Half-width of the image in lambda/D.
        params (tuple): The backend's parameters as sorted (name, value) pairs.
        module (str): Module that registers the backend (imported if needed).
    """
    kind: str
    img_pixel: int
    psf_range: float
    params: tuple = ()
    module: str = 'toycoronagraph.optics'

    @classmethod
    def vortex(cls, charge, img_pixel, psf_range, lyot_mask_size=0.95):
        """
        An ideal vortex coronagraph with a circular aperture and a circular Lyot stop.

        Args:
            charge (int): Vortex charge (topological charge of the focal-plane mask).
            img_pixel (int): Number of pixels on one edge of the image.
            psf_range (float): Half-width of the image, in lambda/D.
            lyot_mask_size (float): Lyot stop diameter relative to the aperture diameter.

        Returns:
            Optics: kind 'vortex' with params (('charge', int), ('lyot_mask_size', float)).
        """
        return cls('vortex', int(img_pixel), float(psf_range),
                   (('charge', int(charge)), ('lyot_mask_size', float(lyot_mask_size))))

    @classmethod
    def lazuli(cls, charge, img_pixel, psf_range, wavelengths, weights, scene_wavelength, scene_diameter,
               pupil_samples, lyot_mask_size=0.95, opd=""):
        """
        The Lazuli space telescope, alone or with an ideal vortex coronagraph, computed with prysm
        (toycoronagraph.lazuli; observe.optics_for builds it from a [[coronagraph]] table).

        Args:
            charge (int): 0 for the telescope alone; else the vortex charge (even).
            img_pixel (int): Number of pixels on one edge of the image.
            psf_range (float): Half-width of the image, in scene lambda/D (scene_wavelength / scene_diameter).
            wavelengths (sequence of float): Wavelengths of the PSF, in m.
            weights (sequence of float): Relative weight of each wavelength (normalised here to sum 1).
            scene_wavelength (float): The scene's wavelength in m ([observation] wavelength); with
                scene_diameter it sets the lambda/D unit of source positions and of psf_range.
            scene_diameter (float): The scene's telescope diameter in m ([observation] diameter).
            pupil_samples (int): Pupil samples across the aperture (the OPD map's size when there is one).
            lyot_mask_size (float): Lyot stop diameter relative to the aperture (used with a charge).
            opd (str): OPD map file (see lazuli.read_opd); "" for none. Stored as an absolute path, with
                the SHA-1 of the file (opd_sha1), so cached PSFs follow the map's content.

        Returns:
            Optics: kind 'lazuli' (module 'toycoronagraph.lazuli'); the wavelengths and weights are
                stored as tuples (hashable).
        """
        digest = ""
        if opd:
            opd = os.path.abspath(opd)
            with open(opd, "rb") as f:
                digest = hashlib.sha1(f.read()).hexdigest()[:16]
        total = float(sum(weights))
        params = dict(charge=int(charge), lyot_mask_size=float(lyot_mask_size), opd=opd, opd_sha1=digest,
                      pupil_samples=int(pupil_samples), scene_diameter=float(scene_diameter),
                      scene_wavelength=float(scene_wavelength), wavelengths=tuple(float(w) for w in wavelengths),
                      weights=tuple(float(w)/total for w in weights))
        return cls('lazuli', int(img_pixel), float(psf_range), tuple(sorted(params.items())), 'toycoronagraph.lazuli')

    def param(self, name):
        """
        The value of one backend parameter.

        Args:
            name (str): Parameter name, e.g. 'charge'.

        Returns:
            The value stored in params (KeyError if there is none).
        """
        return dict(self.params)[name]

    @property
    def pixel_scale(self):
        """float: Image pixel size, in lambda/D per pixel (2*psf_range/img_pixel)."""
        return 2*self.psf_range/self.img_pixel

    @property
    def centre(self):
        """float: Pixel coordinate of the star on both axes, (img_pixel-1)/2."""
        return (self.img_pixel-1)/2

    def backend(self):
        """The backend object of this system (cached per process; see get_backend)."""
        return get_backend(self)

    def settings(self):
        """
        Flat dict of plain values (for .npz files and cache checks).

        Returns:
            dict: optics_kind, optics_module, img_pixel, psf_range, and 'optics_<name>'
                for every backend parameter.
        """
        out = dict(optics_kind=self.kind, optics_module=self.module, img_pixel=self.img_pixel, psf_range=self.psf_range)
        out.update({'optics_'+name: value for name, value in self.params})
        return out

    @classmethod
    def from_settings(cls, settings):
        """
        Rebuild from settings() (or from a saved spatial PSF model).

        Args:
            settings (dict): Holds the keys of settings(); other keys are ignored. Values may
                be numpy arrays (as loaded from .npz): 0-d ones become Python scalars (.item()),
                1-D ones tuples (e.g. the wavelengths of Optics.lazuli).
                optics_module defaults to 'toycoronagraph.optics'.

        Returns:
            Optics: The described system (params sorted by name).
        """
        plain = lambda v: np.asarray(v).item() if np.ndim(v) == 0 else tuple(np.asarray(v).tolist())
        params = tuple(sorted((k[len('optics_'):], plain(v)) for k, v in settings.items()
                              if k.startswith('optics_') and k not in ('optics_kind', 'optics_module')))
        return cls(str(settings['optics_kind']), int(settings['img_pixel']), float(settings['psf_range']), params,
                   str(settings.get('optics_module', 'toycoronagraph.optics')))


# A PSF's far wing at r lambda/D from its source comes from pupil structure D/r across, so it needs about 13
# pupil samples across D per lambda/D of r to be right to about 1%. The vortex model (VortexBackend: 2 N samples
# over 1.5 D, with N = 2 x the half-width x the pixels per lambda/D) gets 8/3 x its pixels per lambda/D samples
# per lambda/D of the image's half-width, so it needs at least MIN_SAMPLING pixels per lambda/D (measured
# against the exact Airy pattern, largest error inside the image edge: 5.0% at 2.36, 2.0% at 4, 1.3% at 5, 0.9%
# at 6 pixels per lambda/D); Lazuli (its own pupil samples) needs MIN_SAMPLING x 8/3 = 13.3 per lambda/D of the
# half-width.
MIN_SAMPLING = 5


class SamplingWarning(UserWarning):
    """The pupil is sampled too coarsely for the PSFs' far wings near the image edge (see MIN_SAMPLING)."""


_BACKENDS = {} # kind -> backend class, filled by register_backend


def register_backend(kind):
    """
    Class decorator registering a backend under the given kind.

    Args:
        kind (str): Name used as Optics.kind.

    Returns:
        callable: The decorator; it returns the class unchanged.
    """
    def register(cls):
        """Store cls in _BACKENDS under kind and return it unchanged."""
        _BACKENDS[kind] = cls
        return cls
    return register


@lru_cache(maxsize=4)
def get_backend(optics):
    """
    The backend of an Optics, built once per process (each parallel worker keeps its own).

    Args:
        optics (Optics): The system; its module is imported if its kind is not registered yet.

    Returns:
        The backend instance (e.g. VortexBackend).

    Raises:
        ValueError: If the kind is still unknown after importing optics.module.
    """
    if optics.kind not in _BACKENDS:
        import importlib
        importlib.import_module(optics.module)
    if optics.kind not in _BACKENDS:
        raise ValueError("unknown optics kind %r (registered: %s)" % (optics.kind, sorted(_BACKENDS)))
    return _BACKENDS[optics.kind](optics)


class PupilPlaneSource:
    """
    A point source behind a pupil-plane system: holds the field leaving the last pupil
    (e.g. after the Lyot stop) and focuses it onto any block of the image pixel grid.

    Attributes:
        throughput (float): The source's light leaving the system relative to the same
            source without the starlight-suppression element.
    """

    def __init__(self, backend, exit_field, throughput):
        """
        Args:
            backend: The backend (needs optics, normal_factor, pupil_grid and prop).
            exit_field (hcipy.Wavefront): Field leaving the last pupil, on backend.pupil_grid.
            throughput (float): See the class attributes.
        """
        self._backend = backend
        self._field = exit_field
        self.throughput = throughput

    @property
    def total(self):
        """float: Total light over the whole (infinite) focal plane (Parseval), in the package's pixel units."""
        b = self._backend
        return self._field.power.sum()/(b.optics.pixel_scale**2*b.normal_factor)

    def block(self, center, shape):
        """
        Intensity on a block of the image pixel grid (pixel size optics.pixel_scale).

        Args:
            center (tuple): (x, y) of the block's centre, in lambda/D from the star.
            shape (tuple): (ny, nx) size of the block, in pixels.

        Returns:
            np.ndarray: (ny, nx) float64 intensity (rows = y, columns = x), in the package's units.
        """
        b = self._backend
        ny, nx = shape
        s = b.optics.pixel_scale
        grid = make_uniform_grid([nx, ny], [nx*s, ny*s], center=list(center))
        return FraunhoferPropagator(b.pupil_grid, grid)(self._field).intensity.to_dict()["values"].reshape(ny, nx)/b.normal_factor

    def frame(self):
        """np.ndarray: (img_pixel, img_pixel) float64 intensity on the whole image (= |field()|^2)."""
        return np.abs(self.field())**2

    def field(self):
        """np.ndarray: (img_pixel, img_pixel) complex field on the whole image (rows = y, columns = x)."""
        b = self._backend
        n = b.optics.img_pixel
        return np.asarray(b.prop(self._field).electric_field).reshape(n, n)/np.sqrt(b.normal_factor)


@register_backend('vortex')
class VortexBackend:
    """
    Ideal vortex coronagraph: circular aperture, vortex focal-plane mask, circular Lyot stop.

    Rotationally symmetric, and a source enters only through a tilt across the pupil, so the
    intensity is band-limited in the source position at D/lambda: 'rotate' mode is valid.

    Attributes:
        optics (Optics): The system (kind 'vortex', params charge and lyot_mask_size).
        pupil_grid: hcipy pupil grid, 2*img_pixel samples across 1.5 aperture diameters.
        aperture: Circular aperture of diameter 1 on pupil_grid (4x supersampled edges).
        prop (FraunhoferPropagator): Pupil to the (img_pixel, img_pixel) image grid in lambda/D.
        lyot_stop (Apodizer): Circular Lyot stop of diameter lyot_mask_size.
        coro (VortexCoronagraph): The vortex mask (with its own internal propagation).
        normal_factor (float): Divisor turning hcipy intensities into the package's units.
    """
    rotationally_symmetric = True
    bandlimit = 1.0

    def __init__(self, optics):
        """
        Args:
            optics (Optics): An Optics of kind 'vortex' (e.g. Optics.vortex(...)).
        """
        self.optics = optics
        n, psf_range = optics.img_pixel, optics.psf_range
        self.pupil_grid = make_pupil_grid(n*2, 1.5)
        self.aperture = evaluate_supersampled(make_circular_aperture(1), self.pupil_grid, 4)
        focal_grid_resolution = n/2/psf_range # pixels per lambda/D (= 1/pixel_scale)
        # Symmetric grid in lambda/D with the star at pixel (n-1)/2, i.e. between the four
        # central pixels for even n (make_focal_grid would always put a pixel on the star)
        self.prop = FraunhoferPropagator(self.pupil_grid, make_uniform_grid([n, n], [2*psf_range, 2*psf_range]))
        self.lyot_stop = Apodizer(evaluate_supersampled(make_circular_aperture(optics.param('lyot_mask_size')), self.pupil_grid, 4))
        self.coro = VortexCoronagraph(self.pupil_grid, optics.param('charge'))
        self.normal_factor = np.pi*(focal_grid_resolution/2)**2

    def wavefront(self, x, y):
        """
        A point source at (x, y) lambda/D, as a tilted wavefront on the aperture.

        Args:
            x, y (float): Source position, in lambda/D from the star (x along columns, y along rows).

        Returns:
            hcipy.Wavefront: The aperture field with the source's tilt, on pupil_grid.
        """
        return Wavefront(self.aperture*np.exp(2j*np.pi*(self.pupil_grid.x*x+self.pupil_grid.y*y)))

    def source(self, x, y):
        """
        A point source at (x, y) lambda/D, propagated through the vortex and the Lyot stop.

        Args:
            x, y (float): Source position, in lambda/D from the star.

        Returns:
            PupilPlaneSource: Its throughput is the light after the Lyot stop with the vortex
                relative to the light after the Lyot stop without it.
        """
        wf = self.wavefront(x, y)
        exit_field = self.lyot_stop(self.coro(wf))
        return PupilPlaneSource(self, exit_field, exit_field.power.sum()/self.lyot_stop(wf).power.sum())

    def unocculted(self):
        """
        The star without the vortex (the Lyot stop kept): the reference of the contrast_peak unit.

        Returns:
            PupilPlaneSource: The star through the aperture and the Lyot stop (throughput 1).
        """
        return PupilPlaneSource(self, self.lyot_stop(self.wavefront(0, 0)), 1.0)
