"""
A scene through one or more coronagraphs.

    result = observe(scene)                              # the parameter file's coronagraphs and settings
    result = observe(scene, ["vortex:2", "vortex:6"], combine="together")
    result.final("vortex_c6")                            # the total; result.final("vortex_c6", "wind"): one component
    result.save()                                        # results/<run>/: FITS per coronagraph, params_used.toml, plots

combine = "separate" (default): each component goes through each coronagraph with its own
method (symmetric ones on the exact circular method, the others on the spatial one), and the
final images are summed, each one kept. "together": the components are summed first and go
through each coronagraph once. Planets (point sources) always get their exact PSF.
components=[...] restricts a run to some components (the spatial kernel size still comes
from the whole scene, so separate runs add up to the combined one).

A Coronagraph holds one optical system on one grid: its PSFs (radial series and spatial
models, cached on disk and in memory) and the final images it has made (cached by image),
so only what changed is computed again. Coronagraph.get returns the same object for the same
optics and settings, so Target (toycoronagraph.main) and observe() share these caches.
optics_for turns a [[coronagraph]] table into its Optics, one builder per type (OPTICS_BUILDERS):
"vortex" (hcipy) and "lazuli" (the Lazuli telescope computed with prysm, toycoronagraph.lazuli).

A Result holds the final images in the scene's internal unit ("Jy/arcsec^2" or
"contrast_pixel"); Result.save writes a run folder that Result.load reads back.
Conventions: (N, N) images with the star at pixel (N-1)/2; rows are +y (as in FITS).
"""
import os
import hashlib
import importlib
import json
import warnings
from functools import cached_property
import numpy as np
from astropy.io import fits
from toycoronagraph import params as P
from toycoronagraph import units as U
from toycoronagraph import selection
from toycoronagraph.optics import Optics, MIN_SAMPLING, SamplingWarning
from toycoronagraph.psf import (cir_psf, cir_psf_planets, radial_psf_series, extended_frame_size, choose_kernel_size,
                                spatial_psf_calculation, exact_point_sources, load_spatial_psfs)
from toycoronagraph.spatial import spatial_psf, image_error_estimate


def _digest(array):
    """
    A cache key for an array's values.

    Args:
        array (array-like): Any array (converted to contiguous float64 first).

    Returns:
        str: The SHA-1 hex digest (40 characters) of its bytes.
    """
    return hashlib.sha1(np.ascontiguousarray(array, dtype=np.float64).tobytes()).hexdigest()


# Layout of the spatial model cache files (model_*.npz); bump it if what they hold changes
_MODEL_CACHE_FORMAT = 1


def _code_digest(*modules):
    """
    A cache key part that changes with the code: the SHA-1 of the modules' source files.

    Args:
        modules (str): Module names, e.g. "toycoronagraph.psf" (imported if needed).

    Returns:
        str: The first 16 hex characters of the SHA-1 of their files' bytes, in the given order.
    """
    digest = hashlib.sha1()
    for name in modules:
        with open(importlib.import_module(name).__file__, "rb") as f:
            digest.update(f.read())
    return digest.hexdigest()[:16]


def _read_model(path, description):
    """
    A spatial model saved by _write_model, if it was made from the same description.

    Args:
        path (str): The .npz file.
        description (str): The JSON description the model must have been saved with.

    Returns:
        dict or None: The model (see psf.spatial_psf_calculation); None if there is no file or its
            description differs.
    """
    if not os.path.exists(path):
        return None
    model = load_spatial_psfs(path)
    return model if model.pop("description", None) == description else None


def _write_model(path, description, model):
    """
    Save a spatial model with its description (written to a temporary file first, so an interrupted
    save leaves no partial model behind).

    Args:
        path (str): The .npz file.
        description (str): The JSON description of what the model was made from (see Coronagraph.spatial_model).
        model (dict): The model (arrays and plain values).
    """
    with open(path + ".part", "wb") as f:
        np.savez(f, description=description, **model)
    os.replace(path + ".part", path)


class Coronagraph:
    """
    One coronagraph on one grid: its PSFs and the final images of the images given to it.

    Attributes:
        optics (Optics): The optical system and image sampling (img_pixel N, psf_range in lambda/D).
        name (str): Name in outputs, e.g. "vortex_c6".
        num_cores, error_target, max_kernel_gb, use_gpu, cache_dir, symmetric, symmetry_tol: See __init__.
        rot_number (int): Circular method: number of rotations (N when 0 was given).
        computed (list): (what, key) of everything actually calculated, e.g. ("radial PSFs", name),
            ("spatial model", (mode, kernel_size)), ("final image", first 8 characters of the image's
            digest); what is not listed came from the caches.
    """
    _instances = {}  # (optics, name, settings) -> Coronagraph (see get)

    def __init__(self, optics, name, num_cores=-1, error_target=0.02, max_kernel_gb=4.0, use_gpu=False,
                 cache_dir="", symmetric=P.AUTO, rot_number=0, symmetry_tol=0.01):
        """
        Args:
            optics (Optics): The optical system and image sampling.
            name (str): Name in outputs and messages.
            num_cores (int): CPU cores for the PSF calculations (joblib: -1 = all).
            error_target (float): Target relative error of spatial-method final images (0.02 = 2%).
            max_kernel_gb (float): Memory limit of the spatial PSF kernels, in GB.
            use_gpu (bool): Spatial method: use CuPy for the FFT convolutions.
            cache_dir (str): Folder of the PSF cache files; "": the working directory.
            symmetric (bool or str): Whether the PSFs are circularly symmetric; "auto": from the
                optics model (or measured), see selection.psf_symmetry.
            rot_number (int): Circular method: number of rotations of the radial PSFs; 0: N.
            symmetry_tol (float): Largest measured image asymmetry still counted as symmetric.
        """
        self.optics, self.name = optics, name
        self.num_cores, self.error_target, self.max_kernel_gb, self.use_gpu = num_cores, error_target, max_kernel_gb, use_gpu
        self.cache_dir, self.symmetric, self.symmetry_tol = cache_dir, symmetric, symmetry_tol
        self.rot_number = rot_number or optics.img_pixel
        # in-memory caches: spatial models by (mode, kernel_size), kernel sizes by (mode, weight digest),
        # final images by (image digest, method, mode, kernel_size, iwa_ignore), and the PSF symmetry
        self._models, self._kernel_sizes, self._finals, self._psf_symmetry = {}, {}, {}, None
        self._radial_psfs = None  # the radial PSFs, memory-mapped once (see radial_psfs)
        self.computed = []  # (what, key): what was actually calculated (the rest came from the caches)

    @classmethod
    def get(cls, optics, name, **settings):
        """
        The same object for the same optics and settings (so its PSFs and final images are reused).

        Args:
            optics (Optics): The optical system and image sampling.
            name (str): Name in outputs.
            settings: The other __init__ arguments (keywords; values must be hashable).

        Returns:
            Coronagraph: The cached instance, created on first use.
        """
        key = (optics, name, tuple(sorted(settings.items())))
        if key not in cls._instances:
            cls._instances[key] = cls(optics, name, **settings)
        return cls._instances[key]

    # ---- PSFs
    def _cache_path(self, name):
        """
        The path of a cache file in cache_dir (the folder is created if needed).

        Args:
            name (str): The file name.

        Returns:
            str: cache_dir/name (./name if cache_dir is "").
        """
        folder = self.cache_dir or "."
        os.makedirs(folder, exist_ok=True)
        return os.path.join(folder, name)

    @property
    def tag(self):
        """str: Part of the cache file names: everything the PSFs depend on, e.g. "vortex_c6_N512_range16_lyot0.95"
        for the vortex (and backends with its parameters), else "<kind>_N<N>_range<psf_range>_<digest>" with
        the first 12 hex characters of the SHA-1 of the optics' params, e.g. "lazuli_N101_range3.16_0123456789ab"."""
        o, params = self.optics, dict(self.optics.params)
        if set(params) == {"charge", "lyot_mask_size"}:
            return "%s_c%d_N%d_range%.6g_lyot%.4g" % (o.kind, params["charge"], o.img_pixel, o.psf_range,
                                                      params["lyot_mask_size"])
        return "%s_N%d_range%.6g_%s" % (o.kind, o.img_pixel, o.psf_range,
                                        hashlib.sha1(repr(o.params).encode()).hexdigest()[:12])

    def radial_psfs_file(self):
        """
        The radial PSF series, computed once and cached on disk.

        PSFs of sources along +x at psf.radial_source_positions, on frames larger than the
        image (psf.extended_frame_size) so that rotated about the star they still cover the
        image corners.

        Returns:
            str: Path of the .npy file, cache_dir/radial_<tag>_frame<size>.npy, holding a
                ((N+1)//2, size, size) float32 array.
        """
        size = extended_frame_size(self.optics.img_pixel)
        path = self._cache_path("radial_%s_frame%d.npy" % (self.tag, size))
        if not os.path.exists(path):
            print("Calculating the radial PSFs of %s (%d frames of %d x %d) -> %s"
                  % (self.name, (self.optics.img_pixel+1)//2, size, size, path))
            np.save(path, radial_psf_series(self.optics, self.num_cores, size))
            self.computed.append(("radial PSFs", self.name))
        return path

    def _model_description(self, mode, kernel_size):
        """
        Everything a spatial model is made from, as JSON: the key of its cache file (see spatial_model).

        The optics in full (repr of the Optics: kind, image size and field, every backend parameter,
        e.g. Lazuli's wavelengths, weights, pupil samples and the OPD map's SHA-1), the mode, kernel size
        and error target (which sets the sample spacing), the radial PSF file it is built from in
        'rotate' mode (name, size and modification time: recomputing it rebuilds the model), and the
        code that builds it (the source of psf.py, spatial.py and the optics backend's module).

        Args:
            mode (str): "rotate" or "direct".
            kernel_size (int): Kernel width in pixels.

        Returns:
            str: JSON text (sorted keys).
        """
        description = dict(format=_MODEL_CACHE_FORMAT, optics=repr(self.optics), mode=mode, kernel_size=int(kernel_size),
                           error_target=float(self.error_target),
                           code=_code_digest("toycoronagraph.psf", "toycoronagraph.spatial", self.optics.module))
        if mode == "rotate":
            radial = os.stat(self.radial_psfs_file())
            description["radial_psfs"] = [os.path.basename(self.radial_psfs_file()), radial.st_size, radial.st_mtime_ns]
        return json.dumps(description, sort_keys=True)

    def radial_psfs(self):
        """
        The radial PSF series (radial_psfs_file), memory-mapped once: the circular method and the planets
        read only the frames they use, from the file (or the system's file cache), not the whole file
        every time (552 MB at 512 px).

        Returns:
            np.memmap: ((N+1)//2, size, size) float32 PSFs, read-only.
        """
        if self._radial_psfs is None:
            self._radial_psfs = np.load(self.radial_psfs_file(), mmap_mode="r")
        return self._radial_psfs

    def spatial_model(self, mode, weight):
        """
        The spatial PSF model whose kernel size meets the error target for an image.

        The kernel size is chosen for `weight` (psf.choose_kernel_size; in 'rotate' mode from the
        cached radial PSFs, without new propagations) and cached per image, in memory and on disk
        (cache_dir/kernel_<tag>_<16 hex characters>.json, keyed by the image, mode and error target);
        a model already built with a larger kernel is reused (it is more accurate). The model itself
        (kernels, sample grid, blend weights, transmission map and error report) is cached on disk,
        cache_dir/model_<tag>_<mode>_K<kernel_size>_<16 hex characters>.npz, keyed by everything it is
        made from (_model_description, also stored in the file and compared when it is read: a file
        made from anything else is built again). Building one also caches the x-axis PSF runs and the
        error report (psf.spatial_error_report) on disk (cache_dir/axis_<tag>_K<kernel_size>.npz,
        cache_dir/error_<tag>_<mode>_K<kernel_size>_t<target>.npz), so a new session estimates nothing again.

        Args:
            mode (str): "rotate" or "direct" (see psf.spatial_psf_calculation).
            weight (np.ndarray): (N, N) image whose light sets the kernel size.

        Returns:
            dict: The model (see psf.spatial_psf_calculation; "kernel_size" is its kernel width in pixels).
        """
        key = (mode, _digest(weight))
        radial = None  # the radial PSFs ('rotate' mode), loaded once if needed
        if key not in self._kernel_sizes:
            name = hashlib.sha1(("%s %s %r" % (key[1], mode, self.error_target)).encode()).hexdigest()[:16]
            path = self._cache_path("kernel_%s_%s.json" % (self.tag, name))
            if os.path.exists(path):
                with open(path) as f:
                    self._kernel_sizes[key] = json.load(f)["kernel_size"]
                print("Kernel size K=%d (chosen for this image before: %s)" % (self._kernel_sizes[key], path))
            else:
                radial = np.load(self.radial_psfs_file()) if mode == "rotate" else None
                self._kernel_sizes[key] = choose_kernel_size(self.optics, weight, self.error_target, self.num_cores,
                                                             radial_psfs=radial)[0]
                with open(path, "w") as f:
                    json.dump({"kernel_size": self._kernel_sizes[key]}, f)
        larger = sorted(k for m, k in self._models if m == mode and k >= self._kernel_sizes[key])
        kernel_size = larger[0] if larger else self._kernel_sizes[key]
        if (mode, kernel_size) not in self._models:
            description = self._model_description(mode, kernel_size)
            path = self._cache_path("model_%s_%s_K%d_%s.npz" % (self.tag, mode, kernel_size,
                                                               hashlib.sha1(description.encode()).hexdigest()[:16]))
            model = _read_model(path, description)
            if model is not None:
                print("Spatial PSF model K=%d: from %s" % (kernel_size, path))
                self._models[(mode, kernel_size)] = model
                return model
            if mode == "rotate" and radial is None:
                radial = np.load(self.radial_psfs_file())
            axis_cache = self._cache_path("axis_%s_K%d.npz" % (self.tag, kernel_size))
            report_cache = self._cache_path("error_%s_%s_K%d_t%g.npz" % (self.tag, mode, kernel_size, self.error_target))
            cached = os.path.exists(report_cache)
            model = spatial_psf_calculation(
                self.optics, self.num_cores, radial, kernel_size=kernel_size, error_target=self.error_target,
                pre_img=weight, mode=mode, axis_cache=axis_cache, max_kernel_gb=self.max_kernel_gb, validate=not cached)
            if cached:
                with np.load(report_cache) as f:
                    model.update({k: f[k] for k in f.files})
                print("Spatial PSF error analysis: from %s" % report_cache)
            else:
                np.savez(report_cache, **{k: v for k, v in model.items() if k.startswith("err_")})
            _write_model(path, description, model)
            print("Spatial PSF model K=%d saved (%.2f GB): %s" % (kernel_size, os.path.getsize(path)/1e9, path))
            self._models[(mode, kernel_size)] = model
            self.computed.append(("spatial model", (mode, kernel_size)))
        return self._models[(mode, kernel_size)]

    @cached_property
    def peak_norm(self):
        """
        float or None: The peak of the star's unocculted PSF (the backend's unocculted(): the optics without
        the starlight-suppression element, the Lyot stop kept), as a fraction of the star's light per pixel:
        the contrast_peak reference. It is the intensity at the star's position, so it does not depend on
        where the pixels fall (for an even N the star is between four pixels). Computed once; None if the
        backend has no unocculted().
        """
        backend = self.optics.backend()
        return float(backend.unocculted().block((0.0, 0.0), (1, 1))[0, 0]) if hasattr(backend, "unocculted") else None

    @cached_property
    def iwa(self):
        """
        float or None: The inner working angle, in lambda/D: the distance from the star at which a point source's
        throughput (its light leaving the system relative to the same source without the starlight-suppression
        element, so 1 far from the star) reaches 0.5. It depends on the optics only, not on the image: it is
        searched along +x beyond the image's edge if needed (to 64 lambda/D), by bisection to 0.01 lambda/D (about
        15 sources, assuming the throughput rises with the distance); computed once. None when a source on the
        star already passes half (no coronagraph: e.g. the Lazuli telescope alone, charge = 0), or none reaches it
        by 64 lambda/D.
        """
        backend = self.optics.backend()
        passes = lambda r: backend.source(r, 0.0).throughput >= 0.5
        if passes(0.0):
            return None
        lo, hi = 0.0, self.optics.psf_range
        while not passes(hi):
            if hi >= 64:
                return None
            lo, hi = hi, 2*hi
        while hi-lo > 0.01:
            mid = (lo+hi)/2
            lo, hi = (lo, mid) if passes(mid) else (mid, hi)
        return (lo+hi)/2

    def psf_symmetry(self):
        """
        Whether the PSFs are circularly symmetric (selection.psf_symmetry, computed once).

        Returns:
            tuple: (bool, reason str).
        """
        if self._psf_symmetry is None:
            self._psf_symmetry = selection.psf_symmetry(self.optics, self.symmetric, num_cores=self.num_cores)
        return self._psf_symmetry

    # ---- images
    def choose(self, image, symmetry=None, method=P.AUTO, spatial_mode=P.AUTO):
        """
        The method for an image (see selection.choose_method).

        Args:
            image (np.ndarray): (N, N) image (only used to measure its symmetry when `symmetry` is None).
            symmetry (tuple): (bool, reason str) if the image's symmetry is known; None: measured.
            method (str): "auto", "circular" or "spatial".
            spatial_mode (str): "auto", "rotate" or "direct".

        Returns:
            tuple: (method str, spatial_mode str or None for circular, messages: list of str, the first
                one the summary).

        Raises:
            ValueError: An unknown method / spatial_mode, or "rotate" for optics that do not allow it.
        """
        image_sym, image_why = symmetry or selection.image_symmetry(image, P.AUTO, self.symmetry_tol)
        psf_sym, psf_why = self.psf_symmetry()
        backend = self.optics.backend()
        # 'rotate' needs rotationally symmetric optics with a known band limit
        return selection.choose_method(method, spatial_mode, image_sym, psf_sym,
                                       bool(backend.rotationally_symmetric and backend.bandlimit), (image_why, psf_why))

    def image(self, image, method, spatial_mode=None, weight=None, iwa_ignore=False, report=True):
        """
        The final image of an image, with the circular or the spatial method (cached).

        The coronagraph is linear, so the image may be in any unit; the final image is in the same one.

        Args:
            image (np.ndarray): (N, N) image, star at pixel (N-1)/2.
            method (str): "circular" (exact for symmetric images and PSFs) or "spatial".
            spatial_mode (str): Spatial method: "rotate" or "direct" (ignored for circular).
            weight (np.ndarray): (N, N) image setting the spatial kernel size; None: `image`.
            iwa_ignore (bool): Leave out the light inside the inner working angle.
            report (bool): Spatial method: print the estimated error of this image.

        Returns:
            tuple: (final image: (N, N) float64 in the unit of `image`, estimated typical relative
                error: float, 0 for the exact circular method).
        """
        weight = image if weight is None else weight
        model = self.spatial_model(spatial_mode, weight) if method == "spatial" else None
        key = (_digest(image), method, spatial_mode, model["kernel_size"] if model else None, iwa_ignore)
        if key not in self._finals:
            if method == "circular":
                final = cir_psf(image, [], [], self.optics.pixel_scale, iwa_ignore, False, self.optics.img_pixel,
                                self.rot_number, self.radial_psfs())
                error = 0.0
            else:
                final = spatial_psf(image, model, iwa_ignore, self.use_gpu, self.num_cores, report=report)
                error = image_error_estimate(model, image, verbose=False).get("typical", 0.0)
            self._finals[key] = (final, error)
            self.computed.append(("final image", key[0][:8]))
        return self._finals[key]

    def point_sources(self, pixels, values):
        """
        Point sources through the coronagraph, each with its exact PSF (psf.exact_point_sources).

        Args:
            pixels (list): (row, column) of each source (floats; rows = +y in a scene, see Grid.pixel).
            values (list): Each source's flux per pixel area (e.g. Jy/arcsec^2 = Jy / pixel_scale^2).

        Returns:
            np.ndarray: (N, N) image of the sources, in the unit of `values`. Sources outside the
                image are left out (with a message).
        """
        return exact_point_sources(self.optics, pixels, values)

    def interpolated_point_sources(self, positions, values, pixel_scale):
        """
        Point sources from the radial PSFs, interpolated in radius and rotated (Target's circular method).

        Args:
            positions (list): (radius arcsec, angle rad) of each source; the orientation of
                tool.arcsec_to_pixel (y increases with the row, as in FITS).
            values (list): Each source's flux per pixel area (e.g. Jy/arcsec^2).
            pixel_scale (float): Pixel size in arcsec.

        Returns:
            np.ndarray: (N, N) image of the sources. Sources beyond the radial PSFs are left out silently.
        """
        return cir_psf_planets(positions, values, pixel_scale, self.optics.img_pixel, self.radial_psfs())


# How to make the image size odd (for optics that need one; the %d are the odd sizes around it)
GRID_SIZE_FIX = "set [grid] size = %d or %d (the image is resampled onto it)"


def _vortex_optics(values, img_pixel, psf_range, wavelength, diameter, size_fix):
    """
    The Optics of a "vortex" coronagraph table (hcipy; see optics.VortexBackend).

    Args:
        values (dict): The table's values with the defaults filled in (charge, lyot_mask_size, ...).
        img_pixel (int): Number of pixels on one edge of the image.
        psf_range (float): Half-width of the image, in lambda/D.
        wavelength (float): The scene's wavelength in m (for the pixel scale in the warning).
        diameter (float): The scene's telescope diameter in m (idem).
        size_fix (str): How to make the image size odd; not needed (any size works).

    Returns:
        Optics: Optics.vortex(charge, img_pixel, psf_range, lyot_mask_size).

    Warns:
        SamplingWarning: fewer than MIN_SAMPLING pixels per lambda/D (the model samples its pupil from the
            image's pixels, so its PSFs' far wings come out too low near the image edge).
    """
    sampling = img_pixel/(2*psf_range)
    if sampling < MIN_SAMPLING:
        warnings.warn(SamplingWarning(
            "vortex (hcipy): %.3g pixels per lambda/D is too coarse for the model's PSFs, whose pupil is sampled from "
            "the image's pixels: their far wings come out a few %% too low near the image edge (5%% at 2.4 pixels "
            "per lambda/D). For about 1%%, use at least %d pixels per lambda/D: [grid] sampling = %d, or pixel_scale "
            "<= %.4g arcsec" % (sampling, MIN_SAMPLING, MIN_SAMPLING, wavelength/diameter*P.RAD_TO_ARCSEC/MIN_SAMPLING)),
            stacklevel=2)
    return Optics.vortex(values["charge"], img_pixel, psf_range, values["lyot_mask_size"])


# What "lazuli" needs beyond the package's requirements (not on PyPI in these versions), and how to install it
LAZULI_PACKAGES = {"prysm": "git+https://github.com/brandondube/prysm.git@5008f47803af6ff69eed474fca049456076e4a7d",
                   "config_um": "git+https://github.com/uasal/config_um.git",
                   "config_um_wcc": "git+https://github.com/uasal/config_um_wcc.git"}


def _lazuli_optics(values, img_pixel, psf_range, wavelength, diameter, size_fix):
    """
    The Optics of a "lazuli" coronagraph table (prysm; see toycoronagraph.lazuli, imported here because it
    needs prysm and the Lazuli config packages).

    Args:
        values (dict): The table's values with the defaults filled in (charge, wavelengths, throughput, opd, ...).
        img_pixel (int): Number of pixels on one edge of the image (odd).
        psf_range (float): Half-width of the image, in scene lambda/D.
        wavelength (float): The scene's wavelength in m.
        diameter (float): The scene's telescope diameter in m.
        size_fix (str): How to make the image size odd, in the error for an even one (see GRID_SIZE_FIX).

    Returns:
        Optics: Optics.lazuli, its "auto" settings resolved (lazuli.optics_from_table).

    Raises:
        ImportError: prysm, config_um or config_um_wcc is not installed (the message gives its pip command,
            LAZULI_PACKAGES).
        ParamError: img_pixel is even, or the pupil sampling is too coarse for the image (see
            lazuli.optics_from_table).
        ValueError: "auto" throughput without light at the wavelengths, or an OPD map that is not valid or
            does not match pupil_samples.
        OSError: the OPD map cannot be read.
    """
    try:  # prysm is imported with the module, the config packages when it reads Lazuli's parameters
        from toycoronagraph import lazuli
        return lazuli.optics_from_table(values, img_pixel, psf_range, wavelength, diameter, size_fix)
    except ModuleNotFoundError as err:
        missing = (err.name or "").split(".")[0]  # the package (err.name may be a module of it)
        if missing not in LAZULI_PACKAGES:
            raise
        raise ImportError('[[coronagraph]] type = "lazuli" needs %s, which is not installed: pip install %s (README, '
                          '"Lazuli PSFs": prysm, config_um and config_um_wcc)' % (missing, LAZULI_PACKAGES[missing])) from err


# Coronagraph type -> function(values, img_pixel, psf_range, wavelength, diameter, size_fix) returning its Optics
OPTICS_BUILDERS = {"vortex": _vortex_optics, "lazuli": _lazuli_optics}


def optics_for(ctype, values, img_pixel, psf_range, wavelength, diameter, size_fix=GRID_SIZE_FIX):
    """
    The Optics of a coronagraph of any type (used by coronagraph_for and by Target).

    Args:
        ctype (str): The coronagraph type (a key of toycoronagraph.params.CORONAGRAPHS and OPTICS_BUILDERS).
        values (dict): Its parameters, every key of its declaration present (defaults filled in).
        img_pixel (int): Number of pixels on one edge of the image.
        psf_range (float): Half-width of the image, in lambda/D of the scene.
        wavelength (float): The scene's wavelength in m.
        diameter (float): The scene's telescope diameter in m.
        size_fix (str): How the user makes the image size odd, for optics that need it (two %d for the
            odd sizes around img_pixel); the default names [grid] size.

    Returns:
        Optics: The optical system and image sampling.

    Raises:
        ParamError: The optics cannot work on this image (a lazuli coronagraph on an even size, or with
            too few pupil samples for its field).
        ValueError, OSError, ImportError: A lazuli coronagraph's settings cannot be resolved (see
            _lazuli_optics).
    """
    return OPTICS_BUILDERS[ctype](values, img_pixel, psf_range, wavelength, diameter, size_fix)


def coronagraph_for(grid, table, params):
    """
    The Coronagraph of a [[coronagraph]] table on a grid.

    Args:
        grid (Grid): The scene grid (size, field, wavelength and diameter set the Optics).
        table (dict, str or SimpleNamespace): A [[coronagraph]] table as a dict, a string
            "<type>:<charge>" (e.g. "vortex:6", "lazuli:0"; "" type: vortex, no charge: the type's
            default), or a table already checked by toycoronagraph.params.
        params (Params): The scene's parameters ([compute] and [method] settings).

    Returns:
        tuple: (the checked table: SimpleNamespace with type, name, values; the Coronagraph).

    Raises:
        ParamError: The table is not valid, or its optics cannot work on the grid (see optics_for).
        ValueError, OSError, ImportError: A lazuli coronagraph's settings cannot be resolved (see optics_for).
    """
    if isinstance(table, str):  # "vortex:6"
        kind, _, charge = table.partition(":")
        table = {"type": kind or "vortex"}
        if charge:
            table["charge"] = int(charge)
    if isinstance(table, dict):
        table = P.make({}, [], [table]).coronagraphs[0]
    values = {p.name: p.default for p in P.CORONAGRAPHS[table.type][1]}
    values.update(table.values)
    optics = optics_for(table.type, values, grid.size, grid.visual_range/2, grid.wavelength, grid.diameter)
    return table, Coronagraph.get(optics, table.name, num_cores=params.get("compute.num_cores"),
                           error_target=params.get("method.error_target"), max_kernel_gb=params.get("compute.max_kernel_gb"),
                           use_gpu=params.get("compute.use_gpu"), cache_dir=params.get("compute.psf_cache_dir"),
                           symmetric=values["symmetric"], rot_number=params.get("method.rot_number"),
                           symmetry_tol=params.get("method.symmetry_tol"))


def observe(scene, coronagraphs=None, combine=None, method=None, components=None, iwa_ignore=False, verbose=True):
    """
    The final images of a scene through one or more coronagraphs.

    Args:
        scene (Scene): Rendered here if it is not yet.
        coronagraphs (list): [[coronagraph]] tables (dicts or checked tables), or strings like
            "vortex:6"; None: the scene's parameters (a charge-2 vortex if they have none).
        combine (str): "separate" or "together"; None: [method] combine.
        method (str): "auto", "circular" or "spatial" for every component; None: each
            component's own `method`, else [method] method.
        components (list): Only these components (names, images or point sources), e.g. one
            for an individual analysis; None: all. The spatial kernel size still comes from
            the whole scene, so individual runs add up exactly to the combined one.
        iwa_ignore (bool): Leave out the light inside the inner working angle (images only).
        verbose (bool): Print the method chosen for each image (and the scene's notes if it is rendered here).

    Returns:
        Result: finals[coronagraph name][component name, or "total" for the images with
            combine = "together"], in scene.internal_unit.

    Raises:
        ParamError: An unknown combine, components that are not in the scene, or a coronagraph that
            cannot work on the scene grid (e.g. "lazuli" on an even size); before any PSF is computed.
    """
    if scene.images is None:
        scene.render(verbose=verbose)
    p, g = scene.params, scene.grid
    combine = combine or p.get("method.combine")
    if combine not in ("separate", "together"):
        raise P.ParamError('combine = %r: use "separate" or "together"' % combine)
    available = list(scene.images) + list(scene.points)
    unknown = [n for n in components or [] if n not in available]
    if unknown:
        raise P.ParamError("components %s are not in the scene; it has: %s" % (", ".join(unknown), ", ".join(available)))
    chosen = lambda name: components is None or name in components  # whether a component is in this run
    images = {n: img for n, img in scene.images.items() if chosen(n)}
    points = {n: pts for n, pts in scene.points.items() if chosen(n)}
    result = Result(g, scene.internal_unit, scene.star, scene.star_flux, combine, images, points, p)
    result.method_used = {"combine": combine, **({"method": method} if method else {})}

    def known_symmetry(name):
        """
        Whether a component is circularly symmetric, when known without measuring it.

        Args:
            name (str): The component's name in the scene.

        Returns:
            tuple or None: (bool, reason str) when the user set it or the module knows it; None: measure it.
        """
        source = scene.component_of[name].source
        sym = source.symmetric()
        if sym is None:
            return None
        set_by_user = source.settings.get("symmetric", P.AUTO) != P.AUTO
        return sym, "target %s: %s" % ("symmetric" if sym else "not symmetric",
                                       "set in the parameters" if set_by_user else "known by the %s module" % source.kind)

    # parts: {name: (image, known symmetry or None, wanted method)} going through each coronagraph
    if combine == "together":
        # the sum is symmetric if every component is known to be; not if one is known not to be; else measured
        known = [known_symmetry(n) for n in images]
        if any(k and not k[0] for k in known):
            symmetric = (False, "target not symmetric: a component is not (%s)" % "; ".join(k[1] for k in known if k and not k[0]))
        elif known and all(known):
            symmetric = (True, "target symmetric: every component is")
        else:
            symmetric = None
        parts = {"total": (sum(images.values()), symmetric, method or p.get("method.method"))} if images else {}
    else:
        own = lambda n: scene.component_of[n].method  # the component's own method setting
        parts = {n: (img, known_symmetry(n), method or (own(n) if own(n) != P.AUTO else p.get("method.method")))
                 for n, img in images.items()}
    # every coronagraph is set up (and its optics checked against the grid) before any PSF is computed
    for table, coro in [coronagraph_for(g, table, p) for table in coronagraphs or p.coronagraphs]:
        result.tables.append(table)
        finals = {}
        for name, (img, symmetric, wanted) in parts.items():
            m, mode, messages = coro.choose(img, symmetric, wanted, p.get("method.spatial_mode"))
            final, error = coro.image(img, m, mode, weight=scene.total, iwa_ignore=iwa_ignore, report=False)
            finals[name] = final
            result.methods[coro.name, name] = (m, mode, " ".join(messages), error)
            if verbose:
                print("%s, %s: %s%s" % (coro.name, name, messages[0], "; estimated error %.2g%%" % (100*error) if error else ""))
        for name, (positions, values) in points.items():
            finals[name] = coro.point_sources([g.pixel(r, theta) for r, theta in positions], values)
            result.methods[coro.name, name] = ("exact", None, "point sources: exact PSFs", 0.0)
        result.finals[coro.name] = finals
        result.coronagraphs[coro.name] = coro
        result.peak_norms[coro.name] = coro.peak_norm
    return result


def video(scene, coronagraph, component, frames, path=None, fps=10, unit="contrast_peak", log=None, log_range=None,
          axis=None, dpi=None, iwa_ignore=False, difference=False, **kwargs):
    """
    A video of a scene through one coronagraph while one component changes, e.g. a CME moving outwards.

    The scene is observed once, every component on its own with its current settings (combine =
    "separate", whatever [method] combine says). Then, per frame, the
    component gets that frame's settings (checked by its module), the scene is rendered again
    and only the component is observed again, with the method chosen the first time and the
    spatial kernel chosen for the first scene (its model comes from the caches: about one
    convolution per frame); the other components' final images are reused. The colour scale
    is fixed by the first frame. The component gets its own settings back at the end. The video
    is lossless (toycoronagraph.video). A faint component (e.g. a CME against a stellar wind) shows
    with difference=True: each frame minus the other components' final images, like the
    base-difference images (minus a pre-event image) that reveal solar CMEs against the corona.

    Args:
        scene (Scene): The scene.
        coronagraph (str or dict): One coronagraph, as in observe() (e.g. "vortex:6").
        component (str): The component that changes (an image, not a point source).
        frames (list): One dict per frame: the component's settings in it, e.g. {"horizontal_height_factor": 2.0}.
        path (str): The .webm file; None: <component>_<coronagraph name>.webm in the working directory.
        fps (int): Frames per second.
        unit (str): Output unit (see Result.in_unit; "contrast_peak": relative to the coronagraph's own peak).
        log (bool): Log colour scale. None here and below: the [display] setting.
        log_range (float): Log scale: decades shown below the first frame's top.
        axis (str): Axis unit (see Grid.axis).
        dpi (int): Dots per inch of the 6.4 x 5 inch frames.
        iwa_ignore (bool): Leave out the light inside the inner working angle.
        difference (bool): Show only the component's final image (the frame minus the others').
        **kwargs: Passed to matplotlib's imshow, e.g. cmap="magma".

    Returns:
        tuple: (path, list of the (N, N) images shown, in `unit`, one per frame: the final image of all
            components summed, or with difference=True of the component alone).
    """
    import matplotlib.pyplot as plt
    import tqdm
    from toycoronagraph.scene import show_image, display_setting
    from toycoronagraph.video import write_video
    result = observe(scene, [coronagraph], combine="separate", iwa_ignore=iwa_ignore, verbose=False)
    name = next(iter(result.finals))
    coro, weight = result.coronagraphs[name], scene.total
    method, mode = result.methods[name, component][:2]
    others = 0.0 if difference else result._sum({n: f for n, f in result.finals[name].items() if n != component})
    part = scene.component_of[component]
    original = part.source
    images = []
    try:
        for settings in tqdm.tqdm(frames, desc="%s through %s" % (component, name), mininterval=2):
            part.source = type(original)(**{**original.settings, **settings})
            final = coro.image(scene.render(verbose=False)[component], method, mode, weight, iwa_ignore, report=False)[0]
            images.append(result.in_unit(others+final, unit, coronagraph=name))
    finally:
        part.source = original
        scene.render(verbose=False)
    path = path or "%s_%s.webm" % (component, name)
    log, log_range, axis, dpi = (display_setting(scene.params, key, value) for key, value in
                                 (("log", log), ("log_range", log_range), ("axis", axis), ("dpi", dpi)))
    fig, ax = plt.subplots(figsize=(6.4, 5))
    show_image(ax, fig, images[0], scene.grid, axis, scene.star.distance, log, log_range, unit, **kwargs)
    shown = ax.images[0]
    title = ax.set_title("")

    def update(i):
        """Frame i: its final image, and its settings in the title."""
        shown.set_data(np.log10(np.where(images[i] > 0, images[i], np.nan)) if log else images[i])
        text = ", ".join("%s = %s" % (k, "%.3g" % v if isinstance(v, float) else v) for k, v in frames[i].items())
        title.set_text("%s%s through %s: %s" % (component, "" if difference else " (and the rest)", name, text))
        return [shown, title]
    write_video(fig, update, len(frames), path, fps, dpi)
    plt.close(fig)
    return path, images


def _scale_top(images, log):
    """
    The top of a colour scale shared by images (see scene.show_image): the largest of their 99.99th
    percentiles of log10 of the positive pixels (log), or of their maxima (linear).

    Args:
        images (list): (N, N) images, or None (left out).
        log (bool): Log scale.

    Returns:
        float or None: The top (log10 of the value on a log scale); None if no image has positive pixels (log)
            or there is no image.
    """
    images = [image for image in images if image is not None]
    if log:
        tops = [np.percentile(np.log10(image[image > 0]), 99.99) for image in images if np.any(image > 0)]
    else:
        tops = [np.max(image) for image in images]
    return max(tops) if tops else None


# FITS extension of the "total" component of a combine = "together" run (the final image of the summed
# images): "FINAL_" + "TOTAL" is taken by the sum of every final image, and no component's "FINAL_<NAME>"
# can be this name (see Result.save)
TOGETHER_EXTENSION = "TOGETHER"


class Result:
    """
    The final images of a scene, and what made them.

    Attributes:
        grid (Grid): The scene grid.
        unit (str): The unit of every image: the scene's internal unit ("Jy/arcsec^2" or "contrast_pixel").
        star (Star): The scene's star.
        star_flux (float or None): The star's flux density in Jy (None: unknown).
        combine (str): "separate" or "together".
        inputs (dict): {component name: (N, N) input image} (the images of this run, not the point sources).
        points (dict): {name: (positions, fluxes)}: point sources, positions a list of (radius arcsec,
            angle rad), fluxes a 1-D array of each one's flux per pixel area in `unit`.
        params (Params or None): The scene's parameters (display defaults, run folder settings).
        finals (dict): {coronagraph name: {component name (or "total" for the images with combine =
            "together"): (N, N) float64 final image in `unit`}}.
        methods (dict): {(coronagraph name, component name): (method, spatial_mode or None, why str,
            estimated typical relative error)}; method "exact" for point sources.
        coronagraphs (dict): {coronagraph name: Coronagraph} (empty after load()).
        peak_norms (dict): {coronagraph name: the peak of the star's unocculted PSF through it, as a fraction of
            the star's light per pixel (Coronagraph.peak_norm; None: unknown)}: the contrast_peak references.
        iwas (dict): {coronagraph name: its inner working angle in arcsec, or None}, filled as iwa() is asked.
        tables (list): The [[coronagraph]] tables used (for params_used.toml).
        method_used (dict): The [method] settings given to observe() (for params_used.toml).
        folder (str): The run folder, once saved.
    """

    def __init__(self, grid, unit, star, star_flux, combine, inputs, points=None, params=None):
        """
        Args:
            grid (Grid): The scene grid.
            unit (str): The images' unit ("Jy/arcsec^2" or "contrast_pixel").
            star (Star): The star.
            star_flux (float or None): The star's flux density in Jy.
            combine (str): "separate" or "together".
            inputs (dict): {component name: (N, N) input image in `unit`}.
            points (dict): {name: (positions, fluxes per pixel area)}; None: no point sources.
            params (Params): The parameters; None: no display / run defaults.
        """
        self.grid, self.unit, self.star, self.star_flux, self.combine = grid, unit, star, star_flux, combine
        self.inputs, self.points, self.params = inputs, points or {}, params
        self.finals, self.methods, self.coronagraphs, self.peak_norms, self.iwas = {}, {}, {}, {}, {}
        self.tables, self.method_used = [], {}  # the coronagraphs and [method] settings used (for params_used.toml)

    def _coronagraph(self, name):
        """
        A coronagraph's name, defaulting to the first one.

        Args:
            name (str or None): A coronagraph name; None: the first in self.finals.

        Returns:
            str: The name.
        """
        return name if name is not None else next(iter(self.finals))

    def _sum(self, images):
        """
        The sum of images on this result's grid.

        Args:
            images (dict): {name: (N, N) image in self.unit}.

        Returns:
            np.ndarray: (N, N) float64 sum, in self.unit; zero if `images` is empty (a run without
                components, or no input images).
        """
        return sum(images.values(), np.zeros((self.grid.size, self.grid.size)))

    def final(self, coronagraph=None, component=None, unit=None, star_flux=None):
        """
        The total final image, or one component's, of a coronagraph.

        Args:
            coronagraph (str): Coronagraph name; None: the first.
            component (str): Component name (or "total" with combine = "together"); None: the sum
                of all of them (images and point sources; zero for a run without components).
            unit (str): Output unit (see in_unit); None / "auto": self.unit.
            star_flux (float or Star): See in_unit.

        Returns:
            np.ndarray: (N, N) final image in `unit`.

        Raises:
            KeyError: An unknown coronagraph or component.
        """
        name = self._coronagraph(coronagraph)
        image = self._sum(self.finals[name]) if component is None else self.finals[name][component]
        return self.in_unit(image, unit, star_flux, name)

    def in_unit(self, image, unit, star_flux=None, coronagraph=None):
        """
        An image of this result in another unit (units.convert_image).

        Args:
            image (np.ndarray): (N, N) image in self.unit.
            unit (str): "Jy/arcsec^2", "Jy/pixel", "contrast_pixel" or "contrast_peak"; None or "auto": unchanged.
            star_flux (float or Star): The star's flux in Jy, or a Star (its flux at the grid's
                wavelength), to turn a contrast result into a brightness (and back); None: self.star_flux.
            coronagraph (str): contrast_peak: the coronagraph whose final image it is (relative to the peak
                of the star's PSF through it without its mask, peak_norms); None: an input image (relative
                to the ideal aperture's PSF peak, Grid.peak_norm).

        Returns:
            np.ndarray: (N, N) image in `unit` (the same array if no conversion is needed).

        Raises:
            UnitError: The conversion needs the star's flux, or the coronagraph's peak, and it is unknown.
        """
        if unit in (None, P.AUTO):
            return image
        flux = U.star_flux_of(star_flux, self.grid.wavelength) if star_flux is not None else self.star_flux
        peak = self.grid.peak_norm if coronagraph is None else self.peak_norms.get(coronagraph)
        return U.convert_image(image, self.unit, unit, self.grid.pixel_area, flux, peak)

    def contributions(self, coronagraph=None):
        """
        Each component's fractional contribution to the final image, pixel by pixel.

        Args:
            coronagraph (str): Coronagraph name; None: the first.

        Returns:
            dict: {component name: (N, N) float64 fraction of the total final image (0 where the total is 0)};
                empty for a run without components.
        """
        finals = self.finals[self._coronagraph(coronagraph)]
        total = self._sum(finals)
        with np.errstate(invalid="ignore", divide="ignore"):
            return {n: np.where(total != 0, f/total, 0.0) for n, f in finals.items()}

    def contrast(self, planet, coronagraph=None):
        """
        The light of every component inside a planet's core.

        The core is the set of pixels where the planet's final image is above half its peak;
        each component's final image is summed over it and multiplied by the pixel area.

        Args:
            planet (str): The planet's component name.
            coronagraph (str): Coronagraph name; None: the first.

        Returns:
            dict: {component name (the planet included): float}, in self.unit x arcsec^2 (Jy for
                "Jy/arcsec^2").
        """
        finals = self.finals[self._coronagraph(coronagraph)]
        core = finals[planet] > 0.5*finals[planet].max()
        return {n: float(f[core].sum()*U.flux_per_pixel(self.unit, self.grid.pixel_area)) for n, f in finals.items()}

    # ---- plots
    def _display(self, key, value):
        """
        A plot setting: the given value, else [display] <key> of the parameters (see scene.display_setting).

        Args:
            key (str): The [display] key: "unit", "log", "log_range", "axis" or "dpi".
            value: The value given to the plot method; None: not given.

        Returns:
            The setting.
        """
        from toycoronagraph.scene import display_setting
        return display_setting(self.params, key, value)

    def _rows(self, coronagraphs, input):
        """
        The rows of a figure: the input (if wanted and there are input images or point sources) and coronagraphs.

        Args:
            coronagraphs (str or list): A coronagraph name or names; None: all.
            input (bool): Include the input row.

        Returns:
            list: (row name, {component: (N, N) image in self.unit}, coronagraph name or None for the input).
        """
        names = list(self.finals) if coronagraphs is None else [coronagraphs] if isinstance(coronagraphs, str) else coronagraphs
        rows = [("input", dict(self.inputs), None)] if input and (self.inputs or self.points) else []
        return rows + [(n, self.finals[n], n) for n in names]

    def _unit_note(self, unit):
        """str: What contrast_peak is relative to, for figure titles ("" for other units)."""
        return (" (input: relative to the ideal aperture's PSF peak; finals: to each coronagraph's PSF peak "
                "without its mask)") if unit == "contrast_peak" else ""

    def plot(self, unit=None, log=None, axis=None, path=None, show=True, dpi=None, log_range=None, star_flux=None,
             coronagraphs=None, components=None, input=True, shared=True, **kwargs):
        """
        Plot the input and the final images.

        Rows: the input (input=True) and each coronagraph; columns: each component and the total
        (when there are several components or none, or combine = "together"; zero without
        components). A single column is drawn as one row of panels: components=["total"] compares
        the coronagraphs side by side. Point sources are marked on the input row (their panel there
        is blank: they have no input image), with their flux in the plot's unit family. With unit =
        "contrast_peak" each row is relative to its own peak (see in_unit). Each column has one
        colour scale (shared=True), so the rows compare directly.

        Args:
            unit (str): Output unit (see in_unit); "auto": self.unit. None here and below: the [display] setting.
            log (bool): Log colour scale.
            axis (str): Axis unit (see Grid.axis).
            path (str): Save the figure there; None: not saved.
            show (bool): Show the figure (else it is closed).
            dpi (int): Resolution of the saved file.
            log_range (float): Log scale: decades shown below the top of the scale.
            star_flux (float or Star): See in_unit.
            coronagraphs (str or list): The coronagraphs shown; None: all.
            components (list): The columns (component names, "total"); None: all, see above.
            input (bool): Show the input row.
            shared (bool): One colour scale per column, set by its brightest panel (log: the 99.99th
                percentile, see scene.show_image); False: one per panel.
            **kwargs: Passed to matplotlib's imshow for every image, e.g. cmap="magma".

        Returns:
            matplotlib.figure.Figure: The figure.
        """
        import matplotlib.pyplot as plt
        from toycoronagraph.scene import show_image, mark_points, point_labels, blank_panel, finish
        unit, log, log_range, axis = (self._display(key, value) for key, value in
                                      (("unit", unit), ("log", log), ("log_range", log_range), ("axis", axis)))
        unit = self.unit if unit == P.AUTO else unit
        names = list(self.inputs) + list(self.points)
        columns = components or names + (["total"] if len(names) != 1 or self.combine == "together" else [])
        rows = self._rows(coronagraphs, input)
        panels = []  # per row, each column's image in `unit` (None: a point source's input panel)
        for _, images, coro in rows:
            picked = [self._sum(images) if col == "total" else images.get(col) for col in columns]
            panels.append([None if image is None else self.in_unit(image, unit, star_flux, coro) for image in picked])
        tops = [_scale_top([row[j] for row in panels], log) if shared else None for j in range(len(columns))]
        flux = U.star_flux_of(star_flux, self.grid.wavelength) if star_flux is not None else self.star_flux
        labels = point_labels(self.points, self.unit, self.grid.pixel_area, unit, flux)
        shape = (1, len(rows)) if len(columns) == 1 else (len(rows), len(columns))
        fig, axes = plt.subplots(*shape, figsize=(4.1*shape[1], 3.3*shape[0]), squeeze=False)
        for i, (row, _, _) in enumerate(rows):
            for j, col in enumerate(columns):
                ax = axes.flat[i*len(columns)+j]
                if panels[i][j] is None:  # a point source in the input row: no image, its position marked
                    factor = blank_panel(ax, self.grid, axis, self.star.distance)
                else:
                    factor = show_image(ax, fig, panels[i][j], self.grid, axis, self.star.distance, log, log_range, unit,
                                        vmax=tops[j], **kwargs)
                if row == "input" and (col in self.points or col == "total"):
                    mark_points(ax, {n: v for n, v in labels.items() if col in (n, "total")}, factor)
                ax.set_title("%s: %s" % (row, col), fontsize=10)
                ax.tick_params(labelsize=8)
        fig.suptitle("%s%s, combine = %s" % (unit, self._unit_note(unit), self.combine), fontsize=10)
        fig.tight_layout()
        return finish(fig, path, show, self._display("dpi", dpi))

    def iwa(self, coronagraph):
        """
        A coronagraph's inner working angle (Coronagraph.iwa: where a source's throughput reaches half), in arcsec.

        Computed when first asked (a few seconds) and kept in iwas; save() writes it and load() reads it back.

        Args:
            coronagraph (str): The coronagraph's name.

        Returns:
            float or None: None without a coronagraph (e.g. the Lazuli telescope alone), or when it is not known (a run
                loaded from files written without it).
        """
        if coronagraph not in self.iwas:
            coro = self.coronagraphs.get(coronagraph)
            iwa = coro.iwa if coro is not None else None
            self.iwas[coronagraph] = None if iwa is None else iwa*self.grid.lambda_d
        return self.iwas[coronagraph]

    def plot_profiles(self, components=False, unit=None, log=None, axis=None, path=None, show=True, dpi=None,
                      star_flux=None, extra=None, coronagraphs=None, iwa=True, log_range=None, **kwargs):
        """
        Azimuthal mean against the distance from the star (see scene.plot_profiles).

        components=False: one panel, the summed input and each coronagraph's total final image.
        components=True: one panel per coronagraph (shared axes), each component's input (dashed)
        and final image (solid, same colour), and the totals (thick).

        Args:
            components (bool): One panel per coronagraph with each component (see above).
            unit (str): Output unit (see in_unit); "auto": self.unit. None here and below: the [display] setting.
            log (bool): Log-log axes (else linear, both from 0).
            axis (str): Distance unit (see Grid.axis).
            path (str): Save the figure there; None: not saved.
            show (bool): Show the figure (else it is closed).
            dpi (int): Resolution of the saved file.
            star_flux (float or Star): See in_unit.
            extra (dict): More curves {label: (radius in arcsec, values in `unit`)}, dashed, on every panel.
            coronagraphs (str or list): The coronagraphs shown; None: all.
            iwa (bool): Mark each coronagraph's inner working angle (see iwa) with a dotted vertical line, in its
                curve's colour (components=False) or grey (one panel per coronagraph).
            log_range (float): Log axes show log_range + scene.PROFILE_DECADES decades below the highest point.
            **kwargs: Passed to matplotlib's plot for every curve, over the plot's own styles, e.g. lw=2.

        Returns:
            matplotlib.figure.Figure: The figure.
        """
        import matplotlib.pyplot as plt
        from toycoronagraph.scene import plot_profiles, finish, PROFILE_DECADES
        unit, log, log_range, axis = (self._display(key, value) for key, value in
                                      (("unit", unit), ("log", log), ("log_range", log_range), ("axis", axis)))
        unit = self.unit if unit == P.AUTO else unit
        factor = self.grid.axis(axis, self.star.distance)[0]
        decades = log_range+PROFILE_DECADES

        def mark_iwa(ax, name, colour):
            """A dotted vertical line at the coronagraph's IWA, labelled "IWA" (none when it has no IWA)."""
            angle = self.iwa(name) if iwa else None
            if angle:
                ax.axvline(angle*factor, color=colour, ls=":", lw=1.2)
                ax.text(angle*factor, 0.98, "IWA ", transform=ax.get_xaxis_transform(), rotation=90, ha="right",
                        va="top", fontsize=7, color=colour)

        convert = lambda image, coro=None: self.in_unit(image, unit, star_flux, coro)
        rows = self._rows(coronagraphs, False)
        if not components:
            curves = {"input": convert(self._sum(self.inputs))} if self.inputs else {}
            curves.update({name: convert(self._sum(finals), name) for name, finals, _ in rows})
            fig, ax = plt.subplots(figsize=(6.4, 4.4))
            plot_profiles(ax, curves, self.grid, axis, self.star.distance, log, unit, extra, log_range=decades, **kwargs)
            colours = {line.get_label(): line.get_color() for line in ax.lines}
            for name, _, _ in rows:
                mark_iwa(ax, name, colours[name])
            ax.set_title("Total of each coronagraph%s" % self._unit_note(unit), fontsize=9, wrap=True)
            panels = [ax]
        else:
            columns = min(3, len(rows))
            shape = (-(-len(rows)//columns), columns)
            fig, axes = plt.subplots(*shape, figsize=(5.2*shape[1], 4.0*shape[0]), squeeze=False, sharex=True, sharey=True)
            colours = {n: "C%d" % (i % 10) for i, n in enumerate(list(self.inputs) + list(self.points))}
            for ax, (name, finals, _) in zip(axes.flat, rows):
                curves, styles = {}, {}
                for n, image in self.inputs.items():
                    curves["%s (input)" % n] = convert(image)
                    styles["%s (input)" % n] = dict(color=colours[n], ls="--", lw=1)
                for n, image in finals.items():
                    curves[n] = convert(image, name)
                    styles[n] = dict(color=colours.get(n, "k"), lw=1.5)
                if len(finals) > 1:
                    curves["total"] = convert(self._sum(finals), name)
                    styles["total"] = dict(color="k", lw=4, alpha=0.3, zorder=1)  # behind, so no curve is hidden
                plot_profiles(ax, curves, self.grid, axis, self.star.distance, log, unit, extra, styles,
                              log_range=decades, legend=ax is axes.flat[0], **kwargs)
                mark_iwa(ax, name, "0.35")
                ax.set_title(name, fontsize=10)
            for ax in axes.flat[len(rows):]:
                ax.set_visible(False)
            fig.suptitle("Each component before (dashed) and after (solid) the coronagraph%s" % self._unit_note(unit),
                         fontsize=10)
            panels = axes.flat[:len(rows)]
        for ax in panels:
            ax.tick_params(labelsize=8)
        fig.tight_layout()
        return finish(fig, path, show, self._display("dpi", dpi))

    def plot_contributions(self, coronagraphs=None, input=True, axis=None, path=None, show=True, dpi=None, **kwargs):
        """
        Each component's share of the image, pixel by pixel: before and after each coronagraph.

        Rows: the input (input=True) and each coronagraph; columns: each component's fraction of
        that row's total image (see contributions; one colour scale, 0 to 1), titled with its share
        of the row's light, and a last panel colouring each pixel by the component that dominates
        it. Point sources have no input image: their input panel is blank, with their position.

        Args:
            coronagraphs (str or list): The coronagraphs shown; None: all.
            input (bool): Show the input row.
            axis (str): Axis unit (see Grid.axis); None: [display] axis.
            path (str): Save the figure there; None: not saved.
            show (bool): Show the figure (else it is closed).
            dpi (int): Resolution of the saved file; None: [display] dpi.
            **kwargs: Passed to matplotlib's imshow for every share map, e.g. cmap="magma".

        Returns:
            matplotlib.figure.Figure: The figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Patch
        from toycoronagraph.scene import blank_panel, mark_points, finish
        axis = self._display("axis", axis)
        factor, label = self.grid.axis(axis, self.star.distance)
        extent = [e*factor for e in self.grid.extent]
        rows = self._rows(coronagraphs, input)
        columns = list(dict.fromkeys(n for _, images, _ in rows for n in images))
        columns += [n for n in self.points if n not in columns]  # point sources (no input image)
        colours = ["C%d" % (i % 10) for i in range(len(columns))]
        fig, axes = plt.subplots(len(rows), len(columns)+1, figsize=(3.2*(len(columns)+1)+1.0, 2.9*len(rows)),
                                 squeeze=False, layout="constrained")
        im = None
        for i, (row, images, _) in enumerate(rows):
            total = self._sum(images)
            light = {n: float(image.sum()) for n, image in images.items()}
            if row == "input":  # point sources hold their whole flux in one pixel's worth of flux per pixel area
                light.update({n: float(np.sum(values)) for n, (_, values) in self.points.items()})
            everything = sum(light.values())
            with np.errstate(invalid="ignore", divide="ignore"):
                fractions = {n: np.where(total != 0, image/total, 0.0) for n, image in images.items()}
            for j, n in enumerate(columns):
                ax = axes[i][j]
                if n in fractions:
                    im = ax.imshow(fractions[n], origin="lower", extent=extent, **{"cmap": "viridis", **kwargs}, vmin=0,
                                   vmax=1)
                else:
                    blank_panel(ax, self.grid, axis, self.star.distance)
                    if n in self.points:
                        mark_points(ax, {n: (self.points[n][0], [n])}, factor)
                share = "%.3g%%" % (100*light[n]/everything) if n in light and everything else "-"
                ax.set_title("%s: %s, %s" % (row, n, share), fontsize=9)
            ax = axes[i][-1]
            if fractions:
                stack = np.array([fractions.get(n, np.zeros_like(total)) for n in columns])
                dominant = np.where(total > 0, np.argmax(stack, axis=0), np.nan)
                ax.imshow(dominant, origin="lower", extent=extent, cmap=ListedColormap(colours), vmin=-0.5,
                          vmax=len(columns)-0.5, interpolation="nearest")
            else:
                blank_panel(ax, self.grid, axis, self.star.distance)
            ax.set_title("%s: dominant component" % row, fontsize=9)
        for ax in axes.flat:
            ax.set(xlabel="x [%s]" % label, ylabel="y [%s]" % label)
            ax.tick_params(labelsize=7)
        if im is not None:
            fig.colorbar(im, ax=axes[:, :-1], location="bottom", label="fraction of the row's image", shrink=0.6,
                         aspect=50)
        axes[0][-1].legend(handles=[Patch(color=c, label=n) for c, n in zip(colours, columns)], fontsize=7,
                           loc="upper right", framealpha=0.8)
        fig.suptitle("Each component's share of the image, before and after the coronagraphs (titles: share of the "
                     "row's light)", fontsize=10)
        return finish(fig, path, show, self._display("dpi", dpi))

    def plot_compare(self, first, second, component=None, unit=None, log=None, axis=None, path=None, show=True,
                     dpi=None, log_range=None, star_flux=None, **kwargs):
        """
        Compare two coronagraphs (two sets of PSFs) on the same image.

        Panels: each one's final image (one colour scale); their relative difference, (second - first)
        / first, where the first is within log_range decades of its peak (elsewhere grey), titled with
        its median and 99th percentile in absolute value; and the azimuthal means of the input and of
        both final images.

        Args:
            first, second (str): Coronagraph names.
            component (str): Component name; None: the total.
            unit (str): Output unit (see in_unit); "auto": self.unit. None here and below: the [display] setting.
            log (bool): Log colour scale and log-log profiles.
            axis (str): Axis / distance unit (see Grid.axis).
            path (str): Save the figure there; None: not saved.
            show (bool): Show the figure (else it is closed).
            dpi (int): Resolution of the saved file.
            log_range (float): Decades shown below the top of the colour scale, and kept in the difference.
            star_flux (float or Star): See in_unit.
            **kwargs: Passed to matplotlib's imshow for the two final images, e.g. cmap="magma".

        Returns:
            matplotlib.figure.Figure: The figure.
        """
        import matplotlib
        import matplotlib.pyplot as plt
        from toycoronagraph.scene import show_image, plot_profiles, finish, PROFILE_DECADES
        unit, log, log_range, axis = (self._display(key, value) for key, value in
                                      (("unit", unit), ("log", log), ("log_range", log_range), ("axis", axis)))
        unit = self.unit if unit == P.AUTO else unit
        what = component or "total"
        a, b = (self.final(name, component, unit, star_flux) for name in (first, second))
        fig, axes = plt.subplots(1, 4, figsize=(18, 3.9))
        top = _scale_top([a, b], log)
        for ax, name, image in zip(axes, (first, second), (a, b)):
            show_image(ax, fig, image, self.grid, axis, self.star.distance, log, log_range, unit, vmax=top, **kwargs)
            ax.set_title("%s: %s" % (name, what), fontsize=10)
        kept = a > a.max()*10.0**-log_range
        with np.errstate(invalid="ignore", divide="ignore"):
            difference = np.where(kept, (b-a)/a, np.nan)
        size = np.abs(difference[kept])
        limit = np.percentile(size, 99) if size.size else 1.0
        factor, label = self.grid.axis(axis, self.star.distance)
        cmap = matplotlib.colormaps["RdBu_r"].with_extremes(bad="0.85")
        im = axes[2].imshow(difference, origin="lower", extent=[e*factor for e in self.grid.extent], cmap=cmap,
                            vmin=-limit, vmax=limit)
        fig.colorbar(im, ax=axes[2], fraction=0.046, label="(%s - %s) / %s" % (second, first, first))
        axes[2].set(xlabel="x [%s]" % label, ylabel="y [%s]" % label)
        axes[2].set_title("relative difference: median %.2g%%, 99%% within %.2g%%"
                          % (100*np.median(size) if size.size else 0, 100*limit), fontsize=10)
        source = self._sum(self.inputs) if component is None else self.inputs.get(component)
        curves = {} if source is None else {"input": self.in_unit(source, unit, star_flux)}
        curves.update({first: a, second: b})
        plot_profiles(axes[3], curves, self.grid, axis, self.star.distance, log, unit,
                      styles={"input": dict(color="k", ls="--", lw=1.2), second: dict(ls=":", lw=2)},
                      log_range=log_range+PROFILE_DECADES)
        axes[3].set_title("azimuthal means", fontsize=10)
        fig.suptitle("Two sets of PSFs on the same image: %s and %s (%s)%s" % (first, second, what, self._unit_note(unit)),
                     fontsize=10)
        fig.tight_layout()
        return finish(fig, path, show, self._display("dpi", dpi))

    # ---- saving
    def save(self, folder=None, overwrite=None):
        """
        Save to a run folder (needs self.params).

        Files written in the folder:
            params_used.toml: the parameters used (Params.used_toml, with the coronagraphs and [method] settings).
            scene.fits: the inputs (see scene.write_scene; primary header also STARFLUX in Jy, 0 if
                unknown, and PARHASH), and a POINTS binary table of the point sources (NAME, RADIUS
                arcsec, ANGLE rad, FLUX per pixel area in self.unit) if there are any.
            <coronagraph>.fits: one per coronagraph; primary header CORONAG, OPTICS and the optics
                parameters (keywords cut to 8 characters, lists such as the wavelengths as comma-separated
                text), COMBINE, BUNIT, STARFLUX, PEAKNORM (peak_norms; 0 if unknown), IWA (arcsec, iwa(); 0: none),
                PARHASH (first 16 hex characters of the SHA-1 of params_used.toml); then FINAL_TOTAL (the sum of the others) and FINAL_<COMPONENT> images
                (with METHOD, MODE, ERRTYP and WHY, cut to 68 characters); with combine = "together", the
                final image of the summed images (component "total") is TOGETHER (TOGETHER_EXTENSION).
            quicklook/finals.png (see plot), and quicklook/contributions.png when there are several
                components (see plot_contributions).

        Args:
            folder (str): The run folder; None: [run] output_dir / [run] name ("auto": the parameter
                file's name, else "scene").
            overwrite (bool): Write into an existing folder; else the new one gets _2, _3, ...;
                None: [run] overwrite.

        Returns:
            str: The folder written (also self.folder).
        """
        from toycoronagraph.scene import write_scene, image_hdu
        p = self.params
        if folder is None:
            name = p.get("run.name")
            if name == P.AUTO:
                name = os.path.splitext(os.path.basename(p.path))[0] if p.path else "scene"
            folder = os.path.join(p.get("run.output_dir"), name)
        overwrite = p.get("run.overwrite") if overwrite is None else overwrite
        base, i = folder, 1
        while os.path.exists(folder) and not overwrite:
            i += 1
            folder = "%s_%d" % (base, i)
        os.makedirs(os.path.join(folder, "quicklook"), exist_ok=True)
        used = p.used_toml(self.tables or None, self.method_used)
        with open(os.path.join(folder, "params_used.toml"), "w") as f:
            f.write(used)
        digest = hashlib.sha1(used.encode()).hexdigest()[:16]
        scene_path = os.path.join(folder, "scene.fits")
        write_scene(scene_path, self.inputs, self.unit, self.grid, self.star, overwrite=True,
                    extra={"STARFLUX": (self.star_flux or 0.0, "Jy (0: unknown)"), "PARHASH": digest})
        if self.points:
            # one table row per source: (name, radius arcsec, angle rad, flux: Jy, or a fraction of the star's)
            per_source = U.flux_per_pixel(self.unit, self.grid.pixel_area)
            rows = [(n, r, theta, f*per_source) for n, (pos, fl) in self.points.items()
                    for (r, theta), f in zip(pos, np.atleast_1d(fl))]
            table = fits.BinTableHDU.from_columns([
                fits.Column(name="NAME", format="32A", array=[r[0] for r in rows]),
                fits.Column(name="RADIUS", format="D", unit="arcsec", array=[r[1] for r in rows]),
                fits.Column(name="ANGLE", format="D", unit="rad", array=[r[2] for r in rows]),
                fits.Column(name="FLUX", format="D", unit="Jy" if self.unit == "Jy/arcsec^2" else "fraction of the star",
                            array=[r[3] for r in rows])], name="POINTS")
            fits.append(scene_path, table.data, table.header)
        for coro, finals in self.finals.items():
            primary = fits.PrimaryHDU()
            h = primary.header
            h["CORONAG"] = coro
            c = self.coronagraphs.get(coro)
            if c is not None:
                h["OPTICS"] = c.optics.kind
                for k, v in c.optics.params:  # FITS keywords have at most 8 characters; lists as text
                    h[k.upper()[:8]] = ", ".join("%.6g" % x for x in v) if isinstance(v, tuple) else v
            h["COMBINE"] = self.combine
            h["BUNIT"] = self.unit
            h["STARFLUX"] = (self.star_flux or 0.0, "Jy (0: unknown)")
            h["PEAKNORM"] = (self.peak_norms.get(coro) or 0.0, "unocculted star PSF peak/pixel (0: unknown)")
            h["IWA"] = (self.iwa(coro) or 0.0, "inner working angle, arcsec (0: none/unknown)")
            h["PARHASH"] = (digest, "sha1 of params_used.toml")
            hdus = [primary, image_hdu("FINAL_TOTAL", self._sum(finals), self.unit, self.grid)]
            for name, image in finals.items():
                together = self.combine == "together" and name == "total"
                hdu = image_hdu(TOGETHER_EXTENSION if together else "FINAL_" + name.upper(), image, self.unit, self.grid)
                m, mode, why, error = self.methods[coro, name]
                hdu.header["METHOD"] = m
                hdu.header["MODE"] = mode or ""
                hdu.header["ERRTYP"] = (error, "estimated typical relative error")
                hdu.header["WHY"] = why[:68]  # what fits on one header card
                hdus.append(hdu)
            fits.HDUList(hdus).writeto(os.path.join(folder, "%s.fits" % coro), overwrite=True)
        self.plot(path=os.path.join(folder, "quicklook", "finals.png"), show=False)
        if len(self.inputs) + len(self.points) > 1:
            self.plot_contributions(path=os.path.join(folder, "quicklook", "contributions.png"), show=False)
        self.folder = folder
        return folder

    @classmethod
    def load(cls, folder):
        """
        A saved result (for plots, contributions, other units), without computing anything again.

        Reads scene.fits, params_used.toml and every other .fits file of the folder (one per
        coronagraph). Component names come back in lower case (from the FITS extension names; TOGETHER
        of a combine = "together" run is "total"); `coronagraphs` stays empty and `tables` / `method_used`
        are not restored.

        Args:
            folder (str): A run folder written by save().

        Returns:
            Result or None: The result; None if the folder has no coronagraph FITS file.
        """
        from toycoronagraph.scene import Grid, Star
        points = {}
        with fits.open(os.path.join(folder, "scene.fits")) as h:
            p0 = h[0].header
            images = [hdu for hdu in h[1:] if hdu.is_image]
            size = p0.get("SCNSIZE") or (int(images[0].header["NAXIS1"]) if images else 0)
            grid = Grid(int(size), p0["PIXSCALE"], p0["WAVELEN"], p0["DIAMETER"])
            inputs = {hdu.name.lower(): hdu.data.astype(float) for hdu in images}
            star = Star(**{k: p0.get(k.upper()[:8]) for k in ("distance", "radius", "temperature", "flux_jy")})
            star_flux = p0.get("STARFLUX") or None  # 0 was written for "unknown"
            if "POINTS" in h:
                for row in h["POINTS"].data:
                    pos, fl = points.setdefault(row["NAME"].strip(), ([], []))
                    pos.append((row["RADIUS"], row["ANGLE"]))
                    fl.append(row["FLUX"])  # per source; turned back into per pixel area below
        params = P.load(os.path.join(folder, "params_used.toml"))
        result = None
        for name in sorted(os.listdir(folder)):
            if not name.endswith(".fits") or name == "scene.fits":
                continue
            with fits.open(os.path.join(folder, name)) as h:
                if result is None:
                    unit = h[0].header["BUNIT"]
                    per_source = U.flux_per_pixel(unit, grid.pixel_area)
                    result = cls(grid, unit, star, star_flux, h[0].header["COMBINE"], inputs,
                                 {n: (pos, np.array(fl)/per_source) for n, (pos, fl) in points.items()}, params)
                coro = h[0].header["CORONAG"]
                result.peak_norms[coro] = h[0].header.get("PEAKNORM") or None  # 0 was written for "unknown"
                result.iwas[coro] = h[0].header.get("IWA") or None
                finals = result.finals[coro] = {}
                for hdu in h[1:]:
                    if hdu.name == "FINAL_TOTAL":  # the sum of the others
                        continue
                    component = "total" if hdu.name == TOGETHER_EXTENSION else hdu.name[6:].lower()  # [6:] strips "FINAL_"
                    finals[component] = hdu.data.astype(float)
                    result.methods[coro, component] = (hdu.header["METHOD"], hdu.header["MODE"] or None,
                                                       hdu.header["WHY"], hdu.header["ERRTYP"])
        return result
