"""
Coronagraphic PSFs and the two ways the package images a target through them.

Circular method (cir_*): for rotationally symmetric optics and targets. A radial series of PSFs
of sources along the +x axis (radial_psf_series, one per pixel of radius, out to the image's
edge; beyond it, out to the corners, the last one shifted) is weighted by the target's radial
profile, summed, and rotated about the star.

Spatial method (spatial_*): for any target. PSF kernels are sampled on a grid whose density
follows how fast the PSF shape changes (spatial.equal_error_grid) and blended per pixel by
spatial.spatial_psf. Kernels come either from direct propagation of every sample ('direct', any
optics) or from x-axis PSFs interpolated in source position and rotated ('rotate', rotationally
symmetric band-limited optics). spatial_error_report checks the model against exact PSFs.

Conventions (see toycoronagraph.optics): the star is at pixel coordinate c = (img_pixel-1)/2 on
both axes; a source at pixel (row, col) is at x = (col-c)*pixel_scale, y = (row-c)*pixel_scale in
lambda/D, with pixel_scale = 2*psf_range/img_pixel lambda/D per pixel. Intensities are in the
package's units (hcipy intensity / normal_factor), the same for frames, kernels and strips.
"""
import os
import numpy as np
from joblib import Parallel, delayed
from scipy.ndimage import map_coordinates, spline_filter
from skimage.transform import rotate
from toycoronagraph.tool import frame_interpolation, arcsec_to_pixel
from toycoronagraph.optics import Optics
from toycoronagraph.spatial import (equal_error_grid, blend_map, positions_to_polar, bandlimited_weights,
                                    spatial_psf)
import matplotlib.pyplot as plt

def radial_source_positions(img_pixel):
    """
    Source positions of the radial PSF series along the +x axis.

    For odd img_pixel the star is on a pixel and the sources are at 0, 1, 2, ... pixels. For even
    img_pixel the star lies between the four central pixels and source i is at sqrt(i^2+0.25)
    pixels (these radii lie between pixel centres).

    Args:
        img_pixel (int): Number of pixels on one edge of the image.

    Returns:
        np.ndarray: ((img_pixel+1)//2,) float64, distance of each radial PSF's source from the star, in pixels.
    """
    i = np.arange((img_pixel+1)//2)
    if img_pixel%2==0:
        return (i**2+0.25)**0.5
    return i.astype(float)

def _source_at_pixel(optics, row, col):
    """
    A point source at pixel coordinates (row, col) of the image.

    Args:
        optics (Optics): The optical system and image sampling.
        row (float): Row of the source, in pixels (need not be an integer).
        col (float): Column of the source, in pixels (need not be an integer).

    Returns:
        tuple: the backend's source object (with .throughput, .total, .block(center, shape) and
            .frame(), e.g. optics.PupilPlaneSource), and its position (x, y) in lambda/D.
    """
    x, y = (col-optics.centre)*optics.pixel_scale, (row-optics.centre)*optics.pixel_scale
    return optics.backend().source(x, y), (x, y)

def _radial_psf_worker(i, optics, frame_size=None):
    """
    Parallel worker: the i-th radial PSF, on the image or on a larger frame centred on the star.

    Args:
        i (int): Index of the source in radial_source_positions(optics.img_pixel).
        optics (Optics): The optical system and image sampling.
        frame_size (int or None): None for the image frame; otherwise the edge, in pixels, of a
            larger frame on the image pixel lattice centred on the star.

    Returns:
        np.ndarray: (img_pixel, img_pixel) float64 if frame_size is None, else
            (frame_size, frame_size) float32; intensity in the package's units.
    """
    src = optics.backend().source(radial_source_positions(optics.img_pixel)[i]*optics.pixel_scale, 0)
    if frame_size is None:
        return src.frame()
    return src.block((0.0, 0.0), (frame_size, frame_size)).astype(np.float32)

def extended_frame_size(img_pixel):
    """
    Size of radial PSF frames that, rotated about the star, still cover the image corners.

    Args:
        img_pixel (int): Number of pixels on one edge of the image.

    Returns:
        int: ceil(sqrt(2)*img_pixel)+8 pixels, plus 1 if needed so it has the parity of img_pixel
            (the frame then shares the image's pixel lattice and centre).
    """
    size = int(np.ceil(np.sqrt(2)*img_pixel))+8
    return size+(size-img_pixel)%2

def radial_psf_series(optics, num_cores, frame_size=None):
    """
    PSFs of sources along the +x axis (at radial_source_positions).

    Args:
        optics (Optics): The optical system and image sampling.
        num_cores (int): Number of CPU cores to use for parallel calculation (joblib n_jobs).
        frame_size (int): None for frames of the image size (float64); otherwise the
            size (same parity as the image) of larger float32 frames centred on the star,
            e.g. extended_frame_size(img_pixel) so rotated PSFs cover the image corners.

    Returns:
        np.ndarray: ((img_pixel+1)//2, size, size) array of PSFs (size = img_pixel or frame_size),
            in the package's intensity units.

    Raises:
        ValueError: If frame_size and img_pixel differ in parity.
    """
    if frame_size is not None and (frame_size-optics.img_pixel)%2:
        raise ValueError("frame_size must have the parity of img_pixel (same pixel lattice)")
    # joblib's workers (unlike multiprocessing on macOS/Windows) do not re-run the
    # calling script, so this also works in scripts without an `if __name__ == "__main__":` guard.
    return np.array(Parallel(n_jobs=num_cores)(
        delayed(_radial_psf_worker)(i, optics, frame_size) for i in range((optics.img_pixel+1)//2)))

def _axis_worker(r, strip, optics):
    """
    Parallel worker: throughput of a source on the +x axis and, optionally, its intensity on a block.

    Args:
        r (float): Distance of the source from the star along +x, in pixels.
        strip (tuple or None): (row_lo, row_hi, col_lo, col_hi), inclusive pixel bounds of a block
            of the image pixel lattice (which may extend beyond the image); None for the throughput only.
        optics (Optics): The optical system and image sampling.

    Returns:
        tuple: throughput (float), and the block's intensity as a (row_hi-row_lo+1, col_hi-col_lo+1)
            float64 array (package units), or None if strip is None.
    """
    src = optics.backend().source(r*optics.pixel_scale, 0)
    if strip is None:
        return src.throughput, None
    row_lo, row_hi, col_lo, col_hi = strip
    c, s = optics.centre, optics.pixel_scale
    center = (((col_lo+col_hi)/2-c)*s, ((row_lo+row_hi)/2-c)*s)
    return src.throughput, src.block(center, (row_hi-row_lo+1, col_hi-col_lo+1))

def coronagraph_throughput(optics, num_cores, radii):
    """
    Throughput of point sources on the +x axis: their total light relative to
    the same sources without the starlight-suppression element (for the
    vortex: the power through the Lyot stop with and without the vortex; by
    Parseval this is the whole focal plane's light, independent of the image frame).

    Args:
        optics (Optics): The optical system and image sampling.
        num_cores (int): Number of CPU cores to use for parallel calculation.
        radii (array-like): Source distances from the star, in pixels.

    Returns:
        np.ndarray: (len(radii),) float64, throughput (0-1) at each radius.
    """
    return np.array([t for t, _ in Parallel(n_jobs=num_cores)(
        delayed(_axis_worker)(r, None, optics) for r in radii)])

def _spatial_psf_worker(row, col, kernel_size, optics):
    """
    Parallel worker: throughput and kernel of a source at pixel (row, col).

    The kernel is propagated directly onto a kernel_size x kernel_size grid
    centred on the source (on the image pixel lattice when the source is on a
    pixel), so it is never truncated by the image frame.

    Args:
        row (float): Row of the source, in pixels (need not be an integer).
        col (float): Column of the source, in pixels (need not be an integer).
        kernel_size (int): Kernel width in pixels (odd for a kernel centred on a pixel).
        optics (Optics): The optical system and image sampling.

    Returns:
        tuple: throughput (float), and the kernel as a (kernel_size, kernel_size) float64 array
            (package units, not divided by the throughput).
    """
    src, center = _source_at_pixel(optics, row, col)
    return src.throughput, src.block(center, (kernel_size, kernel_size))

def _validation_worker(row, col, kernel_size, optics, full_frame=False):
    """
    Parallel worker: exact PSF quantities of a source at (row, col), for spatial_error_report.

    Args:
        row (float): Row of the source, in pixels.
        col (float): Column of the source, in pixels.
        kernel_size (int): Kernel width in pixels (odd).
        optics (Optics): The optical system and image sampling.
        full_frame (bool): Also return the PSF on the whole image.

    Returns:
        tuple: throughput (float), the (kernel_size, kernel_size) float64 PSF centred on the source,
            its total light over the whole focal plane (float, same units as the kernels) and, with
            full_frame, its (img_pixel, img_pixel) float64 PSF on the image.
    """
    src, center = _source_at_pixel(optics, row, col)
    kernel = src.block(center, (kernel_size, kernel_size))
    if not full_frame:
        return src.throughput, kernel, src.total
    return src.throughput, kernel, src.total, src.frame()

def _probe_worker(row, col, optics):
    """
    Parallel worker: the PSF, on the whole image, of a source at pixel (row, col).

    Args:
        row (float): Row of the source, in pixels.
        col (float): Column of the source, in pixels.
        optics (Optics): The optical system and image sampling.

    Returns:
        np.ndarray: (img_pixel, img_pixel) float64 PSF in the package's units.
    """
    return _source_at_pixel(optics, row, col)[0].frame()

# Settings of the spatial PSF grid and of mode='rotate'
_BAND_MARGIN = 1.05  # safety factor on the backend's band limit (intensity vs. source position)
_N_TAPS = 16         # nodes used per band-limited interpolation
_CIR_PIXEL_NODES = 6 # circular method: radial PSFs nearest the star weighted by the pixels' own light
_MARGIN = 16         # pixels of window beyond the rotated kernel, so the spline's boundary handling stays outside it
_MAX_SPACING = 1.0   # largest ring spacing / arc between samples, in lambda/D
_N_THETA_MIN = 8     # fewest samples on a ring
_R_INNER_MAX = 3.0   # cap on the radius of per-pixel samples, in lambda/D
_DPHI = np.radians(2) # angle step used to measure the shape change per radian
_AXIS_CACHE_VERSION = 3  # stored in axis_products' cache; bump when the strip layout changes
# Blending error vs. grid tolerance, measured for the vortex at 16 px per lambda/D:
# the largest blended-PSF shape error is ~ _BLEND_COEF * tol**_BLEND_POWER
_BLEND_COEF, _BLEND_POWER = 1.5, 1.6

def tol_for_error_target(error_target, blend_share=0.2):
    """
    Grid tolerance whose blending error fits blend_share of the image error target.

    Uses the measured law (largest blended-PSF shape error ~ 1.5 tol^1.6, vortex
    at 16 px per lambda/D); the error report checks the result for any optics.
    The kernel crop gets the other 80% of the target (see choose_kernel_size).

    Args:
        error_target (float): Target relative error of the final image.
        blend_share (float): Share of the target given to blending.

    Returns:
        float: tol (L1 shape change of unit-sum kernels between neighbouring samples), clipped to [0.003, 0.05].
    """
    return float(np.clip((blend_share*error_target/_BLEND_COEF)**(1/_BLEND_POWER), 0.003, 0.05))

def _band(optics):
    """
    Band limit of the intensity vs. source position, with the _BAND_MARGIN safety factor.

    Args:
        optics (Optics): Optics whose backend has a bandlimit (cycles per lambda/D).

    Returns:
        float: The band limit in cycles per pixel of source displacement.
    """
    return _BAND_MARGIN*optics.backend().bandlimit*optics.pixel_scale

def _rotate_geometry(img_pixel, kernel_size):
    """
    Geometry of mode='rotate'.

    Args:
        img_pixel (int): Number of pixels on one edge of the image.
        kernel_size (int): Kernel width in pixels (odd).

    Returns:
        tuple: c (float, star position in pixels on both axes), h (int, kernel half-width),
            wh (int, half-width of the rotation window: fits the rotated kernel plus _MARGIN), and
            r_near (float, pixels: the largest source radius whose window fits in image-size radial frames).
    """
    c, h = (img_pixel-1)/2, kernel_size//2
    wh = int(np.ceil(h*np.sqrt(2)))+_MARGIN+2
    # The window's last column b0+wh, with b0 = floor(c+r+0.5), must stay within img_pixel-1
    r_near = img_pixel-1.5-c-wh
    return c, h, wh, r_near

def _same_settings(cached, settings):
    """
    Whether a cached result was computed with the given settings.

    Args:
        cached (dict): Loaded cache (see load_spatial_psfs).
        settings (dict): Expected settings; strings are compared exactly, numbers and tuples of numbers
            (e.g. wavelengths in m) to a relative 1e-5.

    Returns:
        bool: True if every key of settings is in cached with a matching value.
    """
    for k, v in settings.items():
        if k not in cached:
            return False
        if isinstance(v, str):
            if str(cached[k]) != v:
                return False
        elif np.shape(cached[k]) != np.shape(v) or not np.allclose(cached[k], v, rtol=1e-5, atol=0):
            return False
    return True

def axis_products(optics, num_cores, kernel_size, cache_path=None):
    """
    Runs for sources on the +x axis, shared by every spatial PSF grid of
    rotationally symmetric, band-limited optics.

    Sources are spaced at half the Nyquist spacing of the intensity's band limit
    (so band-limited interpolation between them is exact), from the star to
    beyond the frame corners. Each run gives the throughput; runs beyond the
    reach of the radial frames also give a "far strip": the intensity on the
    block of the image pixel lattice (possibly outside the image) that
    interpolation windows around nearby radii need. These depend only on the
    optics, the image sampling and the kernel size, not on the grid tolerance,
    so they can be cached.

    Args:
        optics (Optics): The optical system and image sampling.
        num_cores (int): Number of CPU cores to use for parallel calculation.
        kernel_size (int): Kernel width in pixels (odd).
        cache_path (str): If given, load from / save to this .npz (reused when the settings match,
            otherwise recomputed and overwritten).

    Returns:
        dict: axis_x ((n,) float64, source positions in pixels from the star), throughput ((n,) float64),
            far_x ((m,) float64, the positions with a far strip), far_strips ((m, 2wh+3, width) float64,
            intensity on rows far_row0..far_row0+2wh+2 and columns far_col0[k]..far_col0[k]+width-1 of
            the image pixel lattice), far_col0 ((m,) int), far_row0 (int), and the settings
            (optics.settings(), kernel_size, cache_version); wh is from _rotate_geometry.
    """
    settings = dict(optics.settings(), kernel_size=kernel_size,
                    cache_version=_AXIS_CACHE_VERSION) # bump when the strip layout changes
    if cache_path is not None and os.path.exists(cache_path):
        cached = load_spatial_psfs(cache_path)
        if _same_settings(cached, settings):
            return cached

    c, _, wh, r_near = _rotate_geometry(optics.img_pixel, kernel_size)
    pixel_scale = optics.pixel_scale
    step = 1/(4*_band(optics))          # half the Nyquist spacing, pixels
    reach = (_N_TAPS//2+2)*step         # nodes used by a target lie within this distance
    r_top = np.sqrt(2)*c+_MAX_SPACING/pixel_scale # outermost ring radius
    axis_x = np.arange(0, r_top+reach+2*step, step)
    far_node = axis_x >= r_near-reach-step
    # One extra row on each side: sources at -x use the strips mirrored through the star,
    # which moves the window rows by a pixel for even img_pixel
    far_row0 = int(np.floor(c+0.5))-wh-1
    width = int(np.ceil(2*reach))+2*wh+6
    far_col0 = (np.floor(c+axis_x[far_node]-reach+0.5)-wh-2).astype(int)
    strips = iter([(far_row0, far_row0+2*wh+2, col0, col0+width-1) for col0 in far_col0])
    results = Parallel(n_jobs=num_cores)(
        delayed(_axis_worker)(x, next(strips) if is_far else None, optics)
        for x, is_far in zip(axis_x, far_node)
        )
    far_strips = np.array([s for (_, s), is_far in zip(results, far_node) if is_far]).reshape(-1, 2*wh+3, width)
    data = dict(axis_x=axis_x, throughput=np.array([t for t, _ in results]), far_x=axis_x[far_node],
                far_strips=far_strips, far_col0=far_col0, far_row0=far_row0, **settings)
    if cache_path is not None:
        np.savez(cache_path, **data)
    return data

def _window_builder(radial_psfs, axis, optics, kernel_size):
    """
    Build window(r): the intensity of a source at (r, 0) on a (2wh+1)^2 block
    of the pixel lattice around it, and the source's (row, col) in that block.

    For a fixed detector pixel the intensity is band-limited in the source
    position (at most the backend's bandlimit; for pupil-plane systems D/lambda
    cycles per lambda/D, because the source enters only through a tilt across
    the pupil). So the window is interpolated from x-axis PSFs at other source
    positions with band-limited weights (bandlimited_weights). Within r_near
    these are the radial frames, plus the frames flipped on both axes for
    sources at -x (180-degree symmetry), so r = 0 is interpolated from both
    sides; beyond r_near they are the far strips of axis_products. When the far
    strips start at the star (large kernels relative to the image), sources at
    -x are again included, as mirrored strips.

    Args:
        radial_psfs (np.ndarray): ((img_pixel+1)//2, img_pixel, img_pixel) radial PSFs on the image
            frame (crop extended frames first, see _crop_centre).
        axis (dict): From axis_products (far strips).
        optics (Optics): The optical system and image sampling (rotationally symmetric, band-limited).
        kernel_size (int): Kernel width in pixels (odd); sets wh (see _rotate_geometry).

    Returns:
        callable: window(r) with r (float) the source distance in pixels, returning
            (win, src_row, src_col): the (2wh+1, 2wh+1) float64 intensity block (package units) and the
            source's fractional (row, col) position in it. window raises RuntimeError if the far
            strips do not cover the window at r.
    """
    N = optics.img_pixel
    c, _, wh, r_near = _rotate_geometry(N, kernel_size)
    band = _band(optics)
    x_src = radial_source_positions(N)[:len(radial_psfs)]
    # Interpolation nodes within r_near: the radial frames at +x_src, and (flipped) at -x_src;
    # near_frame[j] is node j's frame and near_flip[j] whether it is used flipped
    near_x = np.concatenate([-x_src[x_src > 0][::-1], x_src])
    near_frame = np.concatenate([np.flatnonzero(x_src > 0)[::-1], np.arange(len(x_src))])
    near_flip = np.concatenate([np.ones(np.sum(x_src > 0), bool), np.zeros(len(x_src), bool)])
    a0 = int(np.floor(c+0.5)) # window rows are centred on the x axis
    far_x, far_j = axis['far_x'], np.arange(len(axis['far_x']))
    far_flip = np.zeros(len(far_x), bool)
    if len(far_x) and far_x[0] == 0:
        far_x = np.concatenate([-axis['far_x'][:0:-1], far_x])
        far_j = np.concatenate([far_j[:0:-1], far_j])
        far_flip = np.concatenate([np.ones(len(far_j)-len(far_flip), bool), far_flip])
    row0 = axis['far_row0']

    def window(r):
        """
        The PSF of a source on the +x axis at r pixels from the star, on a window around it: the
        radial frames (r <= r_near) or the far strips, combined with band-limited weights in r.

        Args:
            r (float): Source distance from the star along +x, in pixels.

        Returns:
            tuple: (win, src_row, src_col): win the (2 wh + 1, 2 wh + 1) float window (rows = y),
                src_row / src_col the source's (fractional) position inside it, in pixels.

        Raises:
            RuntimeError: a far strip does not cover the window.
        """
        b0 = int(np.floor(c+r+0.5)) # window columns are centred on the pixel nearest the source
        win = 0
        if r <= r_near:
            for j, wj in zip(*bandlimited_weights(near_x, r, band, _N_TAPS)):
                frame = radial_psfs[near_frame[j]]
                if near_flip[j]:
                    win = win+wj*frame[N-1-a0-wh:N-a0+wh, N-1-b0-wh:N-b0+wh][::-1, ::-1]
                else:
                    win = win+wj*frame[a0-wh:a0+wh+1, b0-wh:b0+wh+1]
        else:
            for i, wi in zip(*bandlimited_weights(far_x, r, band, _N_TAPS)):
                strip, col0 = axis['far_strips'][far_j[i]], axis['far_col0'][far_j[i]]
                if far_flip[i]: # source at -x: the strip mirrored through the star
                    rows, cols = N-1-a0-wh-row0, N-1-b0-wh-col0
                    part = strip[rows:rows+2*wh+1, cols:cols+2*wh+1][::-1, ::-1]
                else:
                    rows, cols = a0-wh-row0, b0-wh-col0
                    part = strip[rows:rows+2*wh+1, cols:cols+2*wh+1]
                if min(rows, cols) < 0 or part.shape != (2*wh+1, 2*wh+1):
                    raise RuntimeError("far strip does not cover the window at r=%.2f px" % r)
                win = win+wi*part
        return win, wh+(c-a0), wh+(c+r-b0)
    return window

def _rotate_kernels(window, r, phi, kernel_size, num_cores=-1):
    """
    Kernels for sources at polar positions (r, phi), from x-axis windows. Valid
    for rotationally symmetric optics.

    Samples at the same radius share one interpolated window and one spline
    prefilter; each kernel is then the window sampled on the kernel grid
    (centred on the source, rotated by phi) with a quintic spline. At the
    default sampling (16 px per lambda/D) this keeps both the core and the faint
    wings accurate (~1e-8 of the peak, <~4e-5 of the local wing level), where
    Fourier shifts and rotations of a finite window are limited in the wings by
    the window truncation. Spline accuracy falls quickly with coarser sampling:
    at 8 px per lambda/D expect ~1e-6 of the peak and ~1e-3 (up to a few percent
    in the faintest near-star wings) of the local wing level.

    Args:
        window (callable): From _window_builder.
        r (np.ndarray): (n,) source distances from the star, in pixels.
        phi (np.ndarray): (n,) source angles (radians, from +x (columns) towards +y (rows)).
        kernel_size (int): Kernel width in pixels (odd).
        num_cores (int): Number of threads (joblib n_jobs, thread backend).

    Returns:
        np.ndarray: (n, kernel_size, kernel_size) float32 kernels centred on their sources, in the
            package's intensity units (not unit-sum, not divided by the throughput; negatives clipped to 0).
    """
    h = kernel_size//2
    dy, dx = np.mgrid[-h:h+1, -h:h+1]
    radii, group = np.unique(np.round(r, 9), return_inverse=True)

    def one_radius(g):
        """
        The kernels of every source at the g-th distinct radius: one window, spline-filtered once,
        sampled along each source's rotated kernel grid (quintic spline).

        Args:
            g (int): Index into radii.

        Returns:
            list: [(k, kernel)]: k the source index, kernel a (kernel_size, kernel_size) float array
                (negatives clipped to 0).
        """
        win, src_row, src_col = window(radii[g])
        coeffs = spline_filter(win, order=5, mode='mirror')
        out = []
        for k in np.flatnonzero(group == g):
            u = dx*np.cos(phi[k])+dy*np.sin(phi[k])  # kernel offsets rotated back onto the +x axis
            v = -dx*np.sin(phi[k])+dy*np.cos(phi[k])
            kernel = map_coordinates(coeffs, [src_row+v, src_col+u], order=5, mode='mirror', prefilter=False)
            out.append((k, np.maximum(kernel, 0)))
        return out

    kernels = np.empty((len(r), kernel_size, kernel_size), dtype=np.float32)
    for part in Parallel(n_jobs=num_cores, prefer='threads', return_as='generator_unordered')(
            delayed(one_radius)(g) for g in range(len(radii))):
        for k, kernel in part:
            kernels[k] = kernel
    return kernels

def _rates_from_kernels(r_metric, s0, s1):
    """
    Shape-change rates from unit-sum kernels sampled along one direction.

    d_r at each radius is the larger of the two neighbouring finite differences; both rates are then
    smoothed with a running maximum over 3 radii (conservative).

    Args:
        r_metric (np.ndarray): (n,) increasing source radii, in pixels.
        s0 (np.ndarray): (n, K, K) unit-sum kernels at (r_metric, phi).
        s1 (np.ndarray): (n, K, K) unit-sum kernels at (r_metric, phi + _DPHI).

    Returns:
        tuple: d_r ((n,), L1 change per pixel of radius) and d_phi ((n,), L1 change per radian).
    """
    d_r_mid = np.abs(np.diff(s0, axis=0)).sum(axis=(1, 2))/np.diff(r_metric)
    d_r = np.maximum(np.concatenate([[d_r_mid[0]], d_r_mid]), np.concatenate([d_r_mid, [d_r_mid[-1]]]))
    d_phi = np.abs(s1-s0).sum(axis=(1, 2))/_DPHI
    running_max = lambda a: np.maximum.reduce([np.r_[a[1:], a[-1]], a, np.r_[a[0], a[:-1]]])
    return running_max(d_r), running_max(d_phi)

def _unit_sum(kernels):
    """
    Kernels divided by their sums (over the last two axes; sums floored at 1e-300).

    Args:
        kernels (np.ndarray): (..., K, K) kernels.

    Returns:
        np.ndarray: Same shape, each kernel summing to 1 (all-zero kernels stay 0).
    """
    return kernels/np.maximum(kernels.sum(axis=(-2, -1), keepdims=True), 1e-300)

def shape_change_rates(window, optics, kernel_size, num_cores=-1):
    """
    How fast the PSF shape changes with source radius and angle, for
    rotationally symmetric optics (from x-axis PSFs).

    d_r(r) is the L1 difference between unit-sum kernels at neighbouring radii
    per pixel of radius, and d_phi(r) the L1 difference between the kernel at
    (r, 0) and at (r, 2 deg) per radian. Both are smoothed with a running
    maximum (conservative). Radii are 0.5 px apart out to 3 lambda/D (or the
    reach of the radial frames) and at the axis node spacing beyond.

    Args:
        window (callable): From _window_builder.
        optics (Optics): The optical system and image sampling.
        kernel_size (int): Kernel width in pixels (odd).
        num_cores (int): Number of threads.

    Returns:
        tuple: r_metric ((n,) radii in pixels, out to past the image corners), d_r ((n,) L1 change
            per pixel of radius), d_phi ((n,) L1 change per radian); float64.
    """
    c, _, _, r_near = _rotate_geometry(optics.img_pixel, kernel_size)
    pixel_scale = optics.pixel_scale
    step = 1/(4*_band(optics))
    r_top = np.sqrt(2)*c+_MAX_SPACING/pixel_scale
    r_fine = max(r_near, min(3/pixel_scale, r_top)) # 0.5 px steps out to 3 lambda/D at least
    r_metric = np.concatenate([np.arange(0.25, r_fine, 0.5), np.arange(r_fine, r_top+step, step)])
    # Kernels at (r, 0) and (r, _DPHI) for every r, interleaved
    shapes = _unit_sum(_rotate_kernels(window, np.repeat(r_metric, 2), np.tile([0, _DPHI], len(r_metric)),
                                       kernel_size, num_cores))
    return (r_metric,)+_rates_from_kernels(r_metric, shapes[0::2], shapes[1::2])

def shape_change_rates_direct(optics, kernel_size, num_cores=-1, n_dirs=4):
    """
    How fast the PSF shape changes with source radius and angle, for any optics
    (from PSFs computed directly along several directions).

    Same quantities as shape_change_rates, measured along n_dirs directions (off
    the pixel axes) and taking the largest, so the grid is dense enough in every
    direction. Radii are 1 px apart out to 3 lambda/D and 0.25 lambda/D beyond.

    Args:
        optics (Optics): The optical system and image sampling.
        kernel_size (int): Kernel width in pixels (odd).
        num_cores (int): Number of CPU cores for the propagation.
        n_dirs (int): Number of directions.

    Returns:
        tuple: r_metric ((n,) radii in pixels), d_r ((n,) L1 change per pixel of radius) and d_phi
            ((n,) L1 change per radian), the maxima over the directions; float64.
    """
    pixel_scale = optics.pixel_scale
    c = optics.centre
    r_top = np.sqrt(2)*c+_MAX_SPACING/pixel_scale
    r_metric = np.concatenate([np.arange(0.5, 3/pixel_scale, 1.0),
                               np.arange(3/pixel_scale, r_top+0.25/pixel_scale, 0.25/pixel_scale)])
    angles = 0.1+2*np.pi*np.arange(n_dirs)/n_dirs  # offset by 0.1 rad from the pixel axes
    jobs = [(a+da, rr) for a in angles for rr in r_metric for da in (0, _DPHI)]
    results = Parallel(n_jobs=num_cores)(
        delayed(_spatial_psf_worker)(c+rr*np.sin(a), c+rr*np.cos(a), kernel_size, optics) for a, rr in jobs)
    shapes = _unit_sum(np.array([k for _, k in results])).reshape(n_dirs, len(r_metric), 2, kernel_size, kernel_size)
    rates = [_rates_from_kernels(r_metric, shapes[i, :, 0], shapes[i, :, 1]) for i in range(n_dirs)]
    return r_metric, np.max([d for d, _ in rates], axis=0), np.max([p for _, p in rates], axis=0)

def choose_kernel_size(optics, pre_img=None, error_target=0.02, num_cores=-1, n_probe=48, crop_share=0.8,
                       seed=0, verbose=True, radial_psfs=None):
    """
    The smallest kernel size whose crop keeps the image error within a target.

    A kernel of K x K pixels misses each PSF's light that lands on the image
    outside the kernel window. That loss is measured, not modelled: PSFs of the
    actual optics are computed on the whole image at n_probe pixels drawn where
    pre_img is bright, and for every odd K the light on the image outside the
    K x K window around each probe is summed. Weighted by each probe's light on
    the image, this is the crop error of pre_img's final image. The smallest K
    keeping it within crop_share * error_target is returned; the rest of the
    target is left for the other error sources (blending, transmission map),
    which the error report checks afterwards.

    With radial_psfs (circularly symmetric optics), a probe's PSF is the radial PSF rotated to
    it, as the circular method places planets (_cir_planet): no new propagation. Against
    propagated probes the kernel size was the same in 207 of 216 cases (4-16 px per lambda/D,
    charges 2 and 6, rings, arcs, 1/r^3 profiles) and one step (2 px) apart otherwise, where the
    true crop error sits at the threshold (within 1% of it). Probes beyond the last radial PSF
    (the image corners) are still propagated.

    Args:
        optics (Optics): The optical system and image sampling.
        pre_img (np.ndarray): (img_pixel, img_pixel) target image (any unit; its absolute value
            weighs the probes); None, or an all-zero image, weighs every pixel equally.
        error_target (float): Target relative error of the final image (e.g. 0.02 for 2%).
        num_cores (int): Number of CPU cores for the propagation.
        n_probe (int): Number of probe pixels drawn (with replacement; repeats count several times).
        crop_share (float): Share of the error target given to the crop.
        seed (int): Random seed of the probe pixels.
        verbose (bool): Print the choice.
        radial_psfs (np.ndarray or None): ((img_pixel+1)//2, S, S) radial PSF series of circularly
            symmetric optics, on frames of S >= extended_frame_size(img_pixel) (smaller frames would
            lose the wings rotated into the image corners: then every probe is propagated); None:
            propagate every probe.

    Returns:
        tuple: the kernel size (int, odd, pixels; 2*img_pixel-1 if no size meets the target), and a
            dict with kernel_sizes ((img_pixel,) int: 1, 3, ..., 2*img_pixel-1) and crop (their
            estimated crop error, a fraction of the image light).
    """
    N = optics.img_pixel
    weight = np.ones((N, N)) if pre_img is None else np.abs(np.asarray(pre_img, dtype=float))
    if not weight.sum() > 0:
        weight = np.ones((N, N))
    picks = np.random.default_rng(seed).choice(N*N, n_probe, p=weight.ravel()/weight.sum())
    probes, counts = np.unique(picks, return_counts=True)
    c = (N-1)/2
    r = np.hypot(probes//N-c, probes%N-c)
    rotated = np.zeros(len(probes), bool)
    if radial_psfs is not None and radial_psfs.shape[-1] >= extended_frame_size(N):
        rotated = r < radial_source_positions(N)[len(radial_psfs)-1]
    frames = dict(zip(probes[~rotated], Parallel(n_jobs=num_cores)(
        delayed(_probe_worker)(*divmod(p, N), optics) for p in probes[~rotated])))
    frames.update({p: _cir_planet(radial_psfs, N, rp, np.arctan2(p//N-c, p%N-c))
                   for p, rp in zip(probes[rotated], r[rotated])})

    kernel_sizes = np.arange(1, 2*N, 2)
    h = kernel_sizes//2
    lost = np.zeros(len(kernel_sizes))
    light = 0.
    for p, count in zip(probes, counts):
        frame = frames[p]
        row, col = divmod(p, N)
        table = np.pad(frame.cumsum(0).cumsum(1), ((1, 0), (1, 0))) # summed-area table
        r0, r1 = np.maximum(row-h, 0), np.minimum(row+h+1, N)
        c0, c1 = np.maximum(col-h, 0), np.minimum(col+h+1, N)
        inside = table[r1, c1]-table[r0, c1]-table[r1, c0]+table[r0, c0]
        lost += count*(frame.sum()-inside)
        light += count*frame.sum()
    crop = lost/light
    ok = np.flatnonzero(crop <= crop_share*error_target)
    kernel_size = int(kernel_sizes[ok[0]] if len(ok) else kernel_sizes[-1])
    if verbose:
        print("Kernel size K=%d (+-%.1f lambda/D): estimated crop error %.2g%% of the image light "
              "(%.0f%% of the %.2g%% error target; %d probe pixels)"
              % (kernel_size, (kernel_size//2)*optics.pixel_scale, 100*crop[kernel_sizes == kernel_size][0],
                 100*crop_share, 100*error_target, len(probes)))
    return kernel_size, dict(kernel_sizes=kernel_sizes, crop=crop)

def _axis_throughput(optics, axis, radii):
    """
    Throughput at the given source radii, for rotationally symmetric, band-limited optics.

    Band-limited interpolation of the x-axis runs; the throughput is even in x, so the runs are
    mirrored to -x and both sides of the star are used.

    Args:
        optics (Optics): The optical system and image sampling.
        axis (dict): From axis_products (axis_x, throughput).
        radii (array-like): (n,) source distances from the star, in pixels.

    Returns:
        np.ndarray: (n,) float64 throughput.
    """
    band = _band(optics)
    sym_x = np.concatenate([-axis['axis_x'][:0:-1], axis['axis_x']])
    sym_t = np.concatenate([axis['throughput'][:0:-1], axis['throughput']])
    return np.array([w@sym_t[i] for i, w in (bandlimited_weights(sym_x, rk, band, _N_TAPS) for rk in radii)])

def spatial_psf_calculation(optics, num_cores=-1, radial_psfs=None, kernel_size=None, error_target=0.02,
                            pre_img=None, tol=None, mode=None, axis_cache=None, save_path=None,
                            validate=True, n_validate=6, max_kernel_gb=4.0):
    """
    Calculates the spatially-sampled PSF kernels, sample grid and blend weights
    for the spatial convolution method.

    The sample grid (equal_error_grid) is built from how fast the PSF shape
    actually changes, so neighbouring samples differ by about tol; inside
    ~1-3 lambda/D every pixel is its own sample. Each kernel is centred on its
    source and divided by the source's throughput, so it carries the PSF shape
    only; the throughput goes into the transmission map tmap: each pixel's own
    throughput for symmetric, band-limited optics, otherwise the samples'
    throughput blended onto the pixels like the kernels. Their product
    reproduces the PSF at every sample.

    Args:
        optics (Optics): The optical system and image sampling (see toycoronagraph.optics).
        num_cores (int): Number of CPU cores to use for parallel calculation.
        radial_psfs (np.ndarray): ((img_pixel+1)//2, S, S) radial PSF series (
            radial_psf_series; S >= img_pixel, larger frames are cropped); computed if needed and
            not given. Used for rotationally symmetric, band-limited optics, to measure the grid's
            rates and by 'rotate'.
        kernel_size (int): Kernel width in pixels (odd); None chooses the smallest one
            meeting error_target for pre_img (see choose_kernel_size).
        error_target (float): Target relative error of the final image (e.g. 0.02 for 2%).
        pre_img (np.ndarray): (img_pixel, img_pixel) target image, used only to choose the kernel size.
        tol (float): Target PSF-shape change (L1 of unit-sum kernels) between neighbouring
            samples; None derives it from error_target (see tol_for_error_target).
        mode (str): 'rotate' interpolates the x-axis PSFs in source position
            (band-limited weights) and rotates them (quintic spline): fast, only for
            optics that are rotationally symmetric with a known band limit.
            'direct' propagates every sample (any optics, much slower). None picks
            'rotate' when the optics allow it.
        axis_cache (str): .npz path to cache the x-axis runs (see axis_products); None for no cache.
        save_path (str): If given, save the result there with np.savez.
        validate (bool): Check the result against exact PSFs and print the error
            analysis per stage (see spatial_error_report); stored in the result.
        n_validate (int): Samples and pixels checked per radius band.
        max_kernel_gb (float): Refuse to build a kernel cube larger than this (GB).

    Returns:
        dict: psfs ((n, K, K) float32 kernels centred on their sources, divided by their throughput;
            0 where the throughput is <= 1e-12), throughput ((n,) float64), tmap ((img_pixel, img_pixel)
            float64 transmission map), the grid of equal_error_grid (positions ((n, 2) (row, col) pixels),
            rings, ring_radii, ring_theta0, ring_ntheta, r_pix, and the blend weights inds, wts
            ((img_pixel, img_pixel, 4))), the shape-change rates r_metric, d_r, d_phi, the settings
            (optics.settings(), mode, kernel_size, tol, error_target) and, if validated, the err_* entries
            of spatial_error_report.

    Raises:
        ValueError: If mode is unknown, 'rotate' is not allowed for the optics, or kernel_size is even.
        MemoryError: If the kernel cube would exceed max_kernel_gb.
        RuntimeError: If 'direct' kernels are not centred on their sources (tilt convention check).
    """
    backend = optics.backend()
    rotate_ok = bool(backend.rotationally_symmetric and backend.bandlimit)
    mode = ('rotate' if rotate_ok else 'direct') if mode is None else mode
    if mode not in ('direct', 'rotate'):
        raise ValueError("mode must be 'direct' or 'rotate'")
    if mode == 'rotate' and not rotate_ok:
        raise ValueError("mode='rotate' needs rotationally symmetric optics with a known band limit; "
                         "use mode='direct'")
    if kernel_size is None:
        kernel_size, _ = choose_kernel_size(optics, pre_img, error_target, num_cores,
                                            radial_psfs=radial_psfs if mode == 'rotate' else None)
    if tol is None:
        tol = tol_for_error_target(error_target)
    if kernel_size%2 != 1:
        raise ValueError("kernel_size must be odd")

    N, pixel_scale = optics.img_pixel, optics.pixel_scale
    if rotate_ok:
        if radial_psfs is None:
            radial_psfs = radial_psf_series(optics, num_cores)
        radial_psfs = _crop_centre(radial_psfs, N) # larger (extended) frames: their image-size centre
        axis = axis_products(optics, num_cores, kernel_size, axis_cache)
        window = _window_builder(radial_psfs, axis, optics, kernel_size)
        r_metric, d_r, d_phi = shape_change_rates(window, optics, kernel_size, num_cores)
    else:
        r_metric, d_r, d_phi = shape_change_rates_direct(optics, kernel_size, num_cores)
    grid = equal_error_grid(N, r_metric, d_r, d_phi, tol, _MAX_SPACING/pixel_scale,
                            _N_THETA_MIN, _R_INNER_MAX/pixel_scale)
    positions = grid['positions']
    r, phi = positions_to_polar(N, positions)
    cube_gb = len(r)*kernel_size**2*4/1e9
    print("%s mode: %d PSF samples (%d per-pixel within %.2f lambda/D, %d rings), K=%d, kernel cube %.2f GB"
          % (mode, len(r), np.sum(grid['rings'] == 0), grid['r_pix']*pixel_scale, len(grid['ring_radii']),
             kernel_size, cube_gb))
    if cube_gb > max_kernel_gb:
        raise MemoryError("the kernel cube would need %.2f GB (K=%d, %d samples) > max_kernel_gb=%.2f; "
                          "raise max_kernel_gb, or loosen error_target (smaller K) or tol (fewer samples)"
                          % (cube_gb, kernel_size, len(r), max_kernel_gb))

    if mode == 'direct':
        kernels = np.empty((len(r), kernel_size, kernel_size), dtype=np.float32)
        throughput = np.empty(len(r))
        for k, (t, kernel) in enumerate(Parallel(n_jobs=num_cores)(
                delayed(_spatial_psf_worker)(row, col, kernel_size, optics) for row, col in positions)):
            throughput[k], kernels[k] = t, kernel
        # Sanity check: well off-axis, the PSF peak must sit on the source
        h = kernel_size//2
        off_axis = np.flatnonzero(r*pixel_scale > 3)
        peaks = np.array([np.unravel_index(np.argmax(kernels[k]), kernels[k].shape) for k in off_axis])
        if len(peaks) and np.abs(peaks-h).max() > 2:
            raise RuntimeError("Direct PSFs are not centred on their sources; check the tilt convention")
    else:
        kernels = _rotate_kernels(window, r, phi, kernel_size, num_cores)
        throughput = _axis_throughput(optics, axis, r)

    # Kernels carry the PSF shape (normalized in place); the throughput goes into the transmission map
    for k in range(len(r)):
        if throughput[k] > 1e-12:
            kernels[k] /= throughput[k]
        else:
            kernels[k] = 0

    # Transmission map: for symmetric, band-limited optics the throughput of every pixel's own
    # radius (exact, once per distinct radius); otherwise the samples' throughput, blended
    if rotate_ok:
        rows, cols = np.indices((N, N))
        pix_radii, which = np.unique(np.round(np.hypot(rows-optics.centre, cols-optics.centre), 9),
                                     return_inverse=True)
        tmap = _axis_throughput(optics, axis, pix_radii)[which].reshape(N, N)
    else:
        tmap = blend_map(grid['inds'], grid['wts'], throughput)

    data = dict(psfs=kernels, throughput=throughput, tmap=tmap,
                r_metric=r_metric, d_r=d_r, d_phi=d_phi, **grid, **optics.settings(),
                mode=mode, kernel_size=kernel_size, tol=tol, error_target=error_target)
    if validate:
        data.update(spatial_error_report(data, num_cores, n_validate))
    if save_path is not None:
        np.savez(save_path, **data)
    return data

def load_spatial_psfs(path):
    """
    Load spatial PSF data saved by spatial_psf_calculation (or axis_products).

    Args:
        path (str): Path of the .npz file.

    Returns:
        dict: The saved arrays, with 0-d entries (scalar settings) converted back to Python values.
    """
    with np.load(path) as f:
        return {k: (f[k].item() if f[k].ndim == 0 else f[k]) for k in f.files}

def spatial_error_report(model, num_cores=-1, n_per_band=6, seed=0, verbose=True):
    """
    Error analysis of the spatial PSF model, stage by stage, against exact PSFs
    (propagated directly through the model's optics).

    Checks a few random samples in each radius band (per-pixel zone, r_pix-3,
    3-8 and >8 lambda/D; empty bands are dropped), and pixels spread evenly in
    radius over each band, so the image estimate can follow the error with
    radius. Errors are L1 differences of the PSF on the image, relative to the
    PSF's light that lands on the image (light falling outside the image is not
    part of it, so it is not an error):
      1. kernel generation: a sample's kernel x throughput vs the exact PSF at
         the sample, over the kernel window;
      2. blending between samples: the blended kernel (with the exact throughput)
         vs the exact PSF of a source at the pixel;
      3. transmission map: its error, times the light in the kernel;
      4. kernel crop: the PSF's light on the image but outside the kernel window;
      5. cumulated: the pixel's effective PSF (transmission map x blended kernel,
         nothing outside the window) vs the exact PSF, end to end.
    Stage 2 includes the kernel generation error of the neighbouring samples,
    and the cumulated error is at most about the sum of stages 2-4.

    Args:
        model (dict): From spatial_psf_calculation (kernels, grid, tmap and settings).
        num_cores (int): Number of CPU cores for the hcipy runs.
        n_per_band (int): Samples and pixels checked per radius band.
        seed (int): Random seed of the choice of samples and pixels.
        verbose (bool): Print the table (largest error per stage and band).

    Returns:
        dict: err_edges ((n_band+1,) band edges, pixels), err_stages ((5,) stage names), err_max and
            err_median ((5, n_band) relative errors; NaN where a band has no check), err_n ((n_band,)
            pixels checked per band), err_pix_r ((n_pix,) the checked pixels' radii, pixels) and err_pix
            ((5, n_pix) their errors per stage; NaN for the kernel generation, which is checked at samples).
    """
    optics = Optics.from_settings(model)
    N, K = model['img_pixel'], model['kernel_size']
    pixel_scale = optics.pixel_scale
    c = optics.centre
    edges = np.unique([0, model['r_pix'], max(3/pixel_scale, model['r_pix']), 8/pixel_scale, np.inf])
    stages = ["kernel generation (%s)" % model['mode'], "blending between samples", "transmission map",
              "kernel crop (K=%d)" % K, "cumulated (end to end)"]
    h = K//2
    rng = np.random.default_rng(seed)
    rows, cols = np.indices((N, N))
    pix_r = np.hypot(rows-c, cols-c).ravel()
    smp_r, _ = positions_to_polar(N, model['positions'])
    # Checks to run: (band, kind, index, row, col); index is the sample index for "sample",
    # the flat pixel index for "pixel"
    jobs = []
    for b in range(len(edges)-1):
        for kind, radii in [("sample", smp_r), ("pixel", pix_r)]:
            cand = np.flatnonzero((radii >= edges[b]) & (radii < edges[b+1]))
            if kind == "pixel" and len(cand) > n_per_band:
                # spread evenly in radius: a random pixel in each of n_per_band equal-width rings
                r_lo, r_hi = radii[cand].min(), radii[cand].max()
                ring = np.minimum(((radii[cand]-r_lo)/max(r_hi-r_lo, 1e-9)*n_per_band).astype(int), n_per_band-1)
                picks = [rng.choice(cand[ring == k]) for k in range(n_per_band) if np.any(ring == k)]
            else:
                picks = rng.choice(cand, min(n_per_band, len(cand)), replace=False) if len(cand) else []
            for i in picks:
                row, col = model['positions'][i] if kind == "sample" else divmod(i, N)
                jobs.append((b, kind, i, row, col))
    results = Parallel(n_jobs=num_cores)(
        delayed(_validation_worker)(row, col, K, optics, full_frame=(kind == "pixel"))
        for _, kind, _, row, col in jobs)

    errors = [[[] for _ in range(len(edges)-1)] for _ in stages]
    pix_r_checked, pix_errors = [], []
    for (b, kind, i, row, col), res in zip(jobs, results):
        if kind == "sample":
            _, exact, _ = res
            kernel = model['psfs'][i].astype(float)*model['throughput'][i]
            errors[0][b].append(np.abs(kernel-exact).sum()/exact.sum())
            continue
        T, _, _, frame = res  # exact throughput and PSF on the image
        # The kernel window around the pixel, clipped to the image
        r0, r1, c0, c1 = max(row-h, 0), min(row+h+1, N), max(col-h, 0), min(col+h+1, N)
        mix = np.einsum('k,kij->ij', model['wts'][row, col], model['psfs'][model['inds'][row, col]].astype(float))
        mix = mix[r0-row+h:r1-row+h, c0-col+h:c1-col+h]
        exact_box, on_image = frame[r0:r1, c0:c1], frame.sum()
        tmap = model['tmap'][row, col]
        errors[1][b].append(np.abs(T*mix-exact_box).sum()/on_image)
        errors[2][b].append(abs(tmap-T)*mix.sum()/on_image)
        errors[3][b].append((on_image-exact_box.sum())/on_image)
        errors[4][b].append((np.abs(tmap*mix-exact_box).sum()+on_image-exact_box.sum())/on_image)
        pix_r_checked.append(np.hypot(row-c, col-c))
        pix_errors.append([np.nan]+[errors[s][b][-1] for s in range(1, len(stages))])
    stat = lambda f: np.array([[f(e) if len(e) else np.nan for e in stage] for stage in errors])
    report = dict(err_edges=edges, err_stages=np.array(stages), err_max=stat(np.max), err_median=stat(np.median),
                  err_n=np.array([len(e) for e in errors[1]]), err_pix_r=np.array(pix_r_checked),
                  err_pix=np.array(pix_errors).reshape(-1, len(stages)).T)
    if verbose:
        bands = (["<%.2g l/D (per-pixel)" % (edges[1]*pixel_scale)]
                 + ["%.2g-%.2g l/D" % (edges[b]*pixel_scale, edges[b+1]*pixel_scale) for b in range(1, len(edges)-2)]
                 + [">%.2g l/D" % (edges[-2]*pixel_scale)])
        print("Spatial PSF error analysis vs exact PSFs (max over %d samples/pixels per band;"
              " L1 of the PSF on the image, relative to its light on the image):" % n_per_band)
        print("  %-28s" % "stage" + "".join("%22s" % b for b in bands))
        for s, row in zip(stages, report['err_max']):
            print("  %-28s" % s + "".join("%22s" % ("-" if np.isnan(v) else "%.2g%%" % (100*v)) for v in row))
    return report

def spatial_psf_planets(planets_pos, planet_brightness, psf_scale, model):
    """
    Add one or more planets to an image with exact coronagraphic PSFs.

    Each planet's PSF is propagated through the model's optics at its exact
    position (no interpolation), on the full image frame.

    Args:
        planets_pos (list): Planet positions (radius in arcsec, angle in radian from +x towards +y).
        planet_brightness (list): Planet brightnesses, in Jy/arcsec^2 (as Target stores them).
        psf_scale (float): The scale of the PSF, in arcsec per pixel.
        model (dict): From spatial_psf_calculation (gives the optics settings).

    Returns:
        np.ndarray: (img_pixel, img_pixel) float64 image of the planets (Jy/arcsec^2); planets
            outside the image are skipped with a message.
    """
    optics = Optics.from_settings(model)
    pixels = [arcsec_to_pixel(pos, psf_scale, optics.img_pixel) for pos in planets_pos]
    return exact_point_sources(optics, pixels, planet_brightness)

def exact_point_sources(optics, pixels, values):
    """
    Point sources through the optics, each propagated at its exact position on the full frame.

    Args:
        optics (Optics): The optical system and image grid.
        pixels (list): (row, column) of each source on the image, in pixels (need not be integers).
        values (list): Each source's flux per pixel area (e.g. Jy/arcsec^2: Jy / pixel_scale^2).

    Returns:
        np.ndarray: (img_pixel, img_pixel) float64 image of the sources, in the unit of values.
            Sources more than half a pixel outside the image are skipped with a message.
    """
    N = optics.img_pixel
    image = np.zeros((N, N))
    for i, ((row, col), value) in enumerate(zip(pixels, values)):
        if not (-0.5 <= row <= N-0.5 and -0.5 <= col <= N-0.5):
            print("Planet #"+str(i+1)+" is outside the range of the plot")
            continue
        image += value*_source_at_pixel(optics, row, col)[0].frame()
    return image

def spatial_psf_contrast(pre_img, planet_pos, planet_brightness, psf_scale, model, use_gpu=False, ncores=-2,
                         report=True):
    """
    Find the planet brightness and background brightness in the final image, with spatially-varying PSFs.

    Same measurement as cir_psf_contrast: the planet's core (pixels above half
    its peak) is the aperture; the planet and the dust (with and without the
    target inside the IWA) are summed inside it.

    Args:
        pre_img (np.ndarray): (img_pixel, img_pixel) input image, in Jy/arcsec^2.
        planet_pos (array-like): Planet position (radius in arcsec, angle in radian).
        planet_brightness (float): Planet brightness in Jy/arcsec^2.
        psf_scale (float): The scale of the PSF, in arcsec per pixel.
        model (dict): From spatial_psf_calculation.
        use_gpu (bool): Use CuPy for the FFT convolutions.
        ncores (int): Number of threads.
        report (bool): Print the estimated error of the disk image.

    Returns:
        tuple: planet brightness, background brightness, background brightness (ignored dust inside IWA) in Jy
            (floats: sums over the planet's core aperture times psf_scale^2).
    """
    disk_img = spatial_psf(pre_img, model, False, use_gpu, ncores, report=report)
    disk_img_iwa = spatial_psf(pre_img, model, True, use_gpu, ncores, report=False)
    planet_img = spatial_psf_planets([planet_pos], [planet_brightness], psf_scale, model)
    planet_filter = planet_img > 0.5*planet_img.max()
    planet_b = np.sum(planet_img[planet_filter])*psf_scale**2
    dust_b = np.sum(disk_img[planet_filter])*psf_scale**2
    dust_b_iwa = np.sum(disk_img_iwa[planet_filter])*psf_scale**2
    return planet_b, dust_b, dust_b_iwa

def _crop_centre(img, img_pixel):
    """
    The central img_pixel x img_pixel part of a (possibly larger) frame centred on the star.

    Args:
        img (np.ndarray): (..., S, S) frame(s) with S >= img_pixel and S - img_pixel even.
        img_pixel (int): Number of pixels on one edge of the image.

    Returns:
        np.ndarray: (..., img_pixel, img_pixel) view of img (img itself if S == img_pixel).
    """
    h0 = (img.shape[-1]-img_pixel)//2
    return img[..., h0:h0+img_pixel, h0:h0+img_pixel]

def _cir_chunks(pre_img, img_pixel, psfs, iwa):
    """
    The chunk images of the circular method, to be rotated around the star.

    Each radial PSF is weighted by the light it stands for, and the weighted PSFs are summed.
    Beyond the first _CIR_PIXEL_NODES source radii the weights treat the image as a continuous
    profile: the profile at each source radius (on the +row axis through the star, cubic
    interpolation) times the area of the annulus it stands for, 2 pi r dr; this is the most
    accurate for well-sampled targets. Near the star a profile can be too steep to sample, and
    the star's own pixel gets no annulus area, so the first _CIR_PIXEL_NODES weights take the
    pixels' own light instead: each pixel's light split between the two source radii around its
    own, linearly in radius. The radial PSFs reach the image's edge; the rings beyond, out to the
    corners, use the last radial PSF shifted outwards by whole pixels (far from the star a PSF no
    longer changes with the source's position) and the profile along the diagonal. All rules are
    linear in the image (so is the method). Against brute force (8 px per lambda/D, charges 0 and
    2): smooth rings ~1e-4, an r^-3 wind with light on the star's pixel ~1e-3 to 6e-3, an r^-2.34
    exozodi ~6e-4 (L1, relative).

    Args:
        pre_img (np.ndarray): (img_pixel, img_pixel) circularly symmetric input image.
        img_pixel (int): Number of pixels on one edge of the image.
        psfs (np.ndarray): (n, S, S) radial PSFs (see radial_psf_series), S >= img_pixel.
        iwa (int): Index of the first radial PSF kept in chunk_iwa (0 keeps all, see cir_iwa).

    Returns:
        tuple: chunk and chunk_iwa, (S, S) float64 images (unit of pre_img times pixel area).
    """
    c = (img_pixel-1)/2
    x_src = radial_source_positions(img_pixel)[:len(psfs)]
    image = np.asarray(pre_img, dtype=float)
    profile = map_coordinates(image, [c+x_src, np.full(len(x_src), c)], order=3, mode='nearest')
    weights = profile*2*np.pi*x_src*(x_src-np.concatenate([[0], x_src[:-1]]))
    # near the star: the pixels' own light, split linearly in radius between the two source radii around each
    rows, cols = np.indices(image.shape)
    r, light = np.hypot(rows-c, cols-c).ravel(), image.ravel()
    near = r < x_src[min(_CIR_PIXEL_NODES, len(x_src)-1)]
    r, light = r[near], light[near]
    i = np.clip(np.searchsorted(x_src, r, side='right')-1, 0, len(x_src)-2)
    f = np.clip((r-x_src[i])/(x_src[i+1]-x_src[i]), 0, 1)
    split = np.bincount(i, light*(1-f), len(x_src)) + np.bincount(i+1, light*f, len(x_src))
    weights[:_CIR_PIXEL_NODES] = split[:_CIR_PIXEL_NODES]
    chunk = np.zeros(psfs.shape[-2:])
    chunk_iwa = np.zeros(psfs.shape[-2:])
    for i in np.flatnonzero(weights):
        chunk += weights[i]*psfs[i]
        if i >= iwa:
            chunk_iwa += weights[i]*psfs[i]
    # beyond the last radial PSF, out to the corner pixels (sqrt(2) c from the star): rings every pixel,
    # each the last PSF shifted k pixels along +x, the profile sampled along the diagonal
    shifts = np.arange(1, int(np.ceil(np.sqrt(2)*c-x_src[-1]))+1)
    x_out = x_src[-1]+shifts
    profile = map_coordinates(image, [c+x_out/np.sqrt(2)]*2, order=3, mode='nearest')
    outer = np.zeros(psfs.shape[-2:])
    for k, w in zip(shifts, profile*2*np.pi*x_out):  # annulus width: 1 pixel
        outer[:, k:] += w*psfs[-1][:, :-k]
    return chunk+outer, chunk_iwa+outer

def _cir_rotate_sum(chunk, img_pixel, rot_number):
    """
    Average of rot_number rotations of the chunk about the star (cubic), cropped to the image.

    The rotations are not clipped to the chunk's value range, so the result is
    linear in the chunk (and the circular method linear in the target).

    Args:
        chunk (np.ndarray): (S, S) chunk image centred on the star (see _cir_chunks).
        img_pixel (int): Number of pixels on one edge of the image.
        rot_number (int): Number of rotations, evenly spaced over 360 degrees.

    Returns:
        np.ndarray: (img_pixel, img_pixel) float64 image.
    """
    out = np.zeros(chunk.shape)
    for i in range(rot_number):
        out += rotate(chunk, angle=360*i/rot_number, order=3, clip=False)
    return _crop_centre(out/rot_number, img_pixel)

def _cir_planet(psfs, img_pixel, r, planet_angle):
    """
    The PSF of a planet from the radial PSFs: interpolated to its radius, rotated, cropped.

    The two radial PSFs bracketing r are shifted to r and averaged (frame_interpolation), then
    rotated about the star by planet_angle (cubic; clockwise in array order, so that the
    planet moves towards increasing rows, i.e. +y as in FITS).

    Args:
        psfs (np.ndarray): (n, S, S) radial PSFs (see radial_psf_series), S >= img_pixel.
        img_pixel (int): Number of pixels on one edge of the image.
        r (float): Distance of the planet from the star, in pixels.
        planet_angle (float): Angle of the planet, in radian from +x towards +y (+y = increasing row).

    Returns:
        np.ndarray: (img_pixel, img_pixel) PSF of a planet of unit brightness.
    """
    x_src = radial_source_positions(img_pixel)
    psfs_number = int(np.clip(np.searchsorted(x_src[:len(psfs)], r)-1, 0, len(psfs)-2))
    # (copies: frames of memory-mapped PSFs are read-only, which skimage's warp and rotate do not take)
    frames = [np.array(psfs[psfs_number]), np.array(psfs[psfs_number+1])]
    inter_psf = frame_interpolation(x_src[psfs_number:psfs_number+2], frames, r)
    return _crop_centre(rotate(inter_psf, angle=-planet_angle*180.0/np.pi, order=3), img_pixel)

def _radial(psfs):
    """
    Radial PSFs given as a file or as an array (see cir_psf).

    Args:
        psfs (str or np.ndarray): The .npy file, or the PSFs themselves.

    Returns:
        np.ndarray: (n, S, S) radial PSFs (the file read, or the array as it is).
    """
    return np.load(psfs) if isinstance(psfs, (str, os.PathLike)) else psfs

def cir_psf_contrast(pre_img, r, planet_angle, planet_brightness, psf_scale, img_pixel, rot_number, psfs_name):
    """
    Find the planet brightness and background brightness in the final image.

    The planet's core (pixels above half its peak) is the aperture; the planet and the dust
    (with and without the target inside the IWA, see cir_iwa) are summed inside it.

    Args:
        pre_img (np.ndarray): (img_pixel, img_pixel) circularly symmetric input image, in Jy/arcsec^2.
        r (float): the distance between the planet and the origin, in pixel unit.
        planet_angle (float): The angle between the x-axis and the line of the planet and the origin (radian).
        planet_brightness (float): The planet brightness, in Jy/arcsec^2.
        psf_scale (float): Scale of the PSF, in arcsec per pixel.
        img_pixel (int): Number of pixels on one edge of the image.
        rot_number (int): Number of rotations for generating the circular PSF.
        psfs_name (str or np.ndarray): The .npy file of the radial PSFs (frames of the image size, or larger), or
            the PSFs themselves (e.g. memory-mapped once: Coronagraph.radial_psfs).

    Returns:
        tuple: planet brightness, background brightness, background brightness (ignored dust inside IWA) in Jy.
    """
    # Load PSFs from the specified file
    psfs = _radial(psfs_name)

    # Disk images with and without the dust inside the IWA
    chunk_img, chunk_img_iwa = _cir_chunks(pre_img, img_pixel, psfs, cir_iwa(psfs))
    disk_img = _cir_rotate_sum(chunk_img, img_pixel, rot_number)
    disk_img_iwa = _cir_rotate_sum(chunk_img_iwa, img_pixel, rot_number)

    # The planet image
    planet_img = planet_brightness*_cir_planet(psfs, img_pixel, r, planet_angle)
    threshold = 0.5*planet_img.max()
    planet_filter = np.where(planet_img > threshold , 1, 0)

    # Brightness
    planet_b = np.sum(np.multiply(planet_filter, planet_img))*psf_scale**2
    dust_b = np.sum(np.multiply(planet_filter, disk_img))*psf_scale**2
    dust_b_iwa = np.sum(np.multiply(planet_filter, disk_img_iwa))*psf_scale**2

    return planet_b, dust_b, dust_b_iwa

def cir_psf(pre_img, planets_pos, planet_brightness, psf_scale, iwa_ignore, add_planet, img_pixel, rot_number,
            psfs_name):
    """
    Calculates the final image of a circular symmetric pre-image through circular symmetric PSF with added planets.

    Each radial PSF is weighted by the light it stands for (see _cir_chunks: the
    profile times the annulus area, and near the star the pixels' own light), and
    the PSFs are rotated with cubic interpolation (not clipped, so the disk image is
    linear in pre_img). With radial PSF frames larger than the image
    (radial_psf_series with frames of extended_frame_size) the rotated PSFs also fill the image
    corners; with frames of the image size, the corners (beyond the inscribed
    circle) miss part of their light.

    Args:
        pre_img (np.ndarray): (img_pixel, img_pixel) circularly symmetric input image (any unit, e.g. Jy/arcsec^2).
        planets_pos (list): Planet positions (radius in arcsec, angle in radian).
        planet_brightness (list): Planet brightnesses, in the unit of pre_img.
        psf_scale (float): Scale of the PSF, in arcsec per pixel.
        iwa_ignore (bool): Whether to ignore the target inside the IWA (see cir_iwa).
        add_planet (bool): Whether to add planets to the PSF.
        img_pixel (int): Number of pixels on one edge of the image.
        rot_number (int): Number of rotations for generating the circular PSF.
        psfs_name (str or np.ndarray): The .npy file of the radial PSFs (frames of the image size, or larger), or
            the PSFs themselves (e.g. memory-mapped once: Coronagraph.radial_psfs).

    Returns:
        np.ndarray: (img_pixel, img_pixel) float64 final image, in the unit of pre_img. Planets beyond
            the last radial PSF are skipped with a message.
    """
    # Load PSFs from the specified file
    psfs = _radial(psfs_name)

    # Generate the chunk image along the x-axis, then rotate it around the star
    # (with iwa = 0 the "_iwa" chunk keeps the whole target)
    _, chunk_img_iwa = _cir_chunks(pre_img, img_pixel, psfs, cir_iwa(psfs) if iwa_ignore else 0)
    final_img = _cir_rotate_sum(chunk_img_iwa, img_pixel, rot_number)

    # Add planets if requested
    if add_planet:
        x_max = radial_source_positions(img_pixel)[len(psfs)-1]
        for i in range(len(planets_pos)):
            r = planets_pos[i][0]/psf_scale
            if r < x_max:
                final_img += planet_brightness[i]*_cir_planet(psfs, img_pixel, r, planets_pos[i][1])
            else:
                print("Planet #"+str(i+1)+" is outside the range of the plot")

    return final_img

def cir_psf_planets(planets_pos, planet_brightness, psf_scale, img_pixel, psfs_name):
    """
    Add one or more planets to an image through a series of circular symmetric PSFs.

    Args:
        planets_pos (list): Planet positions (radius in arcsec, angle in radian).
        planet_brightness (list): Planet brightnesses, e.g. in Jy/arcsec^2 (the image has the same unit).
        psf_scale (float): The scale of the PSF, in arcsec per pixel.
        img_pixel (int): The number of pixels on one edge of the image.
        psfs_name (str or np.ndarray): The .npy file of the radial PSFs (frames of the image size, or larger), or
            the PSFs themselves (e.g. memory-mapped once: Coronagraph.radial_psfs).

    Returns:
        np.ndarray: (img_pixel, img_pixel) float64 image of the planets. Planets beyond the last
            radial PSF are silently skipped.
    """
    # Load PSFs from the specified file
    psfs = _radial(psfs_name)

    # Initialize the final image
    planet_img = np.zeros([img_pixel, img_pixel])

    x_max = radial_source_positions(img_pixel)[len(psfs)-1]
    for i in range(len(planets_pos)):
        r = planets_pos[i][0]/psf_scale
        if r < x_max:
            planet_img += planet_brightness[i]*_cir_planet(psfs, img_pixel, r, planets_pos[i][1])
    return planet_img

def core_throughput_fun(psfs):
    """
    Calculate the core throughput of a series of PSFs.

    Args:
        psfs (numpy.ndarray): (n, S, S) normalized PSFs (e.g. the radial PSF series).

    Returns:
        numpy.ndarray: (n,) float64 core throughput of the psfs: the sum of each PSF's pixels above half its peak.
    """
    thresholds = 0.5*psfs.max(axis=(1, 2), keepdims=True)
    return np.where(psfs > thresholds, psfs, 0).sum(axis=(1, 2), dtype=np.float64)

def cir_iwa(psfs):
    """
    Calculate the inner working angle (IWA) of a series of circular symmetric PSFs.

    The IWA is the index of the first radial PSF whose core throughput reaches half of the
    largest one (np.searchsorted, assuming the core throughput increases with radius).

    Args:
        psfs (numpy.ndarray): (n, S, S) radial PSFs (see radial_psf_series).

    Returns:
        int: The IWA of the PSFs, as a radial PSF index (~ distance from the star in pixels).
    """
    core_throughput = core_throughput_fun(psfs)
    iwa = np.searchsorted(core_throughput, 0.5*core_throughput.max())
    return iwa

def cir_core_throughput_plot(psfs, plot_dpi=300, pixel_scale=None):
    """
    Plot the core throughput of a series of circular symmetric PSFs, with the IWA marked.

    Args:
        psfs (numpy.ndarray): (n, S, S) radial PSFs (see radial_psf_series).
        plot_dpi (int): Dots per inch (DPI) for the plot.
        pixel_scale (float): Pixel size in lambda/D (e.g. Optics.pixel_scale): the distance axis in
            lambda/D; None: in pixels (the index of each radial PSF).

    Returns:
        None. Saves Core_throughput.png in the current directory and shows the plot.
    """
    core_throughput = core_throughput_fun(psfs)
    iwa = np.searchsorted(core_throughput, 0.5*core_throughput.max())  # same as cir_iwa(psfs)
    scale, unit = (pixel_scale, "lambda/D") if pixel_scale else (1, "pixel")
    fig = plt.figure(dpi=plot_dpi)
    axs = plt.gca()
    axs.set_ylabel("core throughput (light above half the peak)")
    axs.set_xlabel("distance of the source from the star [%s]" % unit)
    plt.plot(np.arange(len(core_throughput))*scale, core_throughput, color='orange')
    # The IWA: where the core throughput reaches half its largest value
    plt.axvline(iwa*scale, color='black', linestyle='--', label="IWA = %.3g %s" % (iwa*scale, unit))
    plt.legend()
    axs.set_title("Core throughput of an off-axis source; the IWA is where it reaches half its largest value",
                  fontsize=9)
    fig.savefig("Core_throughput.png", format='png', bbox_inches='tight')
    plt.show()
