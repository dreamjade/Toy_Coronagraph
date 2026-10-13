"""
Spatially-varying PSF convolution: sample grid, blend weights and convolution engines.

The spatial method computes the PSF at a set of sample positions on the image
(psf.spatial_psf_calculation builds this "model") and gives every image pixel a
blend of up to four of those samples. This module holds the parts that do not
need the optics:
    * bandlimited_weights: interpolation weights for band-limited functions of
      the source position (used by psf for mode 'rotate');
    * positions_to_polar, equal_error_grid, make_blend_weights, blend_map: the
      sample grid (neighbouring samples differ in PSF shape by about tol) and
      the per-pixel blend weights;
    * convolve_with_blended_psfs: each pixel convolved with its own blended PSF
      (the engine behind spatial_psf);
    * make_psf_inds, convolve_with_spatial_psfs: the nearest-sample engine (one
      PSF per pixel, no blending), kept for direct use (tests/stage1_engine.py);
    * spatial_iwa, spatial_psf, image_error_estimate: the final image of a
      target and the estimate of its error.

Conventions: images are (img_pixel, img_pixel) arrays indexed (row, col); the
star sits at pixel coordinate c = (img_pixel-1)/2 on both axes; sample positions
are (row, col) in pixels (not necessarily integers); PSF kernels have odd sizes
and are centred on their middle pixel, which is the source position.
"""
import numpy as np
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree
from joblib import Parallel, delayed

def bandlimited_weights(nodes, x, band, n_taps=16):
    """
    Interpolation weights for a band-limited function, designed in the frequency domain.

    Picks the n_taps nodes nearest to x and finds weights w (least squares) such
    that sum_j w_j exp(2 pi i f nodes_j) = exp(2 pi i f x) for every frequency
    |f| <= band. Any function whose spectrum lies within the band is then
    interpolated exactly (up to the fit residual), for uniform or non-uniform nodes.
    At a node, that node is returned with weight 1.

    Args:
        nodes (np.ndarray): (n_nodes,) sample positions (any unit, e.g. pixels).
        x (float): Position to interpolate at, in the unit of nodes.
        band (float): Highest frequency present, in cycles per unit of the positions.
        n_taps (int): Number of nearest nodes used (all nodes if there are fewer).

    Returns:
        tuple:
            idx (np.ndarray): (m,) int indices into nodes, sorted (m = min(n_taps, n_nodes),
                or 1 when x is on a node).
            w (np.ndarray): (m,) float64 weights; the value at x is w @ values[idx].
    """
    nodes = np.asarray(nodes, dtype=float)
    # On a node, take it exactly: the least-squares weights are only accurate to
    # ~1e-8 of the neighbouring values, which matters where the function is ~0
    # (e.g. a source exactly on the star, at the coronagraph null)
    hit = np.flatnonzero(np.abs(nodes-x) < 1e-9)
    if len(hit):
        return hit[:1], np.ones(1)
    idx = np.sort(np.argsort(np.abs(nodes-x), kind='stable')[:n_taps])
    # Match exp(2 pi i f d) = cos + i sin at 8 frequencies per tap: real parts to 1, imaginary parts to 0
    f = np.linspace(-band, band, 8*len(idx))[:, None]
    d = nodes[idx]-x
    A = np.concatenate([np.cos(2*np.pi*f*d), np.sin(2*np.pi*f*d)])
    b = np.concatenate([np.ones(len(f)), np.zeros(len(f))])
    return idx, np.linalg.lstsq(A, b, rcond=1e-10)[0]

def positions_to_polar(img_pixel, positions):
    """
    Convert pixel positions to polar coordinates about the star.

    The star sits at pixel coordinate c = (img_pixel-1)/2 on both axes (between
    the four central pixels for even img_pixel), the origin of the focal grid
    of the optics (optics.VortexBackend). The angle is measured from the +column
    axis towards the +row axis, the same frame as hcipy (x = column, y = row).

    Args:
        img_pixel (int): Number of pixels on one edge of the image.
        positions (np.ndarray): (n, 2) array of (row, col) pixel positions (float or int).

    Returns:
        tuple:
            r (np.ndarray): (n,) float64 distance from the star, in pixels.
            phi (np.ndarray): (n,) float64 angle, in radians in [-pi, pi].
    """
    c = (img_pixel-1)/2
    dy = np.asarray(positions)[:, 0]-c
    dx = np.asarray(positions)[:, 1]-c
    return np.hypot(dx, dy), np.arctan2(dy, dx)

def equal_error_grid(img_pixel, r_metric, d_r, d_phi, tol, max_spacing, n_theta_min=8, r_inner_max=None):
    """
    Build the PSF sample grid so neighbouring samples differ in PSF shape by about tol.

    The grid follows how fast the (normalized) PSF shape actually changes:
    d_r(r) is the shape change per pixel of radius and d_phi(r) per radian of
    angle (see psf.shape_change_rates). Ring radii are placed at equal steps
    of the accumulated radial change (integral of d_r = tol between rings), up
    to the first ring beyond the image corner, and each ring gets
    ~2*pi*d_phi(r)/tol samples, so the shape difference between neighbours
    is ~tol everywhere. Two floors keep the far field well covered:
    ring spacing and the arc between samples are at most max_spacing pixels,
    and each ring has at least n_theta_min samples.

    Near the star the required spacing drops below a pixel; inside that radius
    (r_pix) every pixel is its own sample, so those pixels get their exact PSF.
    Ring samples sit at their exact polar positions (not snapped to pixels):
    the convolution applies each kernel around the pixels it serves, so a
    sample's position only sets the kernel's shape. All samples of a ring then
    share one radius, and the ring structure stays regular for blending.

    Every pixel gets up to four samples with blend weights (make_blend_weights);
    samples that no pixel uses are dropped.

    Args:
        img_pixel (int): Number of pixels on one edge of the image.
        r_metric (np.ndarray): (m,) increasing radii, in pixels, at which the rates are given.
        d_r (np.ndarray): (m,) shape change per pixel of radius at r_metric.
        d_phi (np.ndarray): (m,) shape change per radian of angle at r_metric.
        tol (float): Target shape change (L1 of unit-sum kernels) between neighbouring samples.
        max_spacing (float): Largest ring spacing and arc between samples, in pixels.
        n_theta_min (int): Fewest samples on a ring.
        r_inner_max (float): Cap on r_pix, in pixels (bounds the number of per-pixel
            samples); None caps it at the image corner only.

    Returns:
        dict: The grid, with n the number of samples kept and n_rings the number of rings:
            positions ((n, 2) float64): sample positions (row, col), in pixels;
            rings ((n,) int): ring label of each sample: 0 for a per-pixel inner sample,
                k+1 for ring k;
            ring_radii ((n_rings,) float64): ring radii, in pixels (ring 0 at r_pix, the last
                beyond the image corner);
            ring_theta0 ((n_rings,) float64): angle of the first sample of each ring, in radians;
            ring_ntheta ((n_rings,) int): number of samples on each ring before pruning;
            r_pix (float): radius, in pixels, inside which every pixel is its own sample;
            inds ((img_pixel, img_pixel, 4) int), wts ((img_pixel, img_pixel, 4) float64):
                blend weights (see make_blend_weights), with inds re-indexed to the kept
                samples (0 where the weight is 0).
            The ring_* arrays describe the grid before pruning (all rings, all samples).
    """
    c = (img_pixel-1)/2
    rows, cols = np.indices((img_pixel, img_pixel))
    pix_r = np.hypot(rows-c, cols-c)
    r_corner = pix_r.max()

    # Inner radius: beyond it both the ring spacing and the arc between samples are >= 1 pixel
    r_cap = r_corner if r_inner_max is None else r_inner_max
    spacing = np.minimum(tol/np.maximum(d_r, 1e-12), max_spacing)
    arc = np.minimum(r_metric*tol/np.maximum(d_phi, 1e-12), max_spacing)
    too_fine = np.flatnonzero((np.minimum(spacing, arc) < 1) & (r_metric < r_cap))
    r_pix = r_metric[min(too_fine[-1]+1, len(r_metric)-1)] if len(too_fine) else r_metric[0]
    r_pix = min(max(r_pix, 1.0), r_cap) # always include the central pixel(s)

    # Ring radii: equal steps of accumulated shape change s(r) (trapezoid integral of d_r), capped at max_spacing,
    # until a ring lies beyond the corner, so every pixel at or beyond r_pix lies between two rings (also when
    # r_pix is capped at the corner radius)
    s = np.concatenate([[0], np.cumsum(0.5*(d_r[1:]+d_r[:-1])*np.diff(r_metric))])
    radii = [r_pix]
    while radii[-1] <= r_corner:
        s_next = np.interp(radii[-1], r_metric, s)+tol
        r_next = np.interp(s_next, s, r_metric) if s_next < s[-1] else np.inf
        radii.append(min(r_next, radii[-1]+max_spacing))
    ring_radii = np.array(radii)
    ring_ntheta = np.maximum.reduce([np.full(len(ring_radii), n_theta_min),
                                     np.ceil(2*np.pi*np.interp(ring_radii, r_metric, d_phi)/tol),
                                     np.ceil(2*np.pi*ring_radii/max_spacing)]).astype(int)
    ring_theta0 = np.pi*(np.arange(len(ring_radii))%2)/ring_ntheta # stagger alternate rings by half a step

    # Samples: every inner pixel (row-major order), then the rings
    inner = pix_r < r_pix
    pos = [np.stack([rows[inner], cols[inner]], axis=1).astype(float)]
    ring_of = [np.zeros(inner.sum(), dtype=int)]
    for k, (rk, nk, t0) in enumerate(zip(ring_radii, ring_ntheta, ring_theta0)):
        theta = t0+2*np.pi*np.arange(nk)/nk
        pos.append(np.stack([c+rk*np.sin(theta), c+rk*np.cos(theta)], axis=1))
        ring_of.append(np.full(nk, k+1))
    grid = dict(positions=np.concatenate(pos), rings=np.concatenate(ring_of), ring_radii=ring_radii,
                ring_theta0=ring_theta0, ring_ntheta=ring_ntheta, r_pix=r_pix)
    grid['inds'], grid['wts'] = make_blend_weights(img_pixel, grid)

    # Drop samples no pixel uses (e.g. the parts of the last ring beyond the frame)
    used = np.zeros(len(grid['positions']), bool)
    used[grid['inds'][grid['wts'] > 0]] = True
    new_index = np.cumsum(used)-1 # old sample index -> index among the kept samples
    grid['inds'] = np.where(grid['wts'] > 0, new_index[grid['inds']], 0)
    grid['positions'], grid['rings'] = grid['positions'][used], grid['rings'][used]
    return grid

def make_blend_weights(img_pixel, grid):
    """
    Blend weights of the PSF samples for every pixel.

    A pixel inside r_pix uses its own sample (weight 1). Any other pixel lies
    between two rings: on each ring the two samples bracketing the pixel's
    angle are weighted linearly in angle, and the two rings linearly in
    radius, so up to four samples with weights summing to 1. A function that
    is linear in radius is reproduced exactly.

    Args:
        img_pixel (int): Number of pixels on one edge of the image.
        grid (dict): Output of equal_error_grid before pruning: rings, ring_radii,
            ring_theta0, ring_ntheta and r_pix (samples listed as built there:
            inner pixels in row-major order, then ring after ring). Needs at
            least two rings if any pixel lies at or beyond r_pix.

    Returns:
        tuple:
            inds (np.ndarray): (img_pixel, img_pixel, 4) int sample indices; slots 0-1
                on the inner ring, 2-3 on the outer ring (inner pixels: slot 0 only).
            wts (np.ndarray): (img_pixel, img_pixel, 4) float64 weights, summing to 1 per pixel.
    """
    c = (img_pixel-1)/2
    rows, cols = np.indices((img_pixel, img_pixel))
    pix_r = np.hypot(rows-c, cols-c)
    pix_phi = np.arctan2(rows-c, cols-c)
    radii, ntheta, theta0 = grid['ring_radii'], grid['ring_ntheta'], grid['ring_theta0']
    n_inner = np.sum(grid['rings'] == 0)
    ring_start = n_inner+np.concatenate([[0], np.cumsum(ntheta)[:-1]]) # index of each ring's first sample

    inds = np.zeros((img_pixel, img_pixel, 4), dtype=int)
    wts = np.zeros((img_pixel, img_pixel, 4))
    inner = pix_r < grid['r_pix']
    inds[inner, 0] = np.arange(n_inner) # inner samples are listed in the same (row-major) order
    wts[inner, 0] = 1.

    outer = ~inner
    pr, pp = pix_r[outer], pix_phi[outer]
    # Ring lo is the last ring at or inside the pixel's radius (pixels beyond the last ring use the last two)
    lo = np.clip(np.searchsorted(radii, pr, side='right')-1, 0, len(radii)-2)
    w_hi = np.clip((pr-radii[lo])/(radii[lo+1]-radii[lo]), 0, 1)
    for slot, k, w_ring in [(0, lo, 1-w_hi), (2, lo+1, w_hi)]:
        # Position along ring k in units of its sample step: between samples j and j+1 (cyclic)
        t = np.mod(pp-theta0[k], 2*np.pi)*ntheta[k]/(2*np.pi)
        j = np.floor(t).astype(int)%ntheta[k]
        frac = t-np.floor(t)
        inds[outer, slot] = ring_start[k]+j
        wts[outer, slot] = w_ring*(1-frac)
        inds[outer, slot+1] = ring_start[k]+(j+1)%ntheta[k]
        wts[outer, slot+1] = w_ring*frac
    return inds, wts

def blend_map(inds, wts, values):
    """
    Blend per-sample values onto the pixels, e.g. the coronagraph transmission map.

    Args:
        inds (np.ndarray): (img_pixel, img_pixel, 4) int sample indices (from equal_error_grid).
        wts (np.ndarray): (img_pixel, img_pixel, 4) float blend weights (from equal_error_grid).
        values (np.ndarray): (n_samples,) one value per sample.

    Returns:
        np.ndarray: (img_pixel, img_pixel) map, sum_k wts[..., k]*values[inds[..., k]].
    """
    return np.sum(np.asarray(values)[inds]*wts, axis=-1)

def make_psf_inds(img_pixel, positions):
    """
    Assign every pixel to its nearest PSF sample position (for convolve_with_spatial_psfs).

    Args:
        img_pixel (int): Number of pixels on one edge of the image.
        positions (np.ndarray): (n, 2) array of (row, col) sample positions, in pixels.

    Returns:
        np.ndarray: (img_pixel, img_pixel) int array of indices into positions.
    """
    rows, cols = np.indices((img_pixel, img_pixel))
    pixels = np.stack([rows.ravel(), cols.ravel()], axis=1)
    _, inds = cKDTree(positions).query(pixels)
    return inds.reshape(img_pixel, img_pixel)

def psf_convolve_cpu(im, psf):
    """
    Convolve an image with a PSF on the CPU (scipy fftconvolve, mode 'same').

    Args:
        im (np.ndarray): (ny, nx) image.
        psf (np.ndarray): (ky, kx) PSF kernel with odd dimensions, centred at its middle pixel.

    Returns:
        np.ndarray: (ny, nx) convolved image.
    """
    # float64: scipy's FFT keeps float32 input (e.g. the stored kernels) in single precision
    return fftconvolve(im, np.asarray(psf, dtype=np.float64), mode='same')

def psf_convolve_gpu(im, psf):
    """
    Convolve an image with a PSF on the GPU using CuPy (fftconvolve, mode 'same').

    Args:
        im (np.ndarray): (ny, nx) image.
        psf (np.ndarray): (ky, kx) PSF kernel with odd dimensions, centred at its middle pixel.

    Returns:
        np.ndarray: (ny, nx) convolved image, copied back to the host.

    Raises:
        ImportError: If CuPy is not installed.
    """
    try:
        import cupy as cp
        from cupyx.scipy.signal import fftconvolve as fftconvolve_gpu
    except ImportError as err:
        raise ImportError("use_gpu=True requires CuPy (https://cupy.dev); install it or set use_gpu=False") from err
    return cp.asnumpy(fftconvolve_gpu(cp.asarray(im), cp.asarray(psf, dtype=cp.float64), mode='same'))

def _get_psf_patch_bounds(psf_inds, npsf, psf_hw):
    """
    Find the padded bounding box of the pixels assigned to each PSF.

    Args:
        psf_inds (np.ndarray): (ny, nx) int array giving the PSF index of each pixel.
        npsf (int): Number of PSFs (length of the PSF cube); every index must be < npsf.
        psf_hw (tuple): (y, x) half-widths of the PSF kernel, in pixels, used to pad each box.

    Returns:
        tuple:
            used_inds (np.ndarray): (n_used,) int PSF indices present in psf_inds, increasing.
            patch_bounds (np.ndarray): (n_used, 4) int inclusive bounds [y1, y2, x1, x2]
                of each padded box, clipped to the image.
    """
    ny, nx = psf_inds.shape
    yhw, xhw = psf_hw
    rows, cols = np.indices(psf_inds.shape)
    flat_inds = psf_inds.ravel()

    # Bounding box of each PSF region (unused PSFs keep ymax = -1)
    ymin = np.full(npsf, ny, dtype=int)
    ymax = np.full(npsf, -1, dtype=int)
    xmin = np.full(npsf, nx, dtype=int)
    xmax = np.full(npsf, -1, dtype=int)
    np.minimum.at(ymin, flat_inds, rows.ravel())
    np.maximum.at(ymax, flat_inds, rows.ravel())
    np.minimum.at(xmin, flat_inds, cols.ravel())
    np.maximum.at(xmax, flat_inds, cols.ravel())

    # Keep the PSFs that are used, then pad by the kernel half-width
    used_inds = np.flatnonzero(ymax >= 0)
    patch_bounds = np.stack([np.maximum(ymin[used_inds]-yhw, 0),
                             np.minimum(ymax[used_inds]+yhw, ny-1),
                             np.maximum(xmin[used_inds]-xhw, 0),
                             np.minimum(xmax[used_inds]+xhw, nx-1)], axis=1)
    return used_inds, patch_bounds

def _fft_cost(h, w, ky, kx):
    """
    Estimate the work of an FFT convolution of an (h, w) image with a (ky, kx) kernel.

    Args:
        h, w (int or np.ndarray): Image height and width, in pixels (arrays for several patches).
        ky, kx (int): Kernel height and width, in pixels.

    Returns:
        float or np.ndarray: n*log2(n) with n = (h+ky-1)*(w+kx-1), the padded FFT size
            (same shape as h and w).
    """
    n = (h+ky-1)*(w+kx-1)
    return n*np.log2(n)

def _choose_spatial_psf_convolution_strategy(im_shape, psf_shape, patch_bounds):
    """
    Choose between full-frame and localized convolution by comparing FFT work estimates.

    Args:
        im_shape (tuple): (ny, nx) shape of the science image.
        psf_shape (tuple): (ky, kx) shape of one PSF kernel.
        patch_bounds (np.ndarray): (n_used, 4) padded bounding boxes from _get_psf_patch_bounds.

    Returns:
        str: 'full' or 'localized' (the cheaper one; 'full' on a tie).
    """
    ky, kx = psf_shape
    full_cost = len(patch_bounds)*_fft_cost(im_shape[0], im_shape[1], ky, kx)
    heights = patch_bounds[:, 1]-patch_bounds[:, 0]+1
    widths = patch_bounds[:, 3]-patch_bounds[:, 2]+1
    localized_cost = np.sum(_fft_cost(heights, widths, ky, kx))
    return 'localized' if localized_cost < full_cost else 'full'

def convolve_with_spatial_psfs(im_in, psfs, psf_inds, patch_bounds=None, used_inds=None, coron_tmap=None, use_gpu=False,
                               ncores=-2, method='auto'):
    """
    Convolve an image where each pixel uses its own (nearest-sample) PSF.

    The image is split by PSF index: the pixels assigned to PSF i are convolved
    with psfs[i] and the results are summed. 'full' convolves the whole frame
    once per PSF; 'localized' only the padded bounding box of each PSF's pixels.

    Args:
        im_in (np.ndarray): (ny, nx) image to be convolved (any real dtype; not modified).
        psfs (np.ndarray): (n_psf, ky, kx) spatially-sampled PSFs (odd ky, kx). Generally,
            each slice is the PSF of a source at one sampled position, recentred so that
            the source falls on the central pixel.
        psf_inds (np.ndarray): (ny, nx) int array; each entry is the index of the slice
            of psfs with which that pixel of im_in is convolved (e.g. make_psf_inds).
        patch_bounds (np.ndarray): (n_used, 4) pre-computed padded bounding boxes
            (from _get_psf_patch_bounds), matching used_inds. If it or used_inds is None,
            both are computed here when method is 'auto' or 'localized'.
        used_inds (np.ndarray): (n_used,) pre-computed PSF indices present in psf_inds.
            None: computed here.
        coron_tmap (np.ndarray): (ny, nx) coronagraph transmission, by which im_in is
            multiplied before convolution; None for no transmission map.
        use_gpu (bool): Use CuPy for the FFT convolutions (then run serially).
        ncores (int): Number of threads (joblib n_jobs convention: -1 all, -2 all but one)
            when use_gpu is False; 1 runs serially.
        method (str): 'auto', 'localized' or 'full' (case-insensitive). 'localized'
            convolves padded bounding boxes around each PSF region, 'full' the whole
            frame for every PSF index, and 'auto' picks the cheaper by an FFT work estimate.

    Returns:
        np.ndarray: (ny, nx) convolved image (float64 for real input).

    Raises:
        ValueError: If method is not one of the three.
    """
    # Copy as float so integer images and in-place tmap scaling are safe
    im = im_in.astype(np.result_type(im_in.dtype, np.float64))
    if coron_tmap is not None:
        im *= coron_tmap

    convolution_fn = psf_convolve_gpu if use_gpu else psf_convolve_cpu

    psf_yhw, psf_xhw = np.ceil(np.array(psfs.shape[-2:])/2.).astype(int)

    method = method.lower()
    if method not in {'auto', 'localized', 'full'}:
        raise ValueError("method must be one of 'auto', 'localized', or 'full'")

    if method in {'auto', 'localized'}:
        if (used_inds is None) or (patch_bounds is None):
            used_inds, patch_bounds = _get_psf_patch_bounds(psf_inds, psfs.shape[0], (psf_yhw, psf_xhw))
        if method == 'auto':
            method = _choose_spatial_psf_convolution_strategy(im.shape, psfs.shape[-2:], patch_bounds)

    if method == 'full':
        if used_inds is None:
            used_inds = np.flatnonzero(np.bincount(psf_inds.ravel(), minlength=psfs.shape[0]))
        if use_gpu or ncores == 1:
            imcon = np.zeros(im.shape, dtype=im.dtype)
            for i in used_inds:
                im_to_convolve = np.where(psf_inds == i, im, 0.)
                imcon += convolution_fn(im_to_convolve, psfs[i])
        else:
            results = Parallel(n_jobs=ncores, prefer='threads')(
                delayed(convolution_fn)(np.where(psf_inds == i, im, 0.), psfs[i])
                for i in used_inds
                )
            imcon = np.sum(results, axis=0)

    else:
        def _process_patch(i, bounds):
            """Convolve the pixels of PSF i inside its padded box; returns the box-sized result."""
            y1_p, y2_p, x1_p, x2_p = bounds
            patch = im[y1_p:y2_p+1, x1_p:x2_p+1]
            mask_patch = psf_inds[y1_p:y2_p+1, x1_p:x2_p+1] == i
            return convolution_fn(np.where(mask_patch, patch, 0.), psfs[i])

        imcon = np.zeros(im.shape, dtype=im.dtype)
        if use_gpu or ncores == 1:
            for i, bounds in zip(used_inds, patch_bounds):
                y1_p, y2_p, x1_p, x2_p = bounds
                imcon[y1_p:y2_p+1, x1_p:x2_p+1] += _process_patch(i, bounds)
        else:
            results = Parallel(n_jobs=ncores, prefer='threads')(
                delayed(_process_patch)(i, bounds) for i, bounds in zip(used_inds, patch_bounds)
                )
            for patch_convolved, (y1_p, y2_p, x1_p, x2_p) in zip(results, patch_bounds):
                imcon[y1_p:y2_p+1, x1_p:x2_p+1] += patch_convolved

    return imcon

def convolve_with_blended_psfs(im_in, psfs, inds, wts, coron_tmap=None, use_gpu=False, ncores=-2):
    """
    Convolve an image where each pixel uses its own blend of spatially-sampled PSFs.

    Pixel p contributes im[p] * coron_tmap[p] * sum_k wts[p, k] * psfs[inds[p, k]],
    centred on p. Samples used by a single pixel (the per-pixel samples near
    the star) are added directly; every other sample convolves the image
    weighted by its blend weights over the patch where they are non-zero.
    Kernels are never cropped by the image edge before being placed: the
    result is accumulated on a padded canvas and cut to the image at the end.

    Args:
        im_in (np.ndarray): (ny, nx) image to be convolved (not modified).
        psfs (np.ndarray): (n_samples, K, K) sample kernels (odd K, centred; e.g. float32: each is
            convolved in float64).
        inds (np.ndarray): (ny, nx, 4) int sample indices of each pixel's blend.
        wts (np.ndarray): (ny, nx, 4) float blend weights of each pixel.
        coron_tmap (np.ndarray): (ny, nx) coronagraph transmission, by which im_in is
            multiplied first; None for no transmission map.
        use_gpu (bool): Use CuPy for the FFT convolutions (then run serially).
        ncores (int): Number of threads (joblib n_jobs convention: -1 all, -2 all but one).

    Returns:
        np.ndarray: (ny, nx) float64 convolved image.
    """
    im = im_in.astype(np.result_type(im_in.dtype, np.float64))
    if coron_tmap is not None:
        im *= coron_tmap
    ny, nx = im.shape
    ksize = psfs.shape[-1]
    h = ksize//2

    # (sample, pixel, value) triplets for every non-zero contribution, grouped by sample:
    # group g (sample sample[starts[g]]) is the slice starts[g]:starts[g+1]
    values = im.ravel()[:, None]*wts.reshape(-1, inds.shape[-1])
    keep = values != 0
    pixel = np.broadcast_to(np.arange(ny*nx)[:, None], keep.shape)[keep]
    sample = inds.reshape(-1, inds.shape[-1])[keep]
    values = values[keep]
    order = np.argsort(sample, kind='stable')
    sample, pixel, values = sample[order], pixel[order], values[order]
    starts = np.concatenate([[0], np.flatnonzero(np.diff(sample))+1, [len(sample)]])

    # Accumulate on a canvas padded by h so every kernel fits; pixel (r, c) is canvas (r+h, c+h)
    canvas = np.zeros((ny+2*h, nx+2*h))
    single = np.flatnonzero(np.diff(starts) == 1)
    for g in single:
        j = starts[g]
        r, c = divmod(pixel[j], nx)
        canvas[r:r+ksize, c:c+ksize] += values[j]*psfs[sample[j]]

    def patch_convolution(g):
        """Convolve group g's weighted pixels (on their bounding box) with its kernel; returns (r0, c0, result)."""
        sl = slice(starts[g], starts[g+1])
        r, c = np.divmod(pixel[sl], nx)
        r0, c0 = r.min(), c.min()
        patch = np.zeros((r.max()-r0+1, c.max()-c0+1))
        np.add.at(patch, (r-r0, c-c0), values[sl])
        # the kernel in float64: an FFT keeps float32 input (the stored kernels) in single precision, which
        # leaves errors of ~1e-8 of the image's light (up to ~1% in pixels 8 decades below the peak)
        kernel = np.asarray(psfs[sample[starts[g]]], dtype=np.float64)
        if use_gpu:
            import cupy as cp
            from cupyx.scipy.signal import fftconvolve as fftconvolve_gpu
            out = cp.asnumpy(fftconvolve_gpu(cp.asarray(patch), cp.asarray(kernel), mode='full'))
        else:
            out = fftconvolve(patch, kernel, mode='full')
        # mode 'full' extends the patch by h on every side: its corner lands on canvas (r0, c0)
        return r0, c0, out

    multi = np.flatnonzero(np.diff(starts) > 1)
    for r0, c0, out in Parallel(n_jobs=1 if use_gpu else ncores, prefer='threads', return_as='generator_unordered')(
            delayed(patch_convolution)(g) for g in multi):
        canvas[r0:r0+out.shape[0], c0:c0+out.shape[1]] += out
    return canvas[h:h+ny, h:h+nx]

def spatial_iwa(model):
    """
    Inner working angle of spatially-sampled PSFs.

    The IWA is the smallest source radius at which the core throughput (the PSF's
    flux above half its peak, PSF = kernel x throughput) reaches half its maximum;
    the core throughput is made non-decreasing in radius first.

    Args:
        model (dict): From psf.spatial_psf_calculation (uses img_pixel, positions,
            psfs and throughput).

    Returns:
        float: The IWA in pixels (comparable to psf.cir_iwa).
    """
    r, _ = positions_to_polar(model['img_pixel'], model['positions'])
    psf = model['psfs']*model['throughput'][:, None, None]
    core = np.where(psf > 0.5*psf.max(axis=(1, 2), keepdims=True), psf, 0).sum(axis=(1, 2))
    order = np.argsort(r)
    r, core = r[order], np.maximum.accumulate(core[order])
    return float(np.interp(0.5*core.max(), core, r)) if core[0] < 0.5*core.max() else float(r[0])

def spatial_psf(pre_img, model, iwa_ignore=False, use_gpu=False, ncores=-2, report=True):
    """
    Calculate the final image of any (not necessarily circularly symmetric)
    pre-image through spatially-varying coronagraphic PSFs.

    Args:
        pre_img (np.ndarray): (img_pixel, img_pixel) input image (e.g. Jy/arcsec^2);
            the final image is in the same unit.
        model (dict): From psf.spatial_psf_calculation (psfs, inds, wts, tmap, img_pixel, ...).
        iwa_ignore (bool): Set the target to 0 inside the IWA (spatial_iwa) first.
        use_gpu (bool): Use CuPy for the FFT convolutions.
        ncores (int): Number of threads (joblib n_jobs convention).
        report (bool): Print the estimated error of this image (see image_error_estimate),
            if the model has an error report.

    Returns:
        np.ndarray: (img_pixel, img_pixel) float64 final image.
    """
    im = np.array(pre_img, dtype=float)
    if iwa_ignore:
        n = model['img_pixel']
        rows, cols = np.indices((n, n))
        im[np.hypot(rows-(n-1)/2, cols-(n-1)/2) < spatial_iwa(model)] = 0
    if report and 'err_max' in model:
        image_error_estimate(model, im)
    return convolve_with_blended_psfs(im, model['psfs'], model['inds'], model['wts'], model['tmap'], use_gpu, ncores)

def image_error_estimate(model, pre_img, verbose=True):
    """
    Estimated error of the final image of pre_img, from the build-time error check.

    Every pixel's light is spread by its effective PSF, whose error (end to end,
    relative to the PSF's total flux) was measured per radius band when the PSFs
    were built (psf.spatial_error_report). Weighting those errors by the flux
    each pixel sends through the coronagraph (pre_img * tmap) gives the relative
    error of this image's light: 'typical' from the checked pixels' errors,
    interpolated in radius (or the band medians, for older models), 'up to' from
    the band maxima. The same weighting per error stage shows what dominates.

    Args:
        model (dict): From psf.spatial_psf_calculation, with its error report (img_pixel,
            tmap, err_edges, err_stages, err_median, err_max; err_pix_r and err_pix if
            present; error_target, optional).
        pre_img (np.ndarray): (img_pixel, img_pixel) input image.
        verbose (bool): Print the estimate (and a warning above model['error_target']).

    Returns:
        dict: Relative errors (fractions, not %): typical (float), upper (float) and
            per_stage ({stage name: typical error}); empty if pre_img has no flux.
    """
    n = model['img_pixel']
    rows, cols = np.indices((n, n))
    r_pix = np.hypot(rows-(n-1)/2, cols-(n-1)/2)
    n_bands = model['err_max'].shape[1]
    band = np.clip(np.searchsorted(model['err_edges'], r_pix, side='right')-1, 0, n_bands-1)
    flux = np.abs(pre_img)*model['tmap']
    total_flux = flux.sum()
    if total_flux == 0:
        return {}
    weight = np.array([flux[band == b].sum() for b in range(n_bands)])/total_flux # flux fraction per band
    median, upper = np.nan_to_num(model['err_median']), np.nan_to_num(model['err_max'])
    per_stage = median@weight
    if len(model.get('err_pix_r', [])) >= 2:
        # the error changes with radius within a band (e.g. the crop near the star): follow it;
        # stage 0 (kernel generation) is only checked at samples, so it keeps the band medians
        order = np.argsort(model['err_pix_r'])
        for s in range(1, len(per_stage)):
            per_stage[s] = (flux*np.interp(r_pix, model['err_pix_r'][order], model['err_pix'][s][order])).sum()/total_flux
    out = dict(typical=per_stage[-1], upper=(upper@weight)[-1],
               per_stage=dict(zip([str(s) for s in model['err_stages']], per_stage)))
    if verbose:
        dominant = int(np.argmax(per_stage[:-1])) # the last stage is the cumulated error
        print("Estimated error of this image (flux-weighted, from the PSF check): typical %.2g%%, up to %.2g%%; largest part: %s (%.2g%%)"
              % (100*out['typical'], 100*out['upper'], model['err_stages'][dominant], 100*per_stage[dominant]))
        target = model.get('error_target')
        if target is not None and out['typical'] > target:
            print("  WARNING: above the %.2g%% error target; a larger kernel_size (or smaller tol) would reduce it" % (100*target))
    return out
