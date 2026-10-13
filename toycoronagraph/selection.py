"""
Choosing between the circular and the spatial method.

The circular method (psf.cir_psf) is exact only when both the target and the PSFs are
circularly symmetric about the star; the spatial method (spatial.spatial_psf) works for
any target and PSFs. Either symmetry can be set by the user ("symmetric = true/false" in
the parameter file) or, with "auto", found here: the target's by measuring it, the PSFs'
from the optics model (or, for a model that does not declare it, by measuring it).
"""
import numpy as np


def image_asymmetry(pre_img, bin_width=0.02):
    """
    How far a target is from being circularly symmetric about the star.

    In a circularly symmetric image, pixels at the same distance from the star have
    the same value. Pixels are grouped in thin radius bins (bin_width pixels) and the
    asymmetry is the light by which they differ from their bin's mean, relative to
    the target's light (inside the inscribed circle). Near the star a bin holds only
    the pixels at exactly the same radius (mirror images of each other), so even a
    steep, barely resolved profile gives no false signal; further out a bin also mixes
    pixels at slightly different radii, testing the full rotational symmetry. No
    interpolation is involved. Inclined or off-centre disks, arcs, clumps and n-fold
    patterns show up; the value is about the error the circular method would make by
    assuming symmetry.

    Light beyond the inscribed circle (the image corners) is reported separately: the
    circular method reads the target only out to the image edge.

    Args:
        pre_img (np.ndarray): (n, n) target image, star at pixel (n-1)/2 (any unit; the
            absolute value is used).
        bin_width (float): Width of the radius bins, in pixels.

    Returns:
        tuple:
            asymmetry (float): Summed |pixel - its bin's mean| inside the inscribed circle,
                relative to the light there (0 for an exactly symmetric target).
            corner_fraction (float): Light beyond the inscribed circle relative to the total.
        Both are 0 for an empty image.
    """
    img = np.abs(np.asarray(pre_img, dtype=float))
    n = img.shape[0]
    c = (n-1)/2
    rows, cols = np.indices(img.shape)
    rr = np.hypot(rows-c, cols-c)
    total = img.sum()
    if total == 0:
        return 0.0, 0.0
    disk = rr < c # inscribed circle
    r, v = rr[disk], img[disk]
    if v.sum() == 0:
        return 0.0, float(img[~disk].sum()/total)
    b = np.floor(r/bin_width).astype(np.int64)
    mean = np.bincount(b, weights=v)[b]/np.bincount(b)[b] # each pixel's bin mean
    return float(np.abs(v-mean).sum()/v.sum()), float(img[~disk].sum()/total)


def measure_psf_symmetry(optics, kernel_size=None, radii_ld=(0.5, 1.0, 2.0, 4.0), angles_deg=(30.0, 75.0, 140.0),
                         num_cores=-1):
    """
    Measure whether the PSFs are circularly symmetric.

    For a few radii and angles, the PSF of a source at (r, phi) (propagated through
    the optics) is compared with the PSF of a source at (r, 0) rotated by phi (quintic
    spline), over a kernel_size window centred on the source. The reference PSF is
    computed on a window twice as wide, so the rotation never reads beyond it.

    Args:
        optics (Optics): The optical system and image sampling.
        kernel_size (int): Comparison window, in pixels (odd); None: about +-3 lambda/D.
        radii_ld (tuple): Source radii, in lambda/D.
        angles_deg (tuple): Source angles, in degrees (the reference is at 0).
        num_cores (int): Number of CPU processes (joblib n_jobs; -1 uses all).

    Returns:
        float: The largest L1 difference, relative to the exact PSF's light in the window.
    """
    from joblib import Parallel, delayed
    from scipy.ndimage import map_coordinates
    from toycoronagraph.psf import _spatial_psf_worker
    s, c = optics.pixel_scale, optics.centre
    if kernel_size is None:
        kernel_size = 2*int(round(3/s))+1
    h = kernel_size//2
    # For each radius (in pixels): the reference at angle 0, then every angle in angles_deg
    jobs = [(r/s, np.radians(a)) for r in radii_ld for a in (0.0,)+tuple(angles_deg)]
    results = Parallel(n_jobs=num_cores)(
        delayed(_spatial_psf_worker)(c+r*np.sin(phi), c+r*np.cos(phi), kernel_size+2*h, optics) for r, phi in jobs)
    dy, dx = np.mgrid[-h:h+1, -h:h+1]
    worst = 0.0
    per_radius = 1+len(angles_deg)
    for i in range(len(radii_ld)):
        ref = results[i*per_radius][1] # source at (r, 0), on a larger window for the rotation
        for j, phi in enumerate(np.radians(angles_deg)):
            # Offset (dx, dy) from the rotated source maps back to (u, v) from the reference source
            u = dx*np.cos(phi)+dy*np.sin(phi)
            v = -dx*np.sin(phi)+dy*np.cos(phi)
            rotated = map_coordinates(ref, [2*h+v, 2*h+u], order=5)
            exact = results[i*per_radius+1+j][1][h:3*h+1, h:3*h+1]
            worst = max(worst, np.abs(rotated-exact).sum()/exact.sum())
    return float(worst)


def psf_symmetry(optics, setting="auto", tol=0.01, num_cores=-1):
    """
    Whether the PSFs are circularly symmetric, and why.

    Args:
        optics (Optics): The optical system and image sampling.
        setting (bool or str): True/False from the user, or "auto" (or None): the optics
            model's declaration (backend attribute rotationally_symmetric), or, if the
            model does not declare it, a measurement (measure_psf_symmetry).
        tol (float): Largest measured relative difference still counted as symmetric.
        num_cores (int): Number of CPU processes (measurement only).

    Returns:
        tuple:
            symmetric (bool): Whether the PSFs are treated as circularly symmetric.
            reason (str): One line saying why (for the method message).
    """
    declared = getattr(optics.backend(), "rotationally_symmetric", None)
    if setting is not None and setting != "auto":
        note = "" if declared is None or declared == bool(setting) else \
            " (the %s model says %s)" % (optics.kind, "symmetric" if declared else "not symmetric")
        return bool(setting), "PSFs %s: set in the parameters%s" % ("symmetric" if setting else "not symmetric", note)
    if declared is not None:
        return bool(declared), "PSFs %s: %s model" % ("symmetric" if declared else "not symmetric", optics.kind)
    diff = measure_psf_symmetry(optics, num_cores=num_cores)
    return diff <= tol, "PSFs %s: measured %.2g%% difference to rotated copies (tol %.2g%%)" % (
        "symmetric" if diff <= tol else "not symmetric", 100*diff, 100*tol)


def image_symmetry(pre_img, setting="auto", tol=0.01):
    """
    Whether the target is circularly symmetric about the star, and why.

    Args:
        pre_img (np.ndarray): (n, n) target image, star at pixel (n-1)/2.
        setting (bool or str): True/False from the user, or "auto" (or None): measured
            with image_asymmetry.
        tol (float): Largest asymmetry and corner fraction still counted as symmetric.

    Returns:
        tuple:
            symmetric (bool): Whether the target is treated as circularly symmetric.
            reason (str): One line saying why (for the method message).
    """
    if setting is not None and setting != "auto":
        return bool(setting), "target %s: set in the parameters" % ("symmetric" if setting else "not symmetric")
    asym, corners = image_asymmetry(pre_img)
    symmetric = asym <= tol and corners <= tol
    pct = lambda v: "%.2g%%" % (100*v) if v >= 1e-8 else "0%" # (round-off level shown as 0)
    detail = "asymmetry %s, light in the corners %s (tol %.2g%%)" % (pct(asym), pct(corners), 100*tol)
    return symmetric, "target %s: %s" % ("symmetric" if symmetric else "not symmetric", detail)


def choose_method(method="auto", spatial_mode="auto", image_sym=True, psf_sym=True, rotate_allowed=True,
                  reasons=()):
    """
    Resolve the method from the parameter settings and the symmetries.

    method "auto": circular if the target and the PSFs are both symmetric, else spatial.
    spatial_mode "auto": rotate if the PSFs are symmetric and the optics model allows it
    (rotationally symmetric with a known band limit), else direct. Explicit settings are
    kept, with a warning when they contradict the symmetries; 'rotate' for optics that
    cannot use it is an error.

    Args:
        method (str): "auto", "circular" or "spatial".
        spatial_mode (str): "auto", "rotate" or "direct".
        image_sym (bool): Whether the target is circularly symmetric.
        psf_sym (bool): Whether the PSFs are circularly symmetric.
        rotate_allowed (bool): Whether the optics model supports mode 'rotate'.
        reasons (tuple): Strings explaining image_sym and psf_sym (for the message).

    Returns:
        tuple:
            method (str): "circular" or "spatial".
            spatial_mode (str or None): "rotate" or "direct"; None for circular.
            messages (list): Message lines: "Method: ..." first, then any WARNING lines.

    Raises:
        ValueError: For an unknown method or spatial_mode, or 'rotate' when not rotate_allowed.
    """
    messages = []
    if method == "auto":
        method = "circular" if (image_sym and psf_sym) else "spatial"
        how = "auto"
    elif method in ("circular", "spatial"):
        how = "set in the parameters"
        if method == "circular" and not (image_sym and psf_sym):
            messages.append("WARNING: the circular method is only exact for a symmetric target and symmetric PSFs")
    else:
        raise ValueError("method must be 'auto', 'circular' or 'spatial'")
    if method == "spatial":
        if spatial_mode == "auto":
            spatial_mode = "rotate" if (psf_sym and rotate_allowed) else "direct"
        elif spatial_mode == "rotate":
            if not rotate_allowed:
                raise ValueError("spatial_mode 'rotate' needs optics that are rotationally symmetric with a known band limit")
            if not psf_sym:
                messages.append("WARNING: spatial_mode 'rotate' assumes symmetric PSFs")
        elif spatial_mode != "direct":
            raise ValueError("spatial_mode must be 'auto', 'rotate' or 'direct'")
    else:
        spatial_mode = None
    messages.insert(0, "Method: %s%s (%s; %s)" % (method, "" if spatial_mode is None else ", " + spatial_mode + " mode",
                                                    how, "; ".join(reasons)))
    return method, spatial_mode, messages
