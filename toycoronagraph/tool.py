"""
Small helpers shared by the package: coordinate conversions, input checks and frame interpolation.

Image convention (as in toycoronagraph.optics): the star sits at pixel coordinate
c = (img_pixel-1)/2 on both axes; sky x increases with the column and sky y with the row
(the FITS convention, as toycoronagraph.scene), so a source at polar position (r, theta)
lands at row c + r*sin(theta), column c + r*cos(theta) (r in pixels).
"""
import numbers
import numpy as np
from skimage.transform import warp, SimilarityTransform
from scipy.spatial import distance

def convert_to_polar(pos_car):
    """
    Converts cartesian coordinates to polar coordinates.

    Args:
        pos_car (array-like): (x, y) of the point, in any length unit (e.g. arcsec).

    Returns:
        numpy.ndarray: (2,) float array [r, theta]: the radius (same unit as the input) and the
            angle (radian, from +x towards +y, in (-pi, pi]).
    """
    r = np.sqrt(pos_car[0]**2 + pos_car[1]**2)
    theta = np.arctan2(pos_car[1], pos_car[0])
    return np.array([r, theta])

def arcsec_to_pixel(pos_polar, psf_scale, img_pixel):
    """
    Converts a polar sky position to a (row, col) pixel position in the image.

    The star is at pixel coordinate (img_pixel-1)/2 on both axes (the origin of
    the focal grid of optics.VortexBackend, and the
    rotation centre used by cir_psf_planets); x increases with the column and y
    with the row (FITS), the same orientation as cir_psf_planets, Target and scene.Grid.pixel.

    Args:
        pos_polar (array-like): The radius (arcsec) and angle (radian, from +x towards +y) of the point.
        psf_scale (float): The scale of the PSF, in arcsec per pixel.
        img_pixel (int): The number of pixels on one edge of the image.

    Returns:
        numpy.ndarray: (2,) float array, the (row, col) position in pixels (need not be integers).
    """
    c = (img_pixel-1)/2
    r = pos_polar[0]/psf_scale
    return np.array([c + r*np.sin(pos_polar[1]), c + r*np.cos(pos_polar[1])])

def is_positive_even_integer(number):
    """
    Checks if the input number is a positive even integer.

    Only Python int and np.int64 are accepted (bools count as int: True is odd, False is not positive).

    Args:
        number (int): The number to check.

    Returns:
        bool: True if the number is a positive even integer, False otherwise.
    """
    if not isinstance(number, (int, np.int64)) or number <= 0:
        return False
    return bool(number % 2 == 0)

def _is_number(value):
    """
    Whether a value is a finite real number.

    Args:
        value: Any object.

    Returns:
        bool: True for a finite Python or numpy int or float (any precision); False for a bool, NaN,
            infinity or anything else.
    """
    return isinstance(value, numbers.Real) and not isinstance(value, bool) and bool(np.isfinite(value))

def is_planet_pos_allowed(pos, mode):
    """
    Checks if the input planet position is valid based on the specified mode.

    Each parameter must be a finite real number (Python or numpy int or float, not a bool).

    Args:
        pos (list or numpy.ndarray): The planet position:
            - "moving": [a, e, pa, inc, t] = semi-major axis (arcsec, > 0), eccentricity (0 <= e < 1),
              position angle (degree), inclination (degree) and time (in units of the orbital period),
              see planet.planet_position;
            - "cartesian": [x, y] in arcsec;
            - "polar": [r, theta] = radius (arcsec, > 0) and angle (radian).
            Anything that is not a 1-D sequence (None, a single number, ...) is invalid.
        mode (str): "moving", "cartesian" or "polar".

    Returns:
        bool or None: True if the planet position is valid, False otherwise; None (after printing
            a message) for an unknown mode.
    """
    sizes = {"moving": 5, "cartesian": 2, "polar": 2}
    if mode not in sizes:
        print("unknown add planet mode")
        return None
    if np.ndim(pos) != 1 or len(pos) != sizes[mode]:
        if mode == "moving":
            print("The motion of a planet is described by five parameters: the semi-major axis, the eccentricity, "
                  "the position angle, the inclination, and the time in units of the period.")
        return False
    if not all(_is_number(p) for p in pos):
        return False
    if mode == "moving":
        return pos[0] > 0 and 0 <= pos[1] < 1
    if mode == "polar":
        return pos[0] > 0
    return True

def frame_interpolation(positions, images, pos):
    """
    Estimate the image of a source at pos from images of the same source at other positions.

    Each image is shifted (bilinear warp) so that its source moves from its own position to
    pos, and the shifted images are averaged with inverse-distance weights (1/|position - pos|).
    If pos coincides with one of the positions, that image is returned unchanged.

    Args:
        positions (array-like): Source position of each image, in pixels: (n,) scalars for
            shifts along x (columns) only, or (n, 2) (x, y) = (column, row) offsets.
        images (sequence of numpy.ndarray): n images of the same shape, e.g. (S, S) float.
        pos (float or array-like): The target position, in the same form as one entry of positions.

    Returns:
        numpy.ndarray: The interpolated image, same shape as each input image (float64).
    """
    pos = np.array(pos).reshape(1, -1)
    positions = np.array(positions).reshape(-1, np.size(pos))
    distances = distance.cdist(positions, pos, 'euclidean').flatten()
    # Translation (x, y) of each image; scalar positions shift along x only
    rel_pos = positions-pos
    if np.size(pos)==1:
        rel_pos = np.concatenate((rel_pos, np.zeros(rel_pos.shape)), axis=1)
    total_weight = 0
    final_image = np.zeros(np.shape(images[0]))
    for i, d in enumerate(distances):
        if d==0:
            return images[i]
        weight = 1/d
        final_image = final_image + weight*warp(images[i], SimilarityTransform(translation=rel_pos[i]))
        total_weight = total_weight + weight
    return final_image/total_weight
