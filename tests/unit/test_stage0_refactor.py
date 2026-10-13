"""
Unit tests for Stage 0 (refactor guard): the radial PSF series of the circular method, and the
small helpers in toycoronagraph/tool.py.

Covers toycoronagraph.psf: radial_source_positions, extended_frame_size, radial_psf_series (image
and extended frames), _crop_centre, the vortex backend's source wavefront, core_throughput_fun / cir_iwa;
and toycoronagraph.tool: convert_to_polar, is_positive_even_integer, is_planet_pos_allowed,
frame_interpolation. (tool.arcsec_to_pixel is tested with the pixel grid in test_stage2_grid.py.)

Small images (N = 16 and 17, 8 px per lambda/D), serial hcipy runs; a few seconds in total.
The brute-force baseline comparison stays in tests/stage0_refactor.py.
"""
from functools import lru_cache
import numpy as np
import pytest
from toycoronagraph import psf as P
from toycoronagraph.optics import Optics
from toycoronagraph.tool import (convert_to_polar, is_positive_even_integer, is_planet_pos_allowed,
                                 frame_interpolation)

CHARGE, LYOT = 2, 0.95
SIZES = [16, 17]  # even (star between the four central pixels) and odd (star on a pixel)


def psf_range(n):
    """Half-width of an n-pixel image in lambda/D, for 8 px per lambda/D."""
    return n/16


def optics(n):
    """The charge-2 vortex used throughout this file, for image size n."""
    return Optics.vortex(CHARGE, n, psf_range(n), LYOT)


@lru_cache(maxsize=None)
def radial(n):
    """Radial PSF series on image-size frames for size n (computed once, serially; read-only)."""
    psfs = P.radial_psf_series(optics(n), 1)
    psfs.setflags(write=False)
    return psfs


# ---------------------------------------------------------------- radial_source_positions

@pytest.mark.parametrize("n", [1, 2, 16, 17, 64, 65])
def test_radial_source_positions_length_and_dtype(n):
    """One source per pixel of radius: (n+1)//2 float64 positions, strictly increasing."""
    x = P.radial_source_positions(n)
    assert x.shape == ((n+1)//2,) and x.dtype == np.float64
    assert np.all(np.diff(x) > 0)


def test_radial_source_positions_odd_n_on_pixels():
    """Odd N: the star is on a pixel and the sources sit 0, 1, 2, ... pixels from it."""
    assert np.array_equal(P.radial_source_positions(17), np.arange(9.))


def test_radial_source_positions_even_n_between_pixels():
    """Even N: source i sits sqrt(i^2 + 0.25) pixels from the star (the first one half a pixel away)."""
    x = P.radial_source_positions(16)
    assert x[0] == 0.5
    assert np.allclose(x, np.sqrt(np.arange(8)**2+0.25), rtol=0, atol=1e-15)


# ---------------------------------------------------------------- extended_frame_size and _crop_centre

@pytest.mark.parametrize("n", [1, 2, 15, 16, 17, 64, 65, 511, 512])
def test_extended_frame_size_parity_and_size(n):
    """ceil(sqrt(2) n) + 8, plus one if needed to share the image's parity (same lattice and centre)."""
    size = P.extended_frame_size(n)
    base = int(np.ceil(np.sqrt(2)*n))+8
    assert isinstance(size, int)
    assert size%2 == n%2
    assert size in (base, base+1)
    # Rotated about the star, the frame still covers the image corners
    assert (size-1)/2 >= np.sqrt(2)*(n-1)/2


@pytest.mark.parametrize("n, s", [(16, 26), (17, 27), (16, 16)])
def test_crop_centre_keeps_the_star_pixel(n, s):
    """_crop_centre returns a view whose centre (n-1)/2 is the larger frame's centre (s-1)/2."""
    big = np.zeros((2, s, s))
    cs, cn = (s-1)//2, (n-1)//2
    big[:, cs, cs] = 1.0           # a marker at (or next to, for even sizes) the frame centre
    out = P._crop_centre(big, n)
    assert out.shape == (2, n, n)
    assert np.shares_memory(out, big)
    assert out[0, cn, cn] == 1.0 and out.sum() == 2.0


# ---------------------------------------------------------------- radial_psf_series

@pytest.mark.parametrize("n", SIZES)
def test_radial_series_shape_and_dtype(n):
    """Image-size frames: ((N+1)//2, N, N) float64, non-negative."""
    psfs = radial(n)
    assert psfs.shape == ((n+1)//2, n, n) and psfs.dtype == np.float64
    assert psfs.min() >= 0


@pytest.mark.parametrize("n", SIZES)
def test_radial_series_frame_i_is_source_on_plus_x(n):
    """Frame i is the backend's frame of a source at radial_source_positions[i] pixels on +x (y = 0)."""
    o = optics(n)
    x = P.radial_source_positions(n)
    for i in [0, len(x)//2, len(x)-1]:
        expected = o.backend().source(x[i]*o.pixel_scale, 0).frame()
        assert np.array_equal(radial(n)[i], expected), i
    # The outermost source pulls the light to the +column side of the star
    last = radial(n)[-1]
    col_centroid = (last.sum(axis=0)*np.arange(n)).sum()/last.sum()
    assert col_centroid > o.centre+1


@pytest.mark.parametrize("n", SIZES)
def test_radial_series_symmetric_about_x_axis(n):
    """Sources on the x axis give frames mirror-symmetric in rows about the star (rows = y)."""
    psfs = radial(n)
    assert np.abs(psfs-psfs[:, ::-1, :]).max() < 1e-4*psfs.max()


def test_radial_series_rejects_frame_size_of_other_parity():
    """A frame_size of the other parity would not share the image lattice: ValueError naming parity."""
    with pytest.raises(ValueError, match="parity"):
        P.radial_psf_series(optics(16), 1, frame_size=25)


@pytest.mark.parametrize("n", SIZES)
def test_extended_frames_share_the_image_lattice(n):
    """Frames of extended_frame_size(N): float32, and their central N x N part is the image frame."""
    ext = P.radial_psf_series(optics(n), 1, P.extended_frame_size(n))
    s = P.extended_frame_size(n)
    assert ext.shape == ((n+1)//2, s, s) and ext.dtype == np.float32
    assert np.abs(P._crop_centre(ext, n)-radial(n)).max() < 1e-6*radial(n).max()


def test_source_wavefront_is_aperture_times_tilt():
    """A source at (x, y) lambda/D enters as the aperture times the phase ramp exp(2 pi i (x u + y v))."""
    b = optics(16).backend()
    flat = np.asarray(b.wavefront(0, 0).electric_field)
    tilted = np.asarray(b.wavefront(0.7, -1.3).electric_field)
    inside = np.abs(flat) > 0
    assert np.allclose(np.abs(tilted), np.abs(flat))
    ramp = np.exp(2j*np.pi*(np.asarray(b.pupil_grid.x)*0.7-np.asarray(b.pupil_grid.y)*1.3))
    assert np.allclose(tilted[inside]/flat[inside], ramp[inside])


# ---------------------------------------------------------------- core throughput and IWA (circular method)

def test_core_throughput_is_computed_per_psf():
    """core_throughput_fun gives each PSF of a stack its own threshold (half its own peak)."""
    psfs = np.random.default_rng(0).random((6, 9, 9))**4
    expected = [np.where(p > 0.5*p.max(), p, 0).sum() for p in psfs]
    assert np.allclose(P.core_throughput_fun(psfs), expected)


def test_core_throughput_sums_pixels_above_half_peak():
    """The core throughput is the sum of the pixels strictly above half the PSF's peak."""
    psf = np.array([[[4.0, 2.0, 2.1], [1.0, 0.0, 3.0]]])
    assert P.core_throughput_fun(psf)[0] == pytest.approx(4.0+2.1+3.0)


def test_cir_iwa_is_first_index_reaching_half_the_largest_core():
    """cir_iwa is the index of the first PSF whose core throughput reaches half the largest."""
    peaks = [0.1, 0.2, 0.45, 0.5, 0.9, 1.0]
    psfs = np.zeros((len(peaks), 5, 5))
    psfs[:, 2, 2] = peaks          # single-pixel PSFs: core throughput = peak
    assert P.cir_iwa(psfs) == 3


# ---------------------------------------------------------------- tool.convert_to_polar

@pytest.mark.parametrize("xy, r, theta", [((1.0, 0.0), 1.0, 0.0), ((0.0, 2.0), 2.0, np.pi/2),
                                           ((-3.0, 0.0), 3.0, np.pi), ((0.0, -1.0), 1.0, -np.pi/2),
                                           ((1.0, 1.0), np.sqrt(2), np.pi/4)])
def test_convert_to_polar(xy, r, theta):
    """(x, y) -> [r, theta], theta from +x towards +y, as a (2,) float array."""
    out = convert_to_polar(xy)
    assert out.shape == (2,) and out.dtype == np.float64
    assert np.allclose(out, [r, theta])


def test_convert_to_polar_accepts_numpy_input():
    """A numpy array (x, y) gives the same answer as a tuple."""
    assert np.allclose(convert_to_polar(np.array([3.0, -4.0])), [5.0, np.arctan2(-4, 3)])


# ---------------------------------------------------------------- tool.is_positive_even_integer

@pytest.mark.parametrize("number, expected", [
    (2, True), (4, True), (np.int64(6), True),
    (0, False), (-2, False), (3, False), (2.0, False), ("2", False), (None, False),
    (np.int32(2), False),          # only Python int and np.int64 are accepted
    (True, False), (False, False), # bools count as int: True is odd, False is not positive
])
def test_is_positive_even_integer(number, expected):
    """Only positive even Python ints / np.int64 pass; the answer is a plain bool."""
    out = is_positive_even_integer(number)
    assert type(out) is bool and out == expected


# ---------------------------------------------------------------- tool.is_planet_pos_allowed

@pytest.mark.parametrize("pos, mode, expected", [
    ([1.0, 0.5, 30, 45, 0.1], "moving", True),
    ([np.float64(1.0), 0.0, np.int64(30), 45.0, 2], "moving", True),
    ([0.0, 0.5, 30, 45, 0.1], "moving", False),   # a must be > 0
    ([1.0, -0.1, 30, 45, 0.1], "moving", False),  # e >= 0
    ([1.0, 1.0, 30, 45, 0.1], "moving", False),   # e < 1
    ([1.0, 0.5, "30", 45, 0.1], "moving", False), # numbers only
    ([1.0, 0.5, 30, 45], "moving", False),        # five parameters
    ([-0.3, 0.2], "cartesian", True),
    ([np.float32(0.3), np.int32(2)], "cartesian", True),  # numpy scalars of any precision
    ([np.nan, 0.2], "cartesian", False),          # finite numbers only
    ([True, 0.2], "cartesian", False),            # not a bool
    (None, "cartesian", False),                   # not a sequence
    (0.3, "polar", False),
    ([0.3, 0.2, 0.1], "cartesian", False),
    ([0.5, -2.0], "polar", True),
    ([0.0, 1.0], "polar", False),                 # r must be > 0
    ([0.5], "polar", False),
])
def test_is_planet_pos_allowed(pos, mode, expected):
    """Each mode checks its parameter count, types and ranges."""
    assert is_planet_pos_allowed(pos, mode) is expected


def test_is_planet_pos_allowed_explains_wrong_moving_length(capsys):
    """A moving planet with the wrong number of parameters prints what the five parameters are."""
    assert is_planet_pos_allowed([1.0, 0.5], "moving") is False
    assert "five parameters" in capsys.readouterr().out


def test_is_planet_pos_allowed_unknown_mode_returns_none(capsys):
    """An unknown mode prints a message and returns None."""
    assert is_planet_pos_allowed([1.0, 0.5], "spherical") is None
    assert "unknown add planet mode" in capsys.readouterr().out


# ---------------------------------------------------------------- tool.frame_interpolation

def blob(shape, row, col, amp=1.0, sigma=2.0):
    """A Gaussian of the given amplitude centred on (row, col) of an image of the given shape."""
    rows, cols = np.indices(shape)
    return amp*np.exp(-((rows-row)**2+(cols-col)**2)/(2*sigma**2))


def interior(img, margin=4):
    """img without a margin-pixel border (the warp shifts zeros in there, against the Gaussian's tails)."""
    return img[margin:-margin, margin:-margin]


def test_frame_interpolation_returns_image_at_known_position():
    """If pos is one of the positions, that image is returned unchanged."""
    images = [blob((20, 20), 10, 8), blob((20, 20), 10, 11)]
    out = frame_interpolation([3.0, 6.0], images, 6.0)
    assert np.array_equal(out, images[1])


def test_frame_interpolation_scalar_positions_shift_along_columns():
    """Scalar positions are x offsets: each image is shifted along the columns to pos (integer shifts are exact)."""
    shape, r0, c0 = (24, 32), 12, 10
    positions = [3.0, 7.0]
    images = [blob(shape, r0, c0+p) for p in positions]
    out = frame_interpolation(positions, images, 5.0)
    assert out.shape == shape and out.dtype == np.float64
    assert np.allclose(interior(out), interior(blob(shape, r0, c0+5)), rtol=0, atol=1e-12)


def test_frame_interpolation_2d_positions_are_column_row():
    """(n, 2) positions are (x, y) = (column, row) offsets."""
    shape, r0, c0 = (32, 32), 12, 12
    positions = [(0.0, 0.0), (4.0, 2.0)]
    images = [blob(shape, r0+y, c0+x) for x, y in positions]
    out = frame_interpolation(positions, images, (2.0, 1.0))
    assert np.allclose(interior(out), interior(blob(shape, r0+1, c0+2)), rtol=0, atol=1e-12)


def test_frame_interpolation_inverse_distance_weights():
    """Images are averaged with weights 1/|position - pos|."""
    shape, r0, c0 = (24, 32), 12, 10
    images = [blob(shape, r0, c0, amp=1.0), blob(shape, r0, c0+3, amp=4.0)]
    out = frame_interpolation([0.0, 3.0], images, 1.0)  # weights 1 and 1/2
    assert np.allclose(interior(out), interior(blob(shape, r0, c0+1, amp=(1.0+0.5*4.0)/1.5)), rtol=0, atol=1e-12)
