"""
Stage 4 unit tests: the image functions of toycoronagraph/psf.py (radial PSF series, circular
method, planets, IWA, contrast, core throughput, kernel size) and of toycoronagraph/spatial.py
(blended convolution, spatial_psf, spatial_iwa, image_error_estimate), on tiny images.

Images are N = 32 (even) and 33 (odd) pixels at 4 px per lambda/D (psf_range = N/8), a charge-2
vortex. The slow brute-force comparisons and figures stay in tests/stage4_images.py.
"""
import io
import contextlib
from functools import lru_cache
import numpy as np
import pytest
import matplotlib.pyplot as plt
import toycoronagraph.psf as P
from toycoronagraph.optics import Optics
from toycoronagraph.scene import Grid
from toycoronagraph.spatial import (convolve_with_blended_psfs, equal_error_grid, spatial_psf, spatial_iwa,
                                    image_error_estimate)
from toycoronagraph.tool import arcsec_to_pixel

CHARGE, LYOT, NCORES = 2, 0.95, 2
SIZES = [32, 33]
SIZE_IDS = ["even", "odd"]
SCALE = 0.01  # arcsec per pixel, for the functions taking positions in arcsec
KERNEL = 31   # spatial kernel width used by the shared model (about +-3.75 lambda/D)

pytestmark = pytest.mark.filterwarnings("ignore:FigureCanvasAgg is non-interactive")
sizes = pytest.mark.parametrize("n", SIZES, ids=SIZE_IDS)


def optics(n):
    """Charge-2 vortex on an n x n image at 4 px per lambda/D."""
    return Optics.vortex(CHARGE, n, n/8, LYOT)


def radii(n):
    """(n, n) distance of every pixel from the star at (n-1)/2, in pixels."""
    c = (n-1)/2
    rows, cols = np.indices((n, n))
    return np.hypot(rows-c, cols-c)


def ring(n):
    """Circularly symmetric ring (peak 1 at 0.3 n pixels from the star)."""
    return np.exp(-0.5*((radii(n)-0.3*n)/(0.06*n))**2)


def arc(n):
    """Non-symmetric target: the right half of the ring plus an off-centre clump."""
    c = (n-1)/2
    rows, cols = np.indices((n, n))
    return ring(n)*(cols > c)+0.5*np.exp(-0.5*((rows-c-0.2*n)**2+(cols-c+0.15*n)**2)/(0.05*n)**2)


@lru_cache(maxsize=None)
def radial(n):
    """Radial PSF series on image-size frames (in memory, computed once per size)."""
    return P.radial_psf_series(optics(n), NCORES)


@lru_cache(maxsize=None)
def radial_ext(n):
    """Radial PSF series on extended frames (in memory, computed once per size)."""
    return P.radial_psf_series(optics(n), NCORES, P.extended_frame_size(n))


@lru_cache(maxsize=None)
def model(n):
    """Spatial PSF model ('rotate', K=KERNEL, coarse grid), built once per size; nothing is written."""
    with contextlib.redirect_stdout(io.StringIO()):
        return P.spatial_psf_calculation(optics(n), NCORES, radial(n), kernel_size=KERNEL, tol=0.05, n_validate=2)


def save_psfs(folder, n, extended=True):
    """Save the radial PSFs of size n to folder/psfs.npy (the cir_* functions read a file); returns the path."""
    path = str(folder/"psfs.npy")
    np.save(path, radial_ext(n) if extended else radial(n))
    return path


def peak(img):
    """(row, col) of the brightest pixel."""
    return np.array(np.unravel_index(np.argmax(img), img.shape))


# ---------------------------------------------------------------- radial PSF series

@sizes
def test_radial_source_positions_on_pixels_for_odd_between_for_even(n):
    """Odd N: sources at 0, 1, 2, ... px; even N: at sqrt(i^2+0.25) px (the star is between pixels)."""
    x = P.radial_source_positions(n)
    i = np.arange((n+1)//2)
    assert x.shape == ((n+1)//2,) and x.dtype == np.float64
    assert np.allclose(x, i if n % 2 else np.sqrt(i**2+0.25))
    assert np.all(np.diff(x) > 0)


@pytest.mark.parametrize("n", [16, 17, 32, 33, 95, 96])
def test_extended_frame_size_keeps_parity_and_covers_corners(n):
    """The extended frame has the image's parity (same lattice) and reaches sqrt(2) N + 8 pixels."""
    size = P.extended_frame_size(n)
    assert isinstance(size, int)
    assert (size-n) % 2 == 0 and size >= np.sqrt(2)*n+8


@sizes
def test_radial_psf_series_shapes_and_dtypes(n):
    """Image-size frames are float64 (N, N); extended frames are float32 (S, S); one per source radius."""
    s = P.extended_frame_size(n)
    assert radial(n).shape == ((n+1)//2, n, n) and radial(n).dtype == np.float64
    assert radial_ext(n).shape == ((n+1)//2, s, s) and radial_ext(n).dtype == np.float32
    assert np.all(radial(n) >= 0)


def test_radial_psf_series_refuses_frame_of_other_parity():
    """A frame size of the other parity would not share the image lattice: ValueError naming the parity."""
    with pytest.raises(ValueError, match="parity"):
        P.radial_psf_series(optics(32), NCORES, frame_size=51)


@sizes
def test_extended_frames_centre_is_the_image_frame(n):
    """The central N x N block of an extended frame is the image-size frame (same lattice and centre)."""
    ext = P._crop_centre(radial_ext(n), n)
    assert np.allclose(ext, radial(n), rtol=1e-4, atol=1e-6*radial(n).max())


@sizes
def test_radial_psf_peaks_on_its_source_along_plus_x(n):
    """Well outside the IWA, radial PSF i peaks on the star's row at column c + x_i (star at (N-1)/2)."""
    c = (n-1)/2
    i = (n+1)//2-4
    x = P.radial_source_positions(n)[i]
    row, col = peak(radial(n)[i])
    assert abs(row-c) <= 0.5 and abs(col-(c+x)) <= 1


def test_crop_centre_returns_central_block():
    """_crop_centre keeps the central img_pixel block of every frame; same-size frames come back unchanged."""
    frames = np.arange(2*7*7).reshape(2, 7, 7)
    assert np.array_equal(P._crop_centre(frames, 3), frames[:, 2:5, 2:5])
    same = P._crop_centre(frames, 7)
    assert np.array_equal(same, frames) and np.shares_memory(same, frames)


# ---------------------------------------------------------------- core throughput and IWA

def test_core_throughput_sums_each_psf_above_half_its_peak():
    """The core throughput of each PSF is the sum of its pixels above half its own peak."""
    psfs = np.random.default_rng(0).random((5, 9, 9))
    expected = [p[p > 0.5*p.max()].sum() for p in psfs]
    assert np.allclose(P.core_throughput_fun(psfs), expected)


def test_cir_iwa_is_first_index_reaching_half_the_largest_core_throughput():
    """cir_iwa returns the first radial index whose core throughput reaches half of the maximum."""
    values = [0.05, 0.1, 0.3, 0.49, 0.5, 0.8, 1.0]
    psfs = np.zeros((len(values), 5, 5))
    psfs[:, 2, 2] = values  # single-pixel PSFs: core throughput = the value
    assert P.cir_iwa(psfs) == 4


@sizes
def test_spatial_and_circular_iwa_agree(n):
    """The spatial model's IWA (pixels) matches the radial series' IWA within 1.5 pixels."""
    assert abs(spatial_iwa(model(n))-P.cir_iwa(radial(n))) <= 1.5


def test_cir_core_throughput_plot_saves_the_curve_with_the_iwa_marked(in_tmp):
    """Core_throughput.png is saved in the working folder; the curve is core_throughput_fun, the line cir_iwa."""
    psfs = radial(32)
    P.cir_core_throughput_plot(psfs, plot_dpi=20)
    assert (in_tmp/"Core_throughput.png").read_bytes().startswith(b"\x89PNG")
    curve, iwa_line = plt.gcf().axes[0].lines
    assert np.allclose(curve.get_ydata(), P.core_throughput_fun(psfs))
    assert np.allclose(iwa_line.get_xdata(), P.cir_iwa(psfs))
    plt.close("all")


# ---------------------------------------------------------------- circular method

def one_hot_psfs(n):
    """(m, S, S) radial PSFs where PSF i is 1 at pixel (0, i): the chunk's first row holds the weights."""
    m, s = (n+1)//2, P.extended_frame_size(n)
    psfs = np.zeros((m, s, s))
    psfs[np.arange(m), 0, np.arange(m)] = 1
    return psfs


def split_weights(img, n, nodes):
    """Each pixel's light split linearly in radius between the two source radii around its own (first `nodes`)."""
    x = P.radial_source_positions(n)
    w = np.zeros(len(x))
    for rk, v in zip(radii(n).ravel(), img.ravel()):
        i = min(max(np.searchsorted(x, rk, side='right')-1, 0), len(x)-2)
        f = min(max((rk-x[i])/(x[i+1]-x[i]), 0), 1)
        w[i] += v*(1-f)
        w[i+1] += v*f
    return w[:nodes]


@pytest.mark.parametrize("n", [64, 65], ids=SIZE_IDS)
def test_cir_chunks_weights(n):
    """Beyond the first _CIR_PIXEL_NODES source radii: profile x annulus area 2 pi x_i (x_i - x_(i-1));
    the first ones: the pixels' light split linearly in radius; chunk_iwa drops i < iwa."""
    from scipy.ndimage import map_coordinates
    k, iwa = P._CIR_PIXEL_NODES, 3
    img = np.exp(-0.5*((radii(n)-0.25*n)/(0.06*n))**2) + np.random.default_rng(1).random((n, n))
    x = P.radial_source_positions(n)
    c = (n-1)/2
    continuum = map_coordinates(img, [c+x, np.full(len(x), c)], order=3, mode='nearest')*2*np.pi*x*np.diff(x, prepend=0)
    chunk, chunk_iwa = P._cir_chunks(img, n, one_hot_psfs(n), iwa)
    assert chunk.dtype == np.float64
    assert np.allclose(chunk[0, k:len(x)], continuum[k:], rtol=1e-12)
    assert np.allclose(chunk[0, :k], split_weights(img, n, k), rtol=1e-12)
    assert np.allclose(chunk_iwa[0, :len(x)], np.where(np.arange(len(x)) >= iwa, chunk[0, :len(x)], 0), rtol=1e-12)


@pytest.mark.parametrize("n", [16, 17], ids=SIZE_IDS)
def test_cir_chunks_keep_the_light_of_the_star_pixel(n):
    """Light on the star's pixel(s) goes to the first radial PSF (regression: it had weight 0)."""
    centre = (radii(n) < 0.75).astype(float)  # the star's pixel (odd n) or the four central ones (even n)
    chunk, _ = P._cir_chunks(centre, n, one_hot_psfs(n), 0)
    assert chunk[0, :P._CIR_PIXEL_NODES].sum() == pytest.approx(centre.sum(), rel=1e-12)
    if n % 2:  # the star sits on a pixel: all its light on the first radial PSF (the source on the star)
        assert chunk[0, 0] == pytest.approx(1, rel=1e-12)


def test_cir_psf_passes_an_on_axis_point_without_a_coronagraph(in_tmp):
    """Through PSFs that let an on-axis source through (here: a delta at the star), a point on the star's
    pixel keeps its light (regression: the circular method dropped the star's pixel)."""
    n = 17
    s = P.extended_frame_size(n)
    psfs = np.zeros(((n+1)//2, s, s))
    for i in range(len(psfs)):
        psfs[i, s//2, s//2+i] = 1  # each source imaged on its own pixel
    np.save("delta_psfs.npy", psfs)
    point = np.zeros((n, n)); point[n//2, n//2] = 1
    final = P.cir_psf(point, [], [], 0.01, False, False, n, n, "delta_psfs.npy")
    # (the cubic rotations of a one-pixel delta spread it slightly: its sum is kept to ~1%)
    assert final.sum() == pytest.approx(1, abs=0.02) and np.argmax(final) == point.argmax()


@pytest.mark.parametrize("size, n", [(9, 5), (10, 6), (11, 11)])
def test_cir_rotate_sum_averages_the_rotations_about_the_star_and_crops(size, n):
    """rot_number = 4 averages the 0/90/180/270 deg turns about the frame centre; 1 is the chunk itself."""
    chunk = np.random.default_rng(1).random((size, size))
    turns = np.mean([np.rot90(chunk, k) for k in range(4)], axis=0)
    out = P._cir_rotate_sum(chunk, n, 4)
    assert out.shape == (n, n) and out.dtype == np.float64
    assert np.allclose(out, P._crop_centre(turns, n), rtol=0, atol=1e-12)
    assert np.allclose(P._cir_rotate_sum(chunk, n, 1), P._crop_centre(chunk, n), rtol=0, atol=1e-12)


@sizes
def test_cir_psf_shape_dtype_and_empty_target(n, in_tmp):
    """cir_psf returns an (N, N) float64 image; an empty target gives an all-zero image."""
    path = save_psfs(in_tmp, n)
    img = P.cir_psf(np.zeros((n, n)), [], [], SCALE, False, False, n, n, path)
    assert img.shape == (n, n) and img.dtype == np.float64
    assert not img.any()


@sizes
def test_cir_psf_scales_with_the_target(n, in_tmp):
    """Scaling the target (also by a negative factor) scales the final image by the same factor."""
    path = save_psfs(in_tmp, n)
    fa = P.cir_psf(ring(n), [], [], SCALE, False, False, n, n, path)
    assert np.allclose(P.cir_psf(3*ring(n), [], [], SCALE, False, False, n, n, path), 3*fa)
    assert np.allclose(P.cir_psf(-ring(n), [], [], SCALE, False, False, n, n, path), -fa)


@sizes
def test_cir_psf_is_additive(n, in_tmp):
    """The coronagraph is linear: the final image of a sum of targets is the sum of their final images."""
    path = save_psfs(in_tmp, n)
    a, b = ring(n), np.exp(-0.5*((radii(n)-0.15*n)/(0.04*n))**2)
    fa, fb = (P.cir_psf(t, [], [], SCALE, False, False, n, n, path) for t in (a, b))
    fab = P.cir_psf(a+b, [], [], SCALE, False, False, n, n, path)
    assert np.abs(fab-fa-fb).max() <= 1e-6*fab.max()


@sizes
def test_cir_psf_of_symmetric_target_is_symmetric_about_the_star(n, in_tmp):
    """A symmetric ring gives an image mirror-symmetric about the star's row and ~invariant under 90 deg turns."""
    img = P.cir_psf(ring(n), [], [], SCALE, False, False, n, n, save_psfs(in_tmp, n))
    assert np.abs(img-np.flipud(img)).max() <= 1e-5*img.max()
    assert np.abs(img-np.rot90(img)).max() <= 1e-2*img.max()
    assert abs(radii(n)[tuple(peak(img))]-0.3*n) <= 2  # the brightest light stays on the ring


@sizes
def test_cir_psf_with_image_size_frames_matches_inside_and_darkens_the_corners(n, in_tmp):
    """Image-size radial frames give the extended-frame image inside the inscribed circle, less light beyond it."""
    small = P.cir_psf(ring(n), [], [], SCALE, False, False, n, n, save_psfs(in_tmp, n, extended=False))
    full = P.cir_psf(ring(n), [], [], SCALE, False, False, n, n, save_psfs(in_tmp, n))
    c = (n-1)/2
    inside, corners = radii(n) < c-1, radii(n) > c+0.5
    assert np.abs(small-full)[inside].max() <= 1e-6*full.max()
    assert small[corners].sum() < 0.6*full[corners].sum()


@sizes
def test_cir_psf_iwa_ignore_drops_the_radial_psfs_inside_the_iwa(n, in_tmp):
    """iwa_ignore equals the full image made with the radial PSFs inside the IWA set to zero."""
    psfs = radial_ext(n)
    iwa = P.cir_iwa(psfs)
    np.save(in_tmp/"inner_off.npy", np.concatenate([np.zeros_like(psfs[:iwa]), psfs[iwa:]]))
    ignored = P.cir_psf(ring(n), [], [], SCALE, True, False, n, n, save_psfs(in_tmp, n))
    expected = P.cir_psf(ring(n), [], [], SCALE, False, False, n, n, str(in_tmp/"inner_off.npy"))
    assert iwa > 0
    assert np.allclose(ignored, expected, rtol=1e-6, atol=1e-12)


@sizes
def test_cir_psf_with_planets_adds_cir_psf_planets(n, in_tmp):
    """cir_psf with add_planet adds exactly the image of cir_psf_planets to the disk image."""
    path = save_psfs(in_tmp, n)
    pos, b = [np.array([0.3*n*SCALE, np.radians(35)])], [2.0]
    with_p = P.cir_psf(ring(n), pos, b, SCALE, False, True, n, n, path)
    without = P.cir_psf(ring(n), pos, b, SCALE, False, False, n, n, path)
    assert np.allclose(with_p-without, P.cir_psf_planets(pos, b, SCALE, n, path))


def test_cir_planet_beyond_the_radial_psfs_is_skipped(in_tmp, capsys):
    """A planet beyond the last radial PSF: cir_psf skips it with a message, cir_psf_planets silently."""
    n = 33
    path = save_psfs(in_tmp, n)
    far = [np.array([0.6*n*SCALE, 0.3])]
    with_p = P.cir_psf(ring(n), far, [1.0], SCALE, False, True, n, n, path)
    assert "Planet #1 is outside the range of the plot" in capsys.readouterr().out
    assert np.array_equal(with_p, P.cir_psf(ring(n), [], [], SCALE, False, False, n, n, path))
    assert not P.cir_psf_planets(far, [1.0], SCALE, n, path).any()
    assert capsys.readouterr().out == ""


@sizes
@pytest.mark.parametrize("method", ["circular", "spatial"])
@pytest.mark.parametrize("angle_deg", [35, 120, -60, 180])
def test_planet_lands_where_arcsec_to_pixel_and_grid_pixel_put_it(n, method, angle_deg, in_tmp):
    """Regression (orientation): a planet at (r, angle) peaks at Grid.pixel / arcsec_to_pixel, rows = +y."""
    pos = np.array([0.3*n*SCALE, np.radians(angle_deg)])
    if method == "circular":
        img = P.cir_psf_planets([pos], [1.0], SCALE, n, save_psfs(in_tmp, n))
    else:
        img = P.spatial_psf_planets([pos], [1.0], SCALE, optics(n).settings())
    expected = arcsec_to_pixel(pos, SCALE, n)
    grid = Grid(size=n, pixel_scale=SCALE, wavelength=1e-6, diameter=2.4)
    assert np.allclose(expected, grid.pixel(pos[0], pos[1]))
    assert np.abs(peak(img)-expected).max() <= 1


@sizes
def test_interpolated_planet_keeps_the_exact_planet_flux(n, in_tmp):
    """The circular method's interpolated planet PSF has the exact PSF's light within 1%."""
    pos = np.array([0.3*n*SCALE, np.radians(35)])
    interpolated = P.cir_psf_planets([pos], [1.0], SCALE, n, save_psfs(in_tmp, n))
    exact = P.exact_point_sources(optics(n), [arcsec_to_pixel(pos, SCALE, n)], [1.0])
    assert abs(interpolated.sum()/exact.sum()-1) < 0.01


@sizes
def test_cir_psf_contrast_returns_planet_and_dust_light(n, in_tmp):
    """cir_psf_contrast gives planet > 0 and 0 <= dust (IWA ignored) <= dust; the planet part is linear."""
    path = save_psfs(in_tmp, n)
    r = 0.3*n
    planet_b, dust_b, dust_iwa = P.cir_psf_contrast(ring(n), r, 0.6, 1.0, SCALE, n, n, path)
    planet_2, dust_2, _ = P.cir_psf_contrast(ring(n), r, 0.6, 2.0, SCALE, n, n, path)
    assert all(np.ndim(v) == 0 for v in (planet_b, dust_b, dust_iwa))
    assert planet_b > 0 and 0 <= dust_iwa <= dust_b
    assert np.isclose(planet_2, 2*planet_b) and np.isclose(dust_2, dust_b)


# ---------------------------------------------------------------- exact point sources

@sizes
def test_exact_point_sources_peak_on_source_and_scale_with_value(n):
    """Each source is propagated at its own (row, col); the image is linear in the values."""
    o = optics(n)
    src = (0.2*n, 0.75*n)
    one = P.exact_point_sources(o, [src], [1.0])
    assert one.shape == (n, n) and one.dtype == np.float64
    assert np.abs(peak(one)-np.array(src)).max() <= 1
    assert np.allclose(P.exact_point_sources(o, [src, src], [1.0, 2.5]), 3.5*one)


def test_exact_point_sources_skips_sources_outside_the_frame(capsys):
    """Sources more than half a pixel outside the image are skipped with a message; -0.5 is still inside."""
    o = optics(16)
    img = P.exact_point_sources(o, [(16.2, 3.0), (-0.5, 3.0)], [1.0, 1.0])
    assert "Planet #1 is outside the range of the plot" in capsys.readouterr().out
    assert np.allclose(img, P.exact_point_sources(o, [(-0.5, 3.0)], [1.0]))
    assert img.any()


# ---------------------------------------------------------------- spatial method

@pytest.mark.parametrize("n", [16, 17], ids=SIZE_IDS)
@pytest.mark.parametrize("ncores", [1, 2])
def test_blended_convolution_matches_direct_sum(n, ncores):
    """Pure NumPy: each pixel adds im * tmap * sum_k w_k kernel_k centred on it (kernels uncropped)."""
    K = 5
    r = np.arange(0.25, 30, 0.5)
    g = equal_error_grid(n, r, 0.07/(1+(r/6)**3)+1e-3, 0.4*(r/5)/(1+(r/5)**2)+0.02, 0.05, 4.0, 8)
    rng = np.random.default_rng(0)
    psfs = rng.random((len(g['positions']), K, K))
    tmap = rng.random((n, n))
    im = rng.random((n, n))*(rng.random((n, n)) > 0.3)
    direct = np.zeros((n+K-1, n+K-1))
    for row in range(n):
        for col in range(n):
            kernel = np.einsum('k,kij->ij', g['wts'][row, col], psfs[g['inds'][row, col]])
            direct[row:row+K, col:col+K] += im[row, col]*tmap[row, col]*kernel
    direct = direct[K//2:K//2+n, K//2:K//2+n]
    out = convolve_with_blended_psfs(im, psfs, g['inds'], g['wts'], tmap, ncores=ncores)
    assert out.dtype == np.float64
    assert np.allclose(out, direct, rtol=1e-10, atol=1e-12*direct.max())
    assert not convolve_with_blended_psfs(np.zeros((n, n)), psfs, g['inds'], g['wts'], tmap).any()


@sizes
def test_spatial_model_contract(n):
    """The model holds float32 (n_samples, K, K) kernels, (N, N, 4) blend weights summing to 1, and tmap."""
    m = model(n)
    ns = len(m['positions'])
    assert m['mode'] == 'rotate' and m['kernel_size'] == KERNEL and m['img_pixel'] == n
    assert m['psfs'].shape == (ns, KERNEL, KERNEL) and m['psfs'].dtype == np.float32
    assert m['throughput'].shape == (ns,) and np.all((m['throughput'] >= 0) & (m['throughput'] <= 1+1e-9))
    assert m['inds'].shape == m['wts'].shape == (n, n, 4)
    assert np.allclose(m['wts'].sum(axis=-1), 1)
    assert m['tmap'].shape == (n, n) and m['tmap'].min() >= 0
    c = int((n-1)//2)
    assert m['tmap'][c, c] < 0.2 < m['tmap'][c, 0]  # the vortex blocks the star, transmits far off-axis


@sizes
def test_spatial_psf_shape_dtype_and_linearity(n):
    """spatial_psf returns an (N, N) float64 image, zero for an empty target, linear in the target."""
    m = model(n)
    a, b = ring(n), arc(n)
    fa = spatial_psf(a, m, report=False, ncores=NCORES)
    assert fa.shape == (n, n) and fa.dtype == np.float64
    assert not spatial_psf(np.zeros((n, n)), m, report=False).any()
    assert np.allclose(spatial_psf(2*a+b, m, report=False, ncores=NCORES),
                       2*fa+spatial_psf(b, m, report=False, ncores=NCORES))


@sizes
def test_spatial_psf_does_not_modify_the_target(n):
    """The target array passed to spatial_psf (also with iwa_ignore) is left unchanged."""
    target = ring(n)
    before = target.copy()
    spatial_psf(target, model(n), iwa_ignore=True, report=False)
    assert np.array_equal(target, before)


@sizes
def test_spatial_and_circular_methods_agree_on_a_symmetric_ring(n, in_tmp):
    """On a symmetric ring, the spatial image matches the circular one within 2% (relative L1)."""
    sp = spatial_psf(ring(n), model(n), report=False)
    cir = P.cir_psf(ring(n), [], [], SCALE, False, False, n, n, save_psfs(in_tmp, n))
    assert np.abs(sp-cir).sum()/cir.sum() < 0.02


@sizes
def test_spatial_psf_iwa_ignore_removes_target_inside_iwa(n):
    """With iwa_ignore, a target inside the IWA gives nothing and one outside the IWA is unchanged."""
    m = model(n)
    iwa = spatial_iwa(m)
    inner = (radii(n) < 0.5*iwa).astype(float)
    outer = ring(n)*(radii(n) >= iwa)
    assert inner.any() and not spatial_psf(inner, m, iwa_ignore=True, report=False).any()
    assert np.allclose(spatial_psf(outer, m, iwa_ignore=True, report=False),
                       spatial_psf(outer, m, iwa_ignore=False, report=False), rtol=1e-12, atol=0)


@sizes
def test_spatial_psf_planets_match_exact_point_sources(n):
    """spatial_psf_planets propagates each planet at arcsec_to_pixel through the model's optics."""
    pos = [np.array([0.3*n*SCALE, np.radians(35)]), np.array([0.2*n*SCALE, np.radians(-120)])]
    img = P.spatial_psf_planets(pos, [1.0, 0.5], SCALE, model(n))
    expected = P.exact_point_sources(optics(n), [arcsec_to_pixel(p, SCALE, n) for p in pos], [1.0, 0.5])
    assert np.allclose(img, expected)


@sizes
def test_spatial_psf_contrast_returns_planet_and_dust_light(n):
    """spatial_psf_contrast gives planet > 0 and 0 <= dust (IWA ignored) <= dust; the planet part is linear."""
    pos = np.array([0.3*n*SCALE, np.radians(35)])
    planet_b, dust_b, dust_iwa = P.spatial_psf_contrast(arc(n), pos, 1.0, SCALE, model(n), ncores=NCORES, report=False)
    planet_2, dust_2, _ = P.spatial_psf_contrast(arc(n), pos, 2.0, SCALE, model(n), ncores=NCORES, report=False)
    assert all(np.ndim(v) == 0 for v in (planet_b, dust_b, dust_iwa))
    assert planet_b > 0 and 0 <= dust_iwa <= dust_b
    assert np.isclose(planet_2, 2*planet_b) and np.isclose(dust_2, dust_b)


# ---------------------------------------------------------------- error report and estimate

@sizes
def test_error_report_structure(n):
    """The model's error report has 5 stages per radius band, the last one the cumulated error."""
    m = model(n)
    nb = len(m['err_edges'])-1
    assert m['err_max'].shape == m['err_median'].shape == (5, nb)
    assert list(m['err_stages'])[-1] == "cumulated (end to end)"
    ok = ~np.isnan(m['err_max'][-1])
    assert ok.any() and np.all(m['err_max'][:, ok] >= 0)


@sizes
def test_image_error_estimate_keys_and_bounds(n):
    """The estimate gives typical <= upper (fractions) and a typical error per stage name."""
    m = model(n)
    est = image_error_estimate(m, arc(n), verbose=False)
    assert set(est) == {"typical", "upper", "per_stage"}
    assert 0 <= est['typical'] <= est['upper'] < 1
    assert list(est['per_stage']) == [str(s) for s in m['err_stages']]


def test_image_error_estimate_of_empty_image_is_empty():
    """An image without flux has no estimate: an empty dict."""
    assert image_error_estimate(model(32), np.zeros((32, 32)), verbose=False) == {}


def test_spatial_psf_prints_the_estimate_only_when_asked(capsys):
    """report=True prints the estimated error of the image; report=False prints nothing."""
    m = model(32)
    spatial_psf(ring(32), m, report=False)
    assert capsys.readouterr().out == ""
    spatial_psf(ring(32), m, report=True)
    assert "Estimated error of this image" in capsys.readouterr().out


# ---------------------------------------------------------------- kernel size and arguments

def test_tol_for_error_target_is_monotonic_and_clipped():
    """The grid tolerance grows with the error target and stays within [0.003, 0.05]."""
    tols = [P.tol_for_error_target(t) for t in (1e-6, 0.005, 0.02, 0.05, 10.0)]
    assert tols[0] == 0.003 and tols[-1] == 0.05
    assert np.all(np.diff(tols) >= 0)


@sizes
def test_choose_kernel_size_is_smallest_odd_size_meeting_the_target(n):
    """K is odd, its crop error is within 80% of the target, and the next smaller size is not."""
    target = 0.05
    with contextlib.redirect_stdout(io.StringIO()):
        K, info = P.choose_kernel_size(optics(n), ring(n), target, NCORES)
    sizes_ = info['kernel_sizes']
    assert isinstance(K, int) and K % 2 == 1
    assert np.array_equal(sizes_, np.arange(1, 2*n, 2))
    assert np.all(np.diff(info['crop']) <= 1e-12) and abs(info['crop'][-1]) < 1e-12  # 2N-1 covers the frame
    i = int(np.flatnonzero(sizes_ == K)[0])
    assert info['crop'][i] <= 0.8*target and (i == 0 or info['crop'][i-1] > 0.8*target)


@sizes
def test_choose_kernel_size_from_radial_psfs_matches_propagated_probes(n, monkeypatch):
    """With the radial PSFs, probes inside the last radial PSF are rotated radial PSFs (not propagated).
    The kernel size is that of propagated probes within one step, and the propagated (true) crop error
    at it stays within the crop budget (an exact crop right at the threshold may tip it by ~1%)."""
    with contextlib.redirect_stdout(io.StringIO()):
        k_exact, exact = P.choose_kernel_size(optics(n), arc(n), 0.05, NCORES)
        calls = []
        probe = P._probe_worker
        monkeypatch.setattr(P, "_probe_worker", lambda *args: calls.append(args) or probe(*args))
        k_rot, rot = P.choose_kernel_size(optics(n), arc(n), 0.05, 1, radial_psfs=radial_ext(n))
    x_max = P.radial_source_positions(n)[-1]
    c = (n-1)/2
    assert all(np.hypot(row-c, col-c) >= x_max for row, col, _ in calls)  # only the corners are propagated
    assert abs(k_rot-k_exact) <= 2 and np.abs(rot['crop']-exact['crop']).max() < 2e-2
    assert exact['crop'][exact['kernel_sizes'] == k_rot][0] <= 1.05*0.8*0.05


def test_choose_kernel_size_without_target_weighs_pixels_equally():
    """pre_img None and an all-zero image both weigh every pixel equally (same choice)."""
    o = optics(16)
    with contextlib.redirect_stdout(io.StringIO()):
        k_none, info_none = P.choose_kernel_size(o, None, 0.05, NCORES, n_probe=8)
        k_zero, info_zero = P.choose_kernel_size(o, np.zeros((16, 16)), 0.05, NCORES, n_probe=8)
    assert k_none == k_zero and np.allclose(info_none['crop'], info_zero['crop'])


@pytest.mark.parametrize("kwargs, error, words", [
    (dict(mode="fast"), ValueError, "mode must be"),
    (dict(mode="rotate", kernel_size=16), ValueError, "odd"),
])
def test_spatial_psf_calculation_rejects_bad_arguments(kwargs, error, words):
    """An unknown mode or an even kernel size is refused before any PSF is computed."""
    with pytest.raises(error, match=words):
        P.spatial_psf_calculation(optics(32), NCORES, radial(32), tol=0.05, **kwargs)


def test_spatial_psf_calculation_refuses_rotate_for_generic_optics():
    """Optics without declared symmetry and band limit cannot use mode 'rotate' (ValueError says 'direct')."""
    generic = Optics('vortex_as_generic', 32, 4.0, optics(32).params, module='extra_backends')
    with pytest.raises(ValueError, match="rotate.*direct"):
        P.spatial_psf_calculation(generic, NCORES, mode='rotate', kernel_size=15)


def test_spatial_psf_calculation_refuses_kernel_cube_above_memory_limit():
    """A kernel cube larger than max_kernel_gb raises MemoryError naming max_kernel_gb."""
    with contextlib.redirect_stdout(io.StringIO()), pytest.raises(MemoryError, match="max_kernel_gb"):
        P.spatial_psf_calculation(optics(32), NCORES, radial(32), kernel_size=15, tol=0.05, max_kernel_gb=1e-9,
                                  validate=False)
