"""
Unit tests for Stage 1: the spatially-varying convolution engine in toycoronagraph/spatial.py.

Covers convolve_with_spatial_psfs (methods 'full' / 'localized' / 'auto', serial and threaded,
transmission map, dtypes, errors), _get_psf_patch_bounds, _choose_spatial_psf_convolution_strategy,
_fft_cost, psf_convolve_cpu / psf_convolve_gpu (without CuPy), and a cross-check of
convolve_with_blended_psfs against the nearest-sample engine for one-hot blend weights.

Pure NumPy/SciPy (no hcipy, no files), even and odd image sizes; about a second in total.
"""
import numpy as np
import pytest
from scipy.signal import fftconvolve, convolve2d
from toycoronagraph.spatial import (convolve_with_spatial_psfs, convolve_with_blended_psfs, psf_convolve_cpu,
                                    psf_convolve_gpu, _get_psf_patch_bounds,
                                    _choose_spatial_psf_convolution_strategy, _fft_cost)

SIZES = [24, 25]  # even and odd image sizes
K = 5             # kernel width (odd)
METHODS = ["full", "localized", "auto"]


def kernels(n_psf, shape=(K, K), seed=0):
    """n_psf distinct odd-sized kernels with a heavy top row (asymmetric, so a flip would show)."""
    psfs = np.random.default_rng(seed).random((n_psf,)+shape)
    psfs[:, 0, :] += 3.0
    return psfs


def reference(im, psfs, psf_inds):
    """Brute-force definition: the sum over PSF indices of full-frame 'same' convolutions of each region."""
    out = np.zeros(im.shape)
    for i in np.unique(psf_inds):
        out += fftconvolve(np.where(psf_inds == i, im, 0.), psfs[i], mode='same')
    return out


def regions(ny, nx):
    """Three irregular regions labelled 0, 1 and 3 (index 2 deliberately unused)."""
    rows, cols = np.indices((ny, nx))
    return np.where(rows < ny//3, 0, np.where(cols < 2*nx//3, 1, 3))


# ---------------------------------------------------------------- convolve_with_spatial_psfs

@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("method", METHODS)
def test_delta_reproduces_kernel_unflipped(n, method):
    """A unit delta at (r, c) gives the kernel itself centred on (r, c), not its flip."""
    psfs = kernels(1)
    im = np.zeros((n, n))
    r, c = 7, n-6
    im[r, c] = 1.0
    out = convolve_with_spatial_psfs(im, psfs, np.zeros((n, n), int), ncores=1, method=method)
    h = K//2
    assert np.allclose(out[r-h:r+h+1, c-h:c+h+1], psfs[0])
    assert np.isclose(out.sum(), psfs[0].sum())


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("ncores", [1, 2])
def test_single_psf_equals_plain_convolution(n, method, ncores):
    """With one PSF everywhere the engine is an ordinary 'same' convolution (checked against convolve2d)."""
    psfs = kernels(1)
    im = np.random.default_rng(1).random((n, n))
    out = convolve_with_spatial_psfs(im, psfs, np.zeros((n, n), int), ncores=ncores, method=method)
    assert np.allclose(out, convolve2d(im, psfs[0], mode='same'))


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("ncores", [1, 2])
def test_regions_match_brute_force_reference(n, method, ncores):
    """Several regions (one PSF index unused): every method equals the per-region sum of convolutions."""
    psfs = kernels(4)
    im = np.random.default_rng(2).random((n, n))
    psf_inds = regions(n, n)
    out = convolve_with_spatial_psfs(im, psfs, psf_inds, ncores=ncores, method=method)
    assert out.shape == (n, n)
    assert np.allclose(out, reference(im, psfs, psf_inds))


@pytest.mark.parametrize("method", METHODS)
def test_non_square_image_and_kernels(method):
    """(ny, nx) images with ny != nx and (ky, kx) kernels with ky != kx (both odd) are handled."""
    ny, nx = 19, 30
    psfs = kernels(4, shape=(3, 7))
    im = np.random.default_rng(3).random((ny, nx))
    psf_inds = regions(ny, nx)
    out = convolve_with_spatial_psfs(im, psfs, psf_inds, ncores=1, method=method)
    assert np.allclose(out, reference(im, psfs, psf_inds))


@pytest.mark.parametrize("n", SIZES)
def test_transmission_map_applied_first_and_input_untouched(n):
    """coron_tmap multiplies the image before convolution, and the caller's image is not modified."""
    psfs = kernels(2)
    im = np.random.default_rng(4).random((n, n))
    im_copy = im.copy()
    psf_inds = (np.indices((n, n))[1] >= n//2).astype(int)
    tmap = np.linspace(0, 1, n)[None, :]*np.ones((n, 1))
    for method in METHODS:
        out = convolve_with_spatial_psfs(im, psfs, psf_inds, coron_tmap=tmap, ncores=1, method=method)
        assert np.allclose(out, reference(im*tmap, psfs, psf_inds)), method
    assert np.array_equal(im, im_copy)


@pytest.mark.parametrize("dtype", [np.int64, np.int32, np.float32])
def test_integer_and_float32_images_give_float64(dtype):
    """Integer and float32 images are convolved in float64 (no truncation, no in-place cast)."""
    psfs = kernels(1)
    im = np.zeros((16, 16), dtype=dtype)
    im[8, 8] = 2
    out = convolve_with_spatial_psfs(im, psfs, np.zeros((16, 16), int), coron_tmap=np.full((16, 16), 0.5), ncores=1)
    assert out.dtype == np.float64
    assert np.allclose(out[6:11, 6:11], psfs[0])  # 2 * 0.5 * kernel
    assert im[8, 8] == 2


def test_method_is_case_insensitive():
    """'FULL' and 'Localized' are accepted like their lower-case names."""
    psfs = kernels(4)
    im = np.random.default_rng(5).random((20, 20))
    psf_inds = regions(20, 20)
    expected = reference(im, psfs, psf_inds)
    for method in ["FULL", "Localized", "AUTO"]:
        assert np.allclose(convolve_with_spatial_psfs(im, psfs, psf_inds, ncores=1, method=method), expected)


def test_unknown_method_raises_value_error():
    """Any other method name is a ValueError listing the valid ones."""
    with pytest.raises(ValueError, match="method must be one of"):
        convolve_with_spatial_psfs(np.zeros((8, 8)), kernels(1), np.zeros((8, 8), int), method='nope')


def test_precomputed_patch_bounds_are_used():
    """Pre-computed (used_inds, patch_bounds) are taken as given: here, deliberately too-small boxes
    make 'localized' drop light outside them, which shows they were not recomputed."""
    psfs = kernels(2)
    n = 20
    im = np.ones((n, n))
    psf_inds = (np.indices((n, n))[1] >= n//2).astype(int)
    used, bounds = _get_psf_patch_bounds(psf_inds, 2, (K//2+1, K//2+1))
    good = convolve_with_spatial_psfs(im, psfs, psf_inds, patch_bounds=bounds, used_inds=used, ncores=1,
                                      method='localized')
    assert np.allclose(good, reference(im, psfs, psf_inds))
    small = bounds.copy()
    small[:, 1] = small[:, 0]+3        # keep only the first four rows of every box
    cropped = convolve_with_spatial_psfs(im, psfs, psf_inds, patch_bounds=small, used_inds=used, ncores=1,
                                         method='localized')
    assert np.all(cropped[4:] == 0) and cropped[:4].sum() > 0


def test_zero_image_gives_zero_image():
    """An all-zero image convolves to an all-zero image (edge case: nothing to spread)."""
    out = convolve_with_spatial_psfs(np.zeros((12, 12)), kernels(4), regions(12, 12), ncores=1)
    assert out.shape == (12, 12) and not out.any()


def test_use_gpu_without_cupy_raises_import_error():
    """use_gpu=True without CuPy installed fails with an ImportError that names CuPy."""
    try:
        import cupy  # noqa: F401
        pytest.skip("CuPy is installed")
    except ImportError:
        pass
    with pytest.raises(ImportError, match="CuPy"):
        psf_convolve_gpu(np.zeros((8, 8)), kernels(1)[0])
    with pytest.raises(ImportError, match="CuPy"):
        convolve_with_spatial_psfs(np.zeros((8, 8)), kernels(1), np.zeros((8, 8), int), use_gpu=True)


# ---------------------------------------------------------------- building blocks

def test_psf_convolve_cpu_is_same_mode_convolution():
    """psf_convolve_cpu is a 'same'-size convolution centred on the kernel's middle pixel."""
    im = np.random.default_rng(6).random((11, 14))
    psf = kernels(1, shape=(3, 5))[0]
    out = psf_convolve_cpu(im, psf)
    assert out.shape == im.shape
    assert np.allclose(out, convolve2d(im, psf, mode='same'))


def test_patch_bounds_pad_and_clip():
    """Each used PSF gets its region's bounding box padded by the half-widths and clipped to the image."""
    psf_inds = np.zeros((10, 12), int)
    psf_inds[2:4, 5:8] = 3   # rows 2-3, cols 5-7
    psf_inds[9, 11] = 1      # corner pixel
    used, bounds = _get_psf_patch_bounds(psf_inds, 5, (2, 1))
    assert list(used) == [0, 1, 3]                      # increasing, unused 2 and 4 left out
    assert bounds.shape == (3, 4) and np.issubdtype(bounds.dtype, np.integer)
    b = dict(zip(used, bounds.tolist()))
    assert b[0] == [0, 9, 0, 11]          # whole frame (padded, then clipped)
    assert b[1] == [7, 9, 10, 11]         # [9-2, 9, 11-1, 11] clipped at the far edges
    assert b[3] == [0, 5, 4, 8]           # [2-2, 3+2, 5-1, 7+1]


def test_strategy_prefers_localized_for_small_regions_and_full_on_tie():
    """Many small regions -> 'localized'; one region covering the frame costs the same -> 'full'."""
    n = 32
    rows, cols = np.indices((n, n))
    _, bounds = _get_psf_patch_bounds((rows//8)*4+cols//8, 16, (K//2+1, K//2+1))
    assert _choose_spatial_psf_convolution_strategy((n, n), (K, K), bounds) == 'localized'
    _, bounds = _get_psf_patch_bounds(np.zeros((n, n), int), 1, (K//2+1, K//2+1))
    assert _choose_spatial_psf_convolution_strategy((n, n), (K, K), bounds) == 'full'


def test_fft_cost_formula_and_vectorization():
    """_fft_cost is n log2 n of the padded size n = (h+ky-1)(w+kx-1), elementwise for arrays."""
    assert _fft_cost(10, 20, 3, 5) == pytest.approx(12*24*np.log2(12*24))
    h, w = np.array([10, 4]), np.array([20, 4])
    assert np.allclose(_fft_cost(h, w, 3, 5), [_fft_cost(10, 20, 3, 5), _fft_cost(4, 4, 3, 5)])


@pytest.mark.parametrize("n", SIZES)
def test_blended_engine_with_one_hot_weights_matches_nearest_engine(n):
    """convolve_with_blended_psfs with weight 1 on one sample per pixel equals convolve_with_spatial_psfs."""
    psfs = kernels(4)
    im = np.random.default_rng(7).random((n, n))
    psf_inds = regions(n, n)
    inds = np.zeros((n, n, 4), int)
    wts = np.zeros((n, n, 4))
    inds[..., 0], wts[..., 0] = psf_inds, 1.0
    tmap = np.random.default_rng(8).random((n, n))
    blended = convolve_with_blended_psfs(im, psfs, inds, wts, coron_tmap=tmap, ncores=2)
    assert blended.dtype == np.float64
    assert np.allclose(blended, convolve_with_spatial_psfs(im, psfs, psf_inds, coron_tmap=tmap, ncores=1))
