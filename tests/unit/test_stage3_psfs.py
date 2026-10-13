"""
Unit tests for Stage 3: optical systems, PSF kernels and the spatial PSF model, plus the
orientation of planets placed through the PSFs.

Covers toycoronagraph/optics.py (Optics, register_backend / get_backend, VortexBackend,
PupilPlaneSource), and from toycoronagraph/psf.py: _source_at_pixel, _spatial_psf_worker,
coronagraph_throughput, tol_for_error_target, _band, _rotate_geometry, _same_settings,
load_spatial_psfs, axis_products (cache reuse), choose_kernel_size (inputs and table),
spatial_psf_calculation (argument checks, None/"auto" choices, output structure, a tiny 'rotate'
model and a tiny 'direct' model), exact_point_sources, spatial_psf_planets, cir_psf_planets
(orientation: +y -> +row), the no-main-guard fix; and spatial.bandlimited_weights.

Small images (N = 32 and 33 at 4 px per lambda/D; spatial PSF models at N = 24 and 25,
8 px per lambda/D, K = 9), serial hcipy runs (one subprocess uses 2 cores); well under a minute in total.
The accuracy checks at the package's 16 px per lambda/D stay in tests/stage3_psfs.py.
"""
import sys
import pickle
import shutil
import subprocess
import dataclasses
from functools import lru_cache
import numpy as np
import pytest
from toycoronagraph import psf as P
from toycoronagraph import optics as O
from toycoronagraph.optics import Optics, VortexBackend, PupilPlaneSource, register_backend
from toycoronagraph.spatial import positions_to_polar, blend_map, bandlimited_weights
from toycoronagraph.tool import arcsec_to_pixel

CHARGE, LYOT = 2, 0.95
SIZES = [32, 33]        # even (star between four pixels) and odd (star on a pixel)
MODEL_SIZES = [24, 25]  # the same, for the (slower) spatial PSF models
K, NCORES = 9, 1


def optics4(n):
    """Charge-2 vortex on an n-pixel image at 4 px per lambda/D (half-width n/8 lambda/D)."""
    return Optics.vortex(CHARGE, n, n/8, LYOT)


def optics8(n):
    """Charge-2 vortex on an n-pixel image at 8 px per lambda/D (half-width n/16 lambda/D)."""
    return Optics.vortex(CHARGE, n, n/16, LYOT)


@lru_cache(maxsize=None)
def radial(n, px_per_ld):
    """Radial PSF series of optics4 / optics8 (px_per_ld 4 or 8), computed once, serially (read-only)."""
    psfs = P.radial_psf_series(optics4(n) if px_per_ld == 4 else optics8(n), NCORES)
    psfs.setflags(write=False)
    return psfs


@pytest.fixture(scope="module")
def spatial_model(tmp_path_factory):
    """Builder of tiny 'rotate' spatial PSF models of optics8(n) (K = 9; mode and tol left to the package
    as None; error_target 0.05; no validation), each built once per module. The model and its x-axis
    cache are saved in a module-wide pytest temp folder.

    Returns:
        callable: get(n) -> (model dict, path of the saved model .npz, path of the x-axis cache .npz).
    """
    folder = tmp_path_factory.mktemp("spatial_models")
    built = {}

    def get(n):
        """The model of size n, its save path and its axis-cache path (built on first use)."""
        if n not in built:
            save, axis = str(folder/("model_N%d.npz" % n)), str(folder/("axis_N%d.npz" % n))
            built[n] = (P.spatial_psf_calculation(optics8(n), NCORES, radial(n, 8), kernel_size=K, error_target=0.05,
                                                  validate=False, axis_cache=axis, save_path=save), save, axis)
        return built[n]
    return get


def peak(img):
    """(row, col) of the brightest pixel of img."""
    return np.array(np.unravel_index(np.argmax(img), img.shape))


@pytest.fixture
def restore_backends():
    """Undo registrations made by a test (the backend registry is global)."""
    saved = dict(O._BACKENDS)
    yield
    O._BACKENDS.clear()
    O._BACKENDS.update(saved)


# ================================================================ optics.Optics

def test_vortex_constructor_types_and_sorted_params():
    """Optics.vortex casts its inputs and stores the params as sorted (name, value) pairs."""
    o = Optics.vortex(np.int64(6), np.int64(64), 4, 1)
    assert o.kind == 'vortex' and o.module == 'toycoronagraph.optics'
    assert type(o.img_pixel) is int and type(o.psf_range) is float
    assert o.params == (('charge', 6), ('lyot_mask_size', 1.0))
    assert type(o.param('charge')) is int and type(o.param('lyot_mask_size')) is float


def test_param_missing_name_raises_key_error():
    """param() of a name the backend does not have is a KeyError."""
    with pytest.raises(KeyError):
        optics4(32).param('diameter')


def test_optics_is_frozen_hashable_and_picklable():
    """Optics is immutable, equal descriptions hash alike, and it survives pickling (parallel workers)."""
    o = optics4(32)
    with pytest.raises(dataclasses.FrozenInstanceError):
        o.img_pixel = 64
    assert o == optics4(32) and hash(o) == hash(optics4(32)) and o != optics4(33)
    assert pickle.loads(pickle.dumps(o)) == o


@pytest.mark.parametrize("n", SIZES)
def test_pixel_scale_and_centre(n):
    """pixel_scale = 2 psf_range / N lambda/D per pixel; the star sits at (N-1)/2."""
    o = optics4(n)
    assert o.pixel_scale == pytest.approx(0.25)
    assert o.centre == (n-1)/2


def test_settings_round_trip_through_npz(in_tmp):
    """settings() holds plain values; from_settings rebuilds the same Optics, also from 0-d arrays of an .npz."""
    o = optics8(33)
    s = o.settings()
    assert s == dict(optics_kind='vortex', optics_module='toycoronagraph.optics', img_pixel=33,
                     psf_range=33/16, optics_charge=2, optics_lyot_mask_size=0.95)
    assert Optics.from_settings(s) == o
    np.savez("settings.npz", **s)
    with np.load("settings.npz") as f:
        assert Optics.from_settings({k: f[k] for k in f.files}) == o


def test_from_settings_ignores_other_keys_and_defaults_module():
    """Keys not describing the optics are ignored; a missing optics_module means toycoronagraph.optics."""
    s = dict(optics8(32).settings(), mode='rotate', kernel_size=9, psfs=np.zeros(3))
    del s['optics_module']
    assert Optics.from_settings(s) == optics8(32)


def test_backend_is_cached_per_optics():
    """Equal Optics share one backend object within a process."""
    assert optics4(32).backend() is optics4(32).backend()
    assert isinstance(optics4(32).backend(), VortexBackend)


def test_unknown_kind_raises_value_error():
    """A kind that no module registers is a ValueError naming it and the registered kinds."""
    with pytest.raises(ValueError, match="unknown optics kind 'no_such_system'.*vortex"):
        Optics('no_such_system', 16, 1.0).backend()


def test_register_backend_returns_class_and_registers_it(restore_backends):
    """register_backend stores the class under its kind and returns it unchanged."""
    class Dummy:
        def __init__(self, optics):
            self.optics = optics
    assert register_backend('unit_dummy')(Dummy) is Dummy
    o = Optics('unit_dummy', 8, 1.0)
    assert isinstance(o.backend(), Dummy) and o.backend().optics == o


def test_backend_module_is_imported_on_demand(in_tmp, monkeypatch, restore_backends):
    """An unregistered kind is found by importing Optics.module (as parallel workers do)."""
    (in_tmp/"unit_lazy_backend.py").write_text(
        "from toycoronagraph.optics import register_backend\n"
        "@register_backend('unit_lazy')\n"
        "class LazyBackend:\n"
        "    def __init__(self, optics):\n"
        "        self.optics = optics\n")
    monkeypatch.syspath_prepend(str(in_tmp))
    try:
        assert 'unit_lazy' not in O._BACKENDS
        backend = Optics('unit_lazy', 8, 1.0, module='unit_lazy_backend').backend()
        assert type(backend).__name__ == 'LazyBackend' and 'unit_lazy' in O._BACKENDS
    finally:
        sys.modules.pop('unit_lazy_backend', None)


# ================================================================ optics.VortexBackend and PupilPlaneSource

@pytest.mark.parametrize("n", SIZES)
def test_vortex_backend_declares_symmetry_grids_and_normal_factor(n):
    """Symmetric, band-limited at D/lambda; 2N pupil samples; N x N image; normal_factor = pi (px per lambda/D / 2)^2."""
    b = optics4(n).backend()
    assert b.rotationally_symmetric is True and b.bandlimit == 1.0
    assert b.pupil_grid.size == (2*n)**2
    assert b.prop(b.wavefront(0, 0)).electric_field.size == n*n
    assert b.normal_factor == pytest.approx(np.pi*(4/2)**2)


@pytest.mark.parametrize("n", SIZES)
def test_source_lands_on_its_pixel_rows_are_y(n):
    """A source at (x, y) = ((col-c) s, (row-c) s) peaks on pixel (row, col): rows are y, columns x."""
    o = optics4(n)
    for row, col in [(5, 26), (27, 6), (13, 24), (24, 13)]:
        src, (x, y) = P._source_at_pixel(o, row, col)
        assert (x, y) == pytest.approx(((col-o.centre)*o.pixel_scale, (row-o.centre)*o.pixel_scale))
        assert isinstance(src, PupilPlaneSource)
        assert tuple(peak(src.frame())) == (row, col)


@pytest.mark.parametrize("n", SIZES)
def test_focal_grid_is_point_symmetric_about_the_centre(n):
    """Regression (even N star centring): sources at (x, y) and (-x, -y) give images mirrored through (N-1)/2."""
    b = optics4(n).backend()
    a, m = b.source(1.0, -0.5).frame(), b.source(-1.0, 0.5).frame()
    assert np.abs(a-m[::-1, ::-1]).max() < 1e-3*a.max()


@pytest.mark.parametrize("n", SIZES)
def test_throughput_depends_only_on_radius(n):
    """The vortex is rotationally symmetric: same throughput at the same distance in any direction."""
    b = optics4(n).backend()
    t = [b.source(2.0, 0).throughput, b.source(0, 2.0).throughput, b.source(-np.sqrt(2), np.sqrt(2)).throughput]
    assert np.ptp(t) < 1e-3


@pytest.mark.parametrize("n", SIZES)
def test_frame_is_field_squared_and_the_centred_image_block(n):
    """frame() = |field()|^2 on the N x N image (complex field), and equals block((0, 0), (N, N))."""
    src = optics4(n).backend().source(1.5, 0.75)
    field, frame = src.field(), src.frame()
    assert field.shape == frame.shape == (n, n) and np.iscomplexobj(field) and frame.dtype == np.float64
    assert np.allclose(frame, np.abs(field)**2, rtol=1e-12, atol=0)
    assert np.allclose(src.block((0.0, 0.0), (n, n)), frame, rtol=1e-10, atol=1e-14*frame.max())


@pytest.mark.parametrize("n", SIZES)
def test_non_square_block_is_rows_by_columns(n):
    """block(center, (ny, nx)) has ny rows (y) and nx columns (x): a 3 x 5 block about pixel (10, 20)
    is that crop of the frame."""
    o = optics4(n)
    src = o.backend().source(1.5, 0.75)
    center = ((20-o.centre)*o.pixel_scale, (10-o.centre)*o.pixel_scale)
    blk = src.block(center, (3, 5))
    assert blk.shape == (3, 5)
    assert np.allclose(blk, src.frame()[9:12, 18:23], rtol=1e-10, atol=1e-14*blk.max())


@pytest.mark.parametrize("n", SIZES)
def test_total_is_the_whole_focal_plane_light(n):
    """total (Parseval) exceeds the light on the image and matches a 4N x 4N block around the source
    (kept within one alias period, 2N/1.5 lambda/D, of the 2N-sample pupil grid)."""
    o = optics4(n)
    src = o.backend().source(1.5, 0.75)
    assert src.frame().sum() < src.total
    assert src.block((1.5, 0.75), (4*n, 4*n)).sum() == pytest.approx(src.total, rel=0.01)


@pytest.mark.parametrize("n", SIZES)
def test_coronagraph_throughput_nulls_the_star_and_rises_off_axis(n):
    """Throughput (0-1) is ~0 on the star and grows with radius; same values as the sources' own."""
    o = optics4(n)
    radii = np.array([0.0, 4.0, 8.0, 12.0])  # 0-3 lambda/D
    T = P.coronagraph_throughput(o, NCORES, radii)
    assert T.shape == (4,) and T.dtype == np.float64
    assert T[0] < 0.01 and np.all(np.diff(T) > 0) and 0.7 < T[-1] <= 1
    assert T[2] == pytest.approx(o.backend().source(8*o.pixel_scale, 0).throughput, rel=1e-12)


@pytest.mark.parametrize("n", SIZES)
def test_spatial_psf_worker_kernel_is_crop_of_frame_and_centred(n):
    """On a pixel the kernel is the frame's crop around it; off the lattice it is centred on the source."""
    o = optics4(n)
    h = K//2
    T, kernel = P._spatial_psf_worker(10, 22, K, o)
    src, _ = P._source_at_pixel(o, 10, 22)
    assert kernel.shape == (K, K) and kernel.dtype == np.float64
    assert T == pytest.approx(src.throughput, rel=1e-12)
    frame = src.frame()
    assert np.allclose(kernel, frame[10-h:10+h+1, 22-h:22+h+1], rtol=1e-10, atol=1e-14*frame.max())
    _, off = P._spatial_psf_worker(10.3, 21.6, K, o)
    assert tuple(peak(off)) == (h, h)


# ================================================================ psf helpers

def test_tol_for_error_target_law_and_clipping():
    """tol = (blend_share * target / 1.5)^(1/1.6), clipped to [0.003, 0.05], increasing with the target."""
    assert P.tol_for_error_target(0.02) == pytest.approx((0.2*0.02/1.5)**(1/1.6))
    assert P.tol_for_error_target(0.02, blend_share=0.5) == pytest.approx((0.5*0.02/1.5)**(1/1.6))
    assert P.tol_for_error_target(1e-8) == 0.003 and P.tol_for_error_target(1.0) == 0.05
    assert isinstance(P.tol_for_error_target(0.02), float)
    targets = np.geomspace(1e-4, 0.5, 20)
    assert np.all(np.diff([P.tol_for_error_target(t) for t in targets]) >= 0)


def test_band_is_cycles_per_pixel_with_margin():
    """_band converts the backend's band limit (cycles per lambda/D) to cycles per pixel, times 1.05."""
    assert P._band(optics8(32)) == pytest.approx(1.05*1.0*optics8(32).pixel_scale)


@pytest.mark.parametrize("n", SIZES)
def test_rotate_geometry_window_fits_the_image_up_to_r_near(n):
    """c = (N-1)/2, h = K//2, wh = ceil(h sqrt 2) + 18; at r_near the window's last column is N-1."""
    c, h, wh, r_near = P._rotate_geometry(n, K)
    assert c == (n-1)/2 and h == K//2 and wh == int(np.ceil(h*np.sqrt(2)))+18
    last_col = lambda r: int(np.floor(c+r+0.5))+wh
    assert last_col(r_near) == n-1 and last_col(r_near+1) > n-1


def test_same_settings_strings_exact_numbers_close():
    """Strings must match exactly, numbers within np.isclose; a missing key is a mismatch."""
    cached = dict(optics_kind='vortex', psf_range=2.0, kernel_size=9)
    assert P._same_settings(cached, dict(optics_kind='vortex', psf_range=2.0+1e-12))
    assert not P._same_settings(cached, dict(optics_kind='vortex_as_generic'))
    assert not P._same_settings(cached, dict(kernel_size=11))
    assert not P._same_settings(cached, dict(tol=0.02))


def test_load_spatial_psfs_turns_0d_entries_into_python_values(in_tmp):
    """Scalars saved in an .npz come back as Python values; arrays stay arrays."""
    np.savez("m.npz", psfs=np.ones((2, 3, 3), np.float32), mode='rotate', kernel_size=9, tol=0.02)
    m = P.load_spatial_psfs("m.npz")
    assert m['mode'] == 'rotate' and isinstance(m['mode'], str)
    assert m['kernel_size'] == 9 and isinstance(m['kernel_size'], int)
    assert isinstance(m['tol'], float)
    assert m['psfs'].shape == (2, 3, 3) and m['psfs'].dtype == np.float32


# ================================================================ spatial.bandlimited_weights

def test_bandlimited_weights_on_a_node_is_that_node():
    """At a node (within 1e-9) the weights are exactly that node with weight 1."""
    idx, w = bandlimited_weights(np.arange(-5., 6.), 2.0, 0.06)
    assert list(idx) == [7] and list(w) == [1.0]
    idx, w = bandlimited_weights([0.0, 1.0, 2.0], 1.0+1e-12, 0.2)
    assert list(idx) == [1] and list(w) == [1.0]


@pytest.mark.parametrize("nodes", [np.arange(-20., 21.), np.concatenate([-np.sqrt(np.arange(20, 0, -1)**2+0.25),
                                                                         np.sqrt(np.arange(20)**2+0.25)])],
                         ids=["uniform odd-N nodes", "non-uniform even-N nodes"])
def test_bandlimited_weights_interpolate_band_limited_functions(nodes):
    """Any function with |f| <= band is interpolated exactly (sorted nearest n_taps nodes)."""
    rng = np.random.default_rng(1)
    band = 0.06
    fr, ph = rng.uniform(-band, band, 15), rng.uniform(0, 2*np.pi, 15)
    fn = lambda t: np.cos(2*np.pi*np.multiply.outer(t, fr)+ph).sum(-1)
    for xq in [0.25, -0.3, 2.4, -7.7, 13.1]:  # between nodes, both sides of the star
        idx, w = bandlimited_weights(nodes, xq, 1.05*band)
        assert len(idx) == 16 and np.all(np.diff(idx) > 0)
        assert abs(w@fn(nodes[idx])-fn(xq)) < 1e-7*np.abs(fn(nodes)).max(), xq


def test_bandlimited_weights_use_all_nodes_when_fewer_than_taps():
    """With fewer nodes than n_taps every node is used; the weights still reproduce a constant
    (to the least-squares fit residual, ~1e-8)."""
    idx, w = bandlimited_weights(np.arange(5.), 1.5, 0.05, n_taps=16)
    assert list(idx) == [0, 1, 2, 3, 4] and w.dtype == np.float64
    assert w.sum() == pytest.approx(1.0, abs=1e-6)


# ================================================================ choose_kernel_size

@pytest.mark.parametrize("n", SIZES)
def test_choose_kernel_size_table_and_choice(n):
    """K is odd and the smallest size whose crop error fits crop_share * error_target; the table lists
    K = 1, 3, ..., 2N-1 with crop errors falling from <= 1 to 0."""
    k, info = P.choose_kernel_size(optics4(n), None, 0.02, NCORES, n_probe=4, verbose=False)
    sizes, crop = info['kernel_sizes'], info['crop']
    assert np.array_equal(sizes, np.arange(1, 2*n, 2))
    assert crop.shape == (n,) and np.all(np.diff(crop) <= 1e-15) and 0 < crop[0] <= 1
    assert abs(crop[-1]) < 1e-12
    assert isinstance(k, int) and k%2 == 1
    assert k == sizes[np.flatnonzero(crop <= 0.8*0.02)[0]]


def test_choose_kernel_size_none_and_zero_image_weigh_pixels_equally():
    """pre_img None, an all-zero image and a uniform image give the same probes and the same answer."""
    o = optics4(32)
    runs = [P.choose_kernel_size(o, img, 0.02, NCORES, n_probe=4, verbose=False)
            for img in [None, np.zeros((32, 32)), np.full((32, 32), 7.0)]]
    for k, info in runs[1:]:
        assert k == runs[0][0] and np.array_equal(info['crop'], runs[0][1]['crop'])


def test_choose_kernel_size_crop_is_the_light_outside_the_window():
    """With one bright pixel (negative values count by magnitude), crop[K] is that pixel's PSF light on the
    image outside the K x K window around it, relative to its light on the image."""
    o = optics4(32)
    img = np.zeros((32, 32))
    img[5, 20] = -3.0
    _, info = P.choose_kernel_size(o, img, 0.02, NCORES, n_probe=4, verbose=False)
    frame = P._probe_worker(5, 20, o)
    for k, crop in zip(info['kernel_sizes'][::4], info['crop'][::4]):
        h = k//2
        inside = frame[max(5-h, 0):5+h+1, max(20-h, 0):20+h+1].sum()
        assert crop == pytest.approx((frame.sum()-inside)/frame.sum(), abs=1e-12)


def test_choose_kernel_size_extreme_targets():
    """A target no size can meet (here negative) gives the largest size 2N-1; a target above 100% gives K = 1."""
    o = optics4(32)
    assert P.choose_kernel_size(o, None, -1.0, NCORES, n_probe=2, verbose=False)[0] == 63
    assert P.choose_kernel_size(o, None, 10.0, NCORES, n_probe=2, verbose=False)[0] == 1


def test_choose_kernel_size_verbose_prints_the_choice(capsys):
    """verbose=True prints the chosen K and its crop error; verbose=False prints nothing."""
    o = optics4(32)
    k, _ = P.choose_kernel_size(o, None, 0.02, NCORES, n_probe=2, verbose=True)
    assert "Kernel size K=%d" % k in capsys.readouterr().out
    P.choose_kernel_size(o, None, 0.02, NCORES, n_probe=2, verbose=False)
    assert capsys.readouterr().out == ""


# ================================================================ spatial_psf_calculation

@pytest.mark.parametrize("make_optics, kwargs, match", [
    (lambda: optics8(32), dict(kernel_size=K, mode='nope'), "mode must be 'direct' or 'rotate'"),
    (lambda: optics8(32), dict(kernel_size=30), "kernel_size must be odd"),
    (lambda: Optics('vortex_as_generic', 32, 2.0, optics8(32).params, module='extra_backends'),
     dict(kernel_size=K, mode='rotate'), "rotationally symmetric"),
], ids=["unknown mode", "even kernel", "rotate needs symmetric optics"])
def test_spatial_psf_calculation_rejects_bad_arguments(make_optics, kwargs, match):
    """Bad arguments are ValueErrors raised before any propagation."""
    with pytest.raises(ValueError, match=match):
        P.spatial_psf_calculation(make_optics(), NCORES, None, validate=False, **kwargs)


def test_spatial_psf_calculation_kernel_size_none_asks_choose_kernel_size(monkeypatch):
    """kernel_size=None hands (optics, pre_img, error_target, num_cores) and, in 'rotate' mode, the
    radial PSFs to choose_kernel_size."""
    class Called(Exception):
        pass

    def fake_choose(*args, **kwargs):
        raise Called(args, kwargs)
    monkeypatch.setattr(P, "choose_kernel_size", fake_choose)
    pre_img, radial_psfs = np.ones((32, 32)), np.zeros((16, 32, 32))
    with pytest.raises(Called) as info:
        P.spatial_psf_calculation(optics8(32), NCORES, radial_psfs, kernel_size=None, error_target=0.07,
                                  pre_img=pre_img, mode="rotate")
    (o, img, target, cores), kwargs = info.value.args
    assert o == optics8(32) and img is pre_img and target == 0.07 and cores == NCORES
    assert kwargs["radial_psfs"] is radial_psfs


@pytest.mark.parametrize("n", MODEL_SIZES)
def test_spatial_model_structure(spatial_model, n):
    """The documented outputs: mode None -> 'rotate' for the vortex, tol None -> tol_for_error_target,
    float32 (n, K, K) kernels, (n,) throughputs, (N, N) tmap, (N, N, 4) blend weights and the settings."""
    d = spatial_model(n)[0]
    n_s = len(d['positions'])
    assert d['mode'] == 'rotate' and d['kernel_size'] == K and d['error_target'] == 0.05
    assert d['tol'] == pytest.approx(P.tol_for_error_target(0.05))
    assert d['psfs'].shape == (n_s, K, K) and d['psfs'].dtype == np.float32
    assert d['throughput'].shape == (n_s,) and d['throughput'].dtype == np.float64
    assert d['tmap'].shape == (n, n) and d['tmap'].dtype == np.float64
    assert d['inds'].shape == d['wts'].shape == (n, n, 4)
    assert d['positions'].shape == (n_s, 2) and d['rings'].shape == (n_s,)
    assert d['r_metric'].shape == d['d_r'].shape == d['d_phi'].shape
    assert Optics.from_settings(d) == optics8(n)
    assert not any(k.startswith('err_') for k in d)  # validate=False


@pytest.mark.parametrize("n", MODEL_SIZES)
def test_spatial_model_kernels_carry_shape_and_tmap_the_throughput(spatial_model, n):
    """kernel x throughput is the exact PSF at the sample (centred on it), and tmap is each pixel's own
    throughput (rotate mode, symmetric optics)."""
    d, o = spatial_model(n)[0], optics8(n)
    r, _ = positions_to_polar(n, d['positions'])
    for target in [3, 9, 15]:
        k = np.argmin(np.abs(r-target))
        T, exact = P._spatial_psf_worker(*d['positions'][k], K, o)
        kernel = d['psfs'][k].astype(float)*d['throughput'][k]
        assert np.abs(kernel-exact).sum()/exact.sum() < 5e-3, target
        assert d['throughput'][k] == pytest.approx(T, abs=1e-3)
    for row, col in [(3, 7), (n//2, n//2), (19, 10)]:
        assert d['tmap'][row, col] == pytest.approx(P._spatial_psf_worker(row, col, K, o)[0], abs=1e-3)


@pytest.mark.parametrize("n", MODEL_SIZES)
def test_spatial_model_kernels_off_axis_peak_at_their_centre(spatial_model, n):
    """Samples beyond 1.5 lambda/D have their kernel's peak on the middle pixel (kernels centred on sources)."""
    d = spatial_model(n)[0]
    r, _ = positions_to_polar(n, d['positions'])
    far = np.flatnonzero(r*optics8(n).pixel_scale > 1.5)
    assert len(far) > 0
    peaks = np.array([peak(d['psfs'][k]) for k in far])
    assert np.all(peaks == K//2)


def test_spatial_model_saved_and_axis_cache_reused(spatial_model, in_tmp, monkeypatch):
    """save_path writes the model (reloadable as it was); axis_products reuses its cache when the
    settings match and recomputes when they differ (kernel_size here)."""
    o = optics8(25)
    d, save_path, axis_path = spatial_model(25)
    f = P.load_spatial_psfs(save_path)
    for key in ['psfs', 'positions', 'rings', 'throughput', 'inds', 'wts', 'tmap', 'ring_radii']:
        assert np.array_equal(f[key], d[key]), key
    assert f['mode'] == 'rotate' and f['kernel_size'] == K and f['img_pixel'] == 25
    shutil.copy(axis_path, "axis.npz")  # work on a copy of the shared cache

    class Recomputed(Exception):
        pass

    def no_runs(*args):
        raise Recomputed
    monkeypatch.setattr(P, "_axis_worker", no_runs)
    cached = P.axis_products(o, NCORES, K, "axis.npz")
    assert cached['kernel_size'] == K and cached['img_pixel'] == 25
    with pytest.raises(Recomputed):
        P.axis_products(o, NCORES, K+2, "axis.npz")


def test_spatial_psf_calculation_refuses_oversized_kernel_cube(spatial_model, in_tmp):
    """A kernel cube above max_kernel_gb is a MemoryError naming the limit (before any kernel is made)."""
    shutil.copy(spatial_model(24)[2], "axis.npz")  # reuse the x-axis runs (a copy of the shared cache)
    with pytest.raises(MemoryError, match="max_kernel_gb"):
        P.spatial_psf_calculation(optics8(24), NCORES, radial(24, 8), kernel_size=K, tol=0.05, validate=False,
                                  axis_cache="axis.npz", max_kernel_gb=1e-9)


def test_generic_optics_take_the_direct_path(monkeypatch):
    """Optics without declared symmetry: mode None -> 'direct', kernels propagated per sample, tmap blended
    from the samples' throughputs. (Rates are stubbed to keep the grid tiny.)"""
    r_metric = np.arange(0.5, 30, 1.0)
    flat = np.full(len(r_metric), 1e-3)
    monkeypatch.setattr(P, "shape_change_rates_direct", lambda *args, **kwargs: (r_metric, flat, flat))
    g = Optics('vortex_as_generic', 24, 1.5, optics8(24).params, module='extra_backends')
    d = P.spatial_psf_calculation(g, NCORES, None, kernel_size=7, tol=0.05, validate=False)
    assert d['mode'] == 'direct' and d['optics_kind'] == 'vortex_as_generic'
    assert np.allclose(d['tmap'], blend_map(d['inds'], d['wts'], d['throughput']))
    for k in [0, len(d['positions'])//2, len(d['positions'])-1]:
        T, exact = P._spatial_psf_worker(*d['positions'][k], 7, g)
        assert d['throughput'][k] == T
        assert np.allclose(d['psfs'][k]*T, exact, rtol=1e-6, atol=1e-7*exact.max())


# ================================================================ planets: exact_point_sources, spatial_psf_planets

@pytest.mark.parametrize("n", SIZES)
def test_exact_point_sources_land_on_their_pixels_and_add_linearly(n):
    """Each source is its exact frame times its value; several sources add up."""
    o = optics4(n)
    one = P.exact_point_sources(o, [(6, 24)], [2.0])
    two = P.exact_point_sources(o, [(6, 24), (25.0, 8.5)], [2.0, 0.5])
    assert one.shape == (n, n) and one.dtype == np.float64
    assert tuple(peak(one)) == (6, 24)
    assert np.allclose(one, 2.0*P._source_at_pixel(o, 6, 24)[0].frame())
    assert np.allclose(two-one, 0.5*P._source_at_pixel(o, 25.0, 8.5)[0].frame())


def test_exact_point_sources_skip_sources_off_the_image(capsys):
    """Sources more than half a pixel outside the image are skipped with a message; -0.5 and N-0.5 are kept."""
    o = optics4(32)
    img = P.exact_point_sources(o, [(-0.5, 10), (12, -0.51), (10, 31.5), (32.0, 3)], [1.0, 1.0, 1.0, 1.0])
    out = capsys.readouterr().out
    assert "Planet #2 is outside" in out and "Planet #4 is outside" in out
    assert "Planet #1" not in out and "Planet #3" not in out
    kept = P.exact_point_sources(o, [(-0.5, 10), (10, 31.5)], [1.0, 1.0])
    assert np.allclose(img, kept)


def test_exact_point_sources_empty_list_gives_zero_image():
    """No sources: an all-zero (N, N) float64 image."""
    img = P.exact_point_sources(optics4(33), [], [])
    assert img.shape == (33, 33) and img.dtype == np.float64 and not img.any()


@pytest.mark.parametrize("n", SIZES)
def test_spatial_psf_planets_use_fits_orientation(n):
    """spatial_psf_planets places a planet at arcsec_to_pixel: +y (angle +90 deg) is above the star (larger row)."""
    o, scale = optics4(n), 0.01
    pos = [0.09, np.pi/2]
    img = P.spatial_psf_planets([pos], [3.0], scale, o.settings())
    expected = arcsec_to_pixel(pos, scale, n)
    assert np.allclose(img, P.exact_point_sources(o, [expected], [3.0]))
    assert peak(img)[0] > o.centre and np.abs(peak(img)-expected).max() <= 0.5


# ================================================================ planets: cir_psf_planets (circular method)

@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("angle_deg", [0, 35, 90, 180, -60])
def test_cir_psf_planets_orientation_matches_exact_psf(in_tmp, n, angle_deg):
    """Regression (orientation): the circular method puts a planet where the exact PSF and
    arcsec_to_pixel put it (angle from +x towards +y = increasing row), with the same flux."""
    o, scale = optics4(n), 0.01
    np.save("psfs.npy", radial(n, 4))
    pos = [0.10, np.radians(angle_deg)]
    img = P.cir_psf_planets([pos], [1.0], scale, n, "psfs.npy")
    exact = P.exact_point_sources(o, [arcsec_to_pixel(pos, scale, n)], [1.0])
    assert img.shape == (n, n)
    assert np.array_equal(peak(img), peak(exact))
    assert np.abs(peak(img)-arcsec_to_pixel(pos, scale, n)).max() <= 1
    assert img.sum() == pytest.approx(exact.sum(), rel=0.03)
    if angle_deg == 90:
        assert peak(img)[0] > o.centre  # +y is above the star (FITS)


def test_cir_psf_planets_skip_planets_beyond_the_last_radial_psf(in_tmp, capsys):
    """A planet beyond the last radial PSF is silently skipped; the others are kept."""
    np.save("psfs.npy", radial(32, 4))
    far = P.cir_psf_planets([[0.2, 0.3]], [1.0], 0.01, 32, "psfs.npy")  # 20 px > 15.5 px
    assert far.shape == (32, 32) and not far.any()
    both = P.cir_psf_planets([[0.2, 0.3], [0.1, 0.3]], [1.0, 1.0], 0.01, 32, "psfs.npy")
    assert np.allclose(both, P.cir_psf_planets([[0.1, 0.3]], [1.0], 0.01, 32, "psfs.npy"))
    assert capsys.readouterr().out == ""


# ================================================================ no main guard (regression)

def test_radial_psf_series_works_in_a_script_without_main_guard(in_tmp):
    """Regression: a script computing the radial PSFs in parallel without `if __name__ == "__main__":`
    must finish (it hung with multiprocessing on macOS)."""
    script = in_tmp/"noguard.py"
    script.write_text("from toycoronagraph.optics import Optics\n"
                      "from toycoronagraph.psf import radial_psf_series\n"
                      "print(radial_psf_series(Optics.vortex(2, 16, 1, 0.95), 2).shape)\n")
    try:
        out = subprocess.run([sys.executable, str(script)], cwd=in_tmp, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        pytest.fail("radial_psf_series hung in a script without a main guard")
    assert out.returncode == 0, out.stderr[-500:]
    assert "(8, 16, 16)" in out.stdout
