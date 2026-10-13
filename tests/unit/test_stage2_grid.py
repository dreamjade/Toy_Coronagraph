"""
Unit tests for Stage 2: the equal-error PSF sample grid and its blend weights
(toycoronagraph/spatial.py), and the sky -> pixel convention of toycoronagraph/tool.py.

Covers spatial.positions_to_polar, equal_error_grid, make_blend_weights, blend_map, make_psf_inds,
and tool.arcsec_to_pixel (FITS orientation: +x -> +column, +y -> +row, star at (N-1)/2).

Synthetic shape-change rates shaped like the measured vortex ones (no hcipy, no files);
even and odd image sizes (32 and 33); about a second in total.
"""
import warnings
import numpy as np
import pytest
from toycoronagraph.spatial import (positions_to_polar, equal_error_grid, make_blend_weights, blend_map,
                                    make_psf_inds)
from toycoronagraph.tool import arcsec_to_pixel, convert_to_polar

SIZES = [32, 33]  # even (star between four pixels) and odd (star on a pixel)
TOL, MAX_SPACING, N_THETA_MIN = 0.05, 4.0, 8

# Synthetic rates (per pixel of radius / per radian), steep near the star like the vortex ones
R_METRIC = np.arange(0.25, 60, 0.5)
D_R = 0.07/(1+(R_METRIC/4)**3)+1e-4
D_PHI = 0.5*2*(R_METRIC/3)/(1+(R_METRIC/3)**2)+0.01


def grid(n, tol=TOL, **kwargs):
    """equal_error_grid of an n-pixel image with the synthetic rates (warnings become errors)."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        return equal_error_grid(n, R_METRIC, D_R, D_PHI, tol, MAX_SPACING, N_THETA_MIN, **kwargs)


def pixel_polar(n):
    """(r, phi) of every pixel of an n x n image about the star, as two (n, n) arrays."""
    rows, cols = np.indices((n, n))
    r, phi = positions_to_polar(n, np.stack([rows.ravel(), cols.ravel()], axis=1))
    return r.reshape(n, n), phi.reshape(n, n)


# ---------------------------------------------------------------- positions_to_polar

@pytest.mark.parametrize("n", SIZES)
def test_positions_to_polar_star_at_centre_and_angle_from_columns_towards_rows(n):
    """The star is at (N-1)/2; phi = 0 along +column, +pi/2 along +row; float64 outputs."""
    c = (n-1)/2
    pos = np.array([[c, c], [c, c+3], [c+3, c], [c, c-3], [c-3, c]])
    r, phi = positions_to_polar(n, pos)
    assert r.dtype == np.float64 and phi.dtype == np.float64
    assert np.allclose(r, [0, 3, 3, 3, 3])
    assert np.allclose(phi[1:], [0, np.pi/2, np.pi, -np.pi/2])


def test_positions_to_polar_accepts_integer_lists():
    """Integer (row, col) positions given as a list work; angles lie in [-pi, pi]."""
    rows, cols = np.indices((9, 9))
    pos = np.stack([rows.ravel(), cols.ravel()], axis=1).tolist()
    r, phi = positions_to_polar(9, pos)
    assert r.shape == (81,) and np.isclose(r.max(), 4*np.sqrt(2))
    assert np.all(np.abs(phi) <= np.pi)


# ---------------------------------------------------------------- equal_error_grid: structure

@pytest.mark.parametrize("n", SIZES)
def test_grid_output_types_and_shapes(n):
    """The documented keys, shapes and dtypes."""
    g = grid(n)
    n_s, n_rings = len(g['positions']), len(g['ring_radii'])
    assert g['positions'].shape == (n_s, 2) and g['positions'].dtype == np.float64
    assert g['rings'].shape == (n_s,) and np.issubdtype(g['rings'].dtype, np.integer)
    assert g['ring_theta0'].shape == g['ring_ntheta'].shape == (n_rings,)
    assert np.issubdtype(g['ring_ntheta'].dtype, np.integer)
    assert isinstance(g['r_pix'], float)
    assert g['inds'].shape == (n, n, 4) and np.issubdtype(g['inds'].dtype, np.integer)
    assert g['wts'].shape == (n, n, 4) and g['wts'].dtype == np.float64


@pytest.mark.parametrize("n", SIZES)
def test_grid_inner_samples_are_the_pixels_inside_r_pix(n):
    """Ring-0 samples are exactly the pixel centres with r < r_pix, in row-major order."""
    g = grid(n)
    pix_r, _ = pixel_polar(n)
    inner = g['positions'][g['rings'] == 0]
    assert np.array_equal(inner, np.argwhere(pix_r < g['r_pix']).astype(float))
    assert np.all(g['rings'][:len(inner)] == 0)  # listed first


@pytest.mark.parametrize("n, n_central", [(32, 4), (33, 1)])
def test_grid_always_includes_central_pixels(n, n_central):
    """Even with no shape change at all, r_pix >= 1: the 4 (even N) or 1 (odd N) central pixels are samples."""
    zero = np.zeros(len(R_METRIC))
    g = equal_error_grid(n, R_METRIC, zero, zero, TOL, MAX_SPACING, N_THETA_MIN)
    assert g['r_pix'] == 1.0
    assert np.sum(g['rings'] == 0) == n_central
    assert np.allclose(np.diff(g['ring_radii'])[:-1], MAX_SPACING)  # spacing capped everywhere


@pytest.mark.parametrize("n", SIZES)
def test_grid_rings_start_at_r_pix_and_reach_past_every_pixel(n):
    """Ring radii start at r_pix, increase by at most max_spacing, and the last reaches the corner."""
    g = grid(n)
    radii = g['ring_radii']
    pix_r, _ = pixel_polar(n)
    assert len(radii) >= 2
    assert np.isclose(radii[0], g['r_pix'])
    assert np.all(np.diff(radii) > 0) and np.all(np.diff(radii) <= MAX_SPACING+1e-9)
    assert radii[-1] >= pix_r.max()


@pytest.mark.parametrize("n", SIZES)
def test_grid_ring_samples_sit_exactly_on_their_ring(n):
    """Ring samples are at their exact polar positions (radius of their ring), not snapped to pixels."""
    g = grid(n)
    on_ring = g['rings'] > 0
    r, _ = positions_to_polar(n, g['positions'][on_ring])
    assert np.allclose(r, g['ring_radii'][g['rings'][on_ring]-1])
    assert not np.allclose(g['positions'][on_ring], np.round(g['positions'][on_ring]))


@pytest.mark.parametrize("n", SIZES)
def test_grid_samples_per_ring_and_stagger(n):
    """Each ring has >= n_theta_min samples, arcs <= max_spacing, an angular shape change <= tol,
    and alternate rings are rotated by half a step."""
    g = grid(n)
    radii, ntheta = g['ring_radii'], g['ring_ntheta']
    assert np.all(ntheta >= N_THETA_MIN)
    assert np.all(2*np.pi*radii/ntheta <= MAX_SPACING+1e-9)
    assert np.all(2*np.pi/ntheta*np.interp(radii, R_METRIC, D_PHI) <= TOL*(1+1e-9))
    assert np.allclose(g['ring_theta0'], np.pi*(np.arange(len(radii))%2)/ntheta)


@pytest.mark.parametrize("n", SIZES)
def test_grid_ring_spacing_follows_equal_shape_change(n):
    """Between rings the accumulated radial shape change is tol, unless the spacing is capped (or the last ring)."""
    g = grid(n)
    radii = g['ring_radii']
    s = np.concatenate([[0], np.cumsum(0.5*(D_R[1:]+D_R[:-1])*np.diff(R_METRIC))])
    ds = np.diff(np.interp(radii, R_METRIC, s))
    capped = np.isclose(np.diff(radii), MAX_SPACING)
    assert np.all(ds <= TOL*(1+1e-6))
    assert np.allclose(ds[~capped][:-1], TOL, rtol=1e-6)


@pytest.mark.parametrize("n", SIZES)
def test_grid_r_inner_max_caps_r_pix(n):
    """r_inner_max bounds the per-pixel zone."""
    free, capped = grid(n), grid(n, r_inner_max=3.0)
    assert free['r_pix'] > 3.0 and capped['r_pix'] == 3.0
    assert np.sum(capped['rings'] == 0) < np.sum(free['rings'] == 0)


@pytest.mark.parametrize("n", SIZES)
def test_grid_finer_tolerance_gives_more_samples(n):
    """A smaller tol gives more samples and a larger (or equal) per-pixel zone."""
    coarse, fine = grid(n, tol=0.08), grid(n, tol=0.03)
    assert len(fine['positions']) > len(coarse['positions'])
    assert fine['r_pix'] >= coarse['r_pix']


# ---------------------------------------------------------------- equal_error_grid: blend weights

@pytest.mark.parametrize("n", SIZES)
def test_grid_blend_weights_are_valid(n):
    """Weights are finite, >= 0 and sum to 1 per pixel; indices point at kept samples, 0 where unused,
    and every kept sample is used by some pixel."""
    g = grid(n)
    inds, wts, n_s = g['inds'], g['wts'], len(g['positions'])
    assert np.all(np.isfinite(wts)) and np.all(wts >= 0)
    assert np.allclose(wts.sum(-1), 1)
    assert inds.min() >= 0 and inds.max() < n_s
    assert np.all(inds[wts == 0] == 0)
    assert len(np.unique(inds[wts > 0])) == n_s


@pytest.mark.parametrize("n", SIZES)
def test_grid_inner_pixels_use_their_own_sample(n):
    """A pixel inside r_pix uses exactly its own sample, with weight 1 in slot 0."""
    g = grid(n)
    pix_r, _ = pixel_polar(n)
    rows, cols = np.nonzero(pix_r < g['r_pix'])
    assert np.all(g['wts'][rows, cols, 0] == 1) and np.all(g['wts'][rows, cols, 1:] == 0)
    own = g['positions'][g['inds'][rows, cols, 0]]
    assert np.array_equal(own, np.stack([rows, cols], axis=1).astype(float))


@pytest.mark.parametrize("n", SIZES)
def test_blend_map_reproduces_a_function_linear_in_radius(n):
    """Blending the samples' radii gives every pixel's own radius."""
    g = grid(n)
    r, _ = positions_to_polar(n, g['positions'])
    pix_r, _ = pixel_polar(n)
    assert np.abs(blend_map(g['inds'], g['wts'], r)-pix_r).max() < 1e-9


@pytest.mark.parametrize("n", SIZES)
def test_outer_pixels_blend_the_two_bracketing_rings(n):
    """Outside r_pix, the samples a pixel uses lie on the rings at and beyond its radius."""
    g = grid(n)
    r, _ = positions_to_polar(n, g['positions'])
    pix_r, _ = pixel_polar(n)
    outer = pix_r >= g['r_pix']
    used_r = np.where(g['wts'] > 0, r[g['inds']], np.nan)
    assert np.all(np.nanmin(used_r[outer], axis=-1) <= pix_r[outer]+1e-9)
    assert np.all(np.nanmax(used_r[outer], axis=-1) >= pix_r[outer]-1e-9)


@pytest.mark.parametrize("n", SIZES)
def test_grid_weights_finite_when_per_pixel_zone_reaches_the_corner(n):
    """Regression guard: steep rates everywhere (r_pix = corner radius) must still give finite weights summing to 1."""
    steep = np.full(len(R_METRIC), 10.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the bug also raises 'invalid value encountered in divide'
        g = equal_error_grid(n, R_METRIC, steep, steep, TOL, MAX_SPACING, N_THETA_MIN)
    assert np.all(np.isfinite(g['wts'])) and np.allclose(g['wts'].sum(-1), 1)


# ---------------------------------------------------------------- make_blend_weights on a hand-made grid

def hand_grid(n=17):
    """An unpruned grid as equal_error_grid builds it: per-pixel samples inside r_pix = 2.5, then
    four rings (radii 2.5, 5, 8, 12 px; 8, 12, 16, 24 samples; first samples at angle 0)."""
    c = (n-1)/2
    pix_r, _ = pixel_polar(n)
    n_inner = int(np.sum(pix_r < 2.5))
    radii, ntheta = np.array([2.5, 5., 8., 12.]), np.array([8, 12, 16, 24])
    rings = np.concatenate([np.zeros(n_inner, int)]+[np.full(k, i+1) for i, k in enumerate(ntheta)])
    return dict(rings=rings, ring_radii=radii, ring_theta0=np.zeros(4), ring_ntheta=ntheta, r_pix=2.5), c, n_inner


def test_make_blend_weights_pixel_on_a_sample_takes_only_that_sample():
    """A pixel exactly at a ring sample (radius 5 px, angle 0) gets weight 1 on that sample, in slot 0."""
    g, c, n_inner = hand_grid()
    inds, wts = make_blend_weights(17, g)
    row, col = int(c), int(c)+5
    assert wts[row, col, 0] == pytest.approx(1.0) and np.allclose(wts[row, col, 1:], 0)
    assert inds[row, col, 0] == n_inner+8  # first sample of the second ring (after 8 samples of the first)


def test_make_blend_weights_slots_inner_ring_then_outer_ring():
    """Slots 0-1 hold the two samples of the ring at or inside the pixel, slots 2-3 those of the next ring."""
    g, c, n_inner = hand_grid()
    inds, wts = make_blend_weights(17, g)
    starts = n_inner+np.concatenate([[0], np.cumsum(g['ring_ntheta'])])  # first index of each ring, and the end
    ring_of = np.searchsorted(starts, inds, side='right')                 # 1-based ring of each slot
    pix_r, _ = pixel_polar(17)
    between = (pix_r >= 5) & (pix_r < 8)                                  # between rings 2 and 3
    assert np.all(ring_of[between][:, :2] == 2) and np.all(ring_of[between][:, 2:] == 3)
    assert np.allclose(wts.sum(-1), 1) and np.all(wts >= 0)


def test_blend_map_is_weighted_sum_of_sample_values():
    """blend_map(inds, wts, values) = sum_k wts[..., k] * values[inds[..., k]] (values may be a list)."""
    rng = np.random.default_rng(0)
    inds = rng.integers(0, 10, (5, 6, 4))
    wts = rng.random((5, 6, 4))
    values = rng.random(10)
    expected = sum(wts[..., k]*values[inds[..., k]] for k in range(4))
    out = blend_map(inds, wts, list(values))
    assert out.shape == (5, 6) and np.allclose(out, expected)


@pytest.mark.parametrize("n", SIZES)
def test_make_psf_inds_assigns_the_nearest_sample(n):
    """Every pixel gets the index of its nearest sample position."""
    p = grid(n)['positions']
    inds = make_psf_inds(n, p)
    assert inds.shape == (n, n) and np.issubdtype(inds.dtype, np.integer)
    rows, cols = np.indices((n, n))
    d_own = np.hypot(rows-p[inds, 0], cols-p[inds, 1])
    d_all = np.hypot(rows[..., None]-p[:, 0], cols[..., None]-p[:, 1]).min(axis=-1)
    assert np.allclose(d_own, d_all)


# ---------------------------------------------------------------- tool.arcsec_to_pixel (FITS orientation)

@pytest.mark.parametrize("n", SIZES)
def test_arcsec_to_pixel_fits_orientation(n):
    """Star at (N-1)/2; +x -> +column, +y -> +row (FITS), -x -> -column, -y -> -row."""
    s, c = 0.1, (n-1)/2
    assert np.allclose(arcsec_to_pixel([0.0, 1.0], s, n), [c, c])
    assert np.allclose(arcsec_to_pixel([1.0, 0.0], s, n), [c, c+10])
    assert np.allclose(arcsec_to_pixel([1.0, np.pi/2], s, n), [c+10, c])
    assert np.allclose(arcsec_to_pixel([0.5, np.pi], s, n), [c, c-5])
    assert np.allclose(arcsec_to_pixel([0.5, -np.pi/2], s, n), [c-5, c])


def test_arcsec_to_pixel_returns_fractional_row_col():
    """The (2,) float result need not be an integer pixel; the radius is divided by the pixel scale."""
    out = arcsec_to_pixel([0.25, np.pi/6], 0.02, 32)
    assert out.shape == (2,) and out.dtype == np.float64
    assert np.allclose(out, [15.5+12.5*0.5, 15.5+12.5*np.sqrt(3)/2])


@pytest.mark.parametrize("n", SIZES)
def test_arcsec_to_pixel_round_trips_with_positions_to_polar(n):
    """positions_to_polar inverts arcsec_to_pixel: same radius (in pixels) and angle, any quadrant."""
    s = 0.05
    for r, theta in [(0.3, 0.2), (0.6, 2.0), (0.45, -2.5), (0.7, -0.9)]:
        pr, pphi = positions_to_polar(n, arcsec_to_pixel([r, theta], s, n)[None, :])
        assert np.isclose(pr[0], r/s) and np.isclose(pphi[0], theta)


@pytest.mark.parametrize("n", SIZES)
def test_cartesian_sky_position_lands_at_row_y_col_x(n):
    """A planet at cartesian (x, y) arcsec, through convert_to_polar, lands at (c + y/s, c + x/s)."""
    s, c = 0.05, (n-1)/2
    for x, y in [(0.3, 0.2), (-0.4, 0.1), (0.25, -0.5)]:
        assert np.allclose(arcsec_to_pixel(convert_to_polar([x, y]), s, n), [c+y/s, c+x/s])
