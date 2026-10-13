"""
Unit tests for stage F3 (tests/scene_stage3_wind.py): the stellar wind in Thomson-scattered light.

Covers toycoronagraph/sources/thomson.py: wind_constant (the K_w formula and its scalings), _beyond
and pixel_integrals (the exact pixel integral of theta^-3: the total, square-annulus sums, the star
pixel(s), symmetry about the star at (N-1)/2, any pixel unit, supersampling, the numerical fallback
for a star larger than half a pixel), profile_factor (Thomson phase function, beta law), and the Wind
module (parameters and errors, render = K_w x pixel integrals, notes, the options mu_e, r_in,
phase_function and velocity_law), in a Scene (contrast or Jy/arcsec^2, strict grid, also with a FITS
image) and through observe() (known symmetric: the circular method). Small grids: 1-18 px.
"""
import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits
from toycoronagraph import params as P
from toycoronagraph.scene import Scene, Star, Grid
from toycoronagraph.sources.images import FitsImage
from toycoronagraph.sources.thomson import Wind, wind_constant, pixel_integrals, profile_factor, _beyond
from toycoronagraph.observe import observe

SIZES = [15, 16]                       # odd (a pixel on the star) and even (the star between four pixels)
EPS_ERI = Star(distance=3.22, radius=0.73)                  # contrast scenes (no stellar flux)
MDOT, V = 30*2e-14, 400.0                                   # wind_image.ipynb's eps Eri wind
LAMBDA, D = 500e-9, 6.0
ARCSEC = (1*u.arcsec).to_value(u.rad)
# theta_in of eps Eri (rad): R_star / d
THETA_STAR = (0.73*u.R_sun/(3.22*u.pc)).decompose().value


def total_integral(n, pixel, theta_in):
    """The integral of theta^-3 over an n x n grid of `pixel` minus the disc theta < theta_in: 2 pi/theta_in - 4 sqrt2/A."""
    return 2*np.pi/theta_in - 4*np.sqrt(2)/(n*pixel/2)


def distances(n):
    """(n, n) distances of the pixel centres from the star at (n-1)/2, in pixels."""
    c = (n-1)/2
    rows, cols = np.indices((n, n))
    return np.hypot(rows-c, cols-c)


def supersampled(n, theta_in, row, col, sub=500):
    """One pixel's integral of theta^-3 (pixel = 1, theta >= theta_in) by the midpoint rule on sub x sub points."""
    c = (n-1)/2
    t = (np.arange(sub)+0.5)/sub-0.5
    r = np.hypot(((col-c)+t)[None, :], ((row-c)+t)[:, None])
    return np.sum(np.where(r >= theta_in, r**-3.0, 0.0))/sub**2


def wind_grid(n, pixel_scale=0.01):
    """A grid of n px of pixel_scale arcsec at 500 nm, 6 m (the notebook's telescope)."""
    return Grid(n, pixel_scale, LAMBDA, D)


# ---------------------------------------------------------------- K_w

def test_wind_constant_matches_the_formula():
    """K_w = sigma_T Mdot / (32 pi v mu_e m_u d), from SI constants written out."""
    sigma_t, m_sun, year, m_u, parsec = 6.6524587321e-29, 1.988409870698051e30, 365.25*86400, 1.66053906660e-27, 3.0856775814913673e16
    expected = sigma_t*MDOT*m_sun/year/(32*np.pi*V*1e3*1.0*m_u*3.22*parsec)
    assert wind_constant(MDOT, V, 3.22) == pytest.approx(expected, rel=1e-6)
    assert wind_constant(MDOT, V, 3.22) == pytest.approx(3.79e-16, rel=0.01)   # the notebook's eps Eri value


def test_wind_constant_scalings():
    """K_w is proportional to Mdot and to 1 / (v mu_e d)."""
    k = wind_constant(MDOT, V, 3.22)
    assert wind_constant(2*MDOT, V, 3.22) == pytest.approx(2*k, rel=1e-12)
    assert wind_constant(MDOT, 2*V, 3.22) == pytest.approx(k/2, rel=1e-12)
    assert wind_constant(MDOT, V, 6.44) == pytest.approx(k/2, rel=1e-12)
    assert wind_constant(MDOT, V, 3.22, mu_e=1.17) == pytest.approx(k/1.17, rel=1e-12)


# ---------------------------------------------------------------- the exact pixel integral

def test_beyond_is_the_quadrant_integral():
    """_beyond(x, y) = integral of r^-3 over [x, inf) x [y, inf): 1/y on an axis; corner sums give rectangles."""
    from scipy import integrate
    assert _beyond(0.0, 2.0) == pytest.approx(0.5, rel=1e-15) and _beyond(4.0, 0.0) == pytest.approx(0.25, rel=1e-15)
    a, b, c, d = 1.0, 2.0, 0.5, 3.0
    rect = _beyond(a, c) - _beyond(a, d) - _beyond(b, c) + _beyond(b, d)
    numeric = integrate.dblquad(lambda y, x: np.hypot(x, y)**-3, a, b, c, d, epsabs=1e-13, epsrel=1e-12)[0]
    assert rect == pytest.approx(numeric, rel=1e-10)
    assert np.allclose(_beyond(np.array([1.0, 2.0]), np.array([[0.0], [1.0]])),
                       [[2/2, 2/4], [2/(2+np.sqrt(2)), 2/(3+np.sqrt(5))]], rtol=1e-15)   # broadcast


@pytest.mark.parametrize("n", SIZES)
def test_pixel_integrals_shape_dtype_and_flag(n):
    """pixel_integrals returns ((n, n) float64, bool fallback); positive and finite everywhere for a small star."""
    img, fallback = pixel_integrals(n, 1.0, 0.2)
    assert img.shape == (n, n) and img.dtype == np.float64 and fallback is False
    assert np.all(np.isfinite(img)) and np.all(img > 0)


@pytest.mark.parametrize("n", [15, 16, 17, 18])
@pytest.mark.parametrize("theta_in", [0.01, 0.2, 0.45])
def test_total_is_the_analytic_integral(n, theta_in):
    """The image sums to 2 pi / theta_in - 4 sqrt2 / A (A: the half width) to round-off."""
    img, _ = pixel_integrals(n, 1.0, theta_in)
    assert img.sum() == pytest.approx(total_integral(n, 1.0, theta_in), rel=1e-13)


@pytest.mark.parametrize("n", [1, 2, 3])
def test_smallest_grids_keep_the_analytic_total(n):
    """Edge case: one pixel on the star (n = 1) or four around it (n = 2) still sum to the analytic integral."""
    img, fallback = pixel_integrals(n, 1.0, 0.2)
    assert img.shape == (n, n) and fallback is False
    assert img.sum() == pytest.approx(total_integral(n, 1.0, 0.2), rel=1e-13)


@pytest.mark.parametrize("n", [15, 16, 17, 18])
def test_square_annulus_sums(n):
    """Pixels between two centred squares of half widths A1 < A2 sum to 4 sqrt2 (1/A1 - 1/A2)."""
    img, _ = pixel_integrals(n, 1.0, 0.2)
    total = img.sum()
    for k in range(1, n//2):  # the central block of side n - 2k (half width A1 = n/2 - k)
        inner = img[k:n-k, k:n-k].sum()
        a1, a2 = n/2-k, n/2
        assert total-inner == pytest.approx(4*np.sqrt(2)*(1/a1-1/a2), rel=1e-11), k


def test_star_pixel_of_an_odd_grid():
    """Odd n: the star pixel holds 2 pi / theta_in - 8 sqrt2 (its square minus the disc), and is the brightest."""
    img, _ = pixel_integrals(15, 1.0, 0.01)
    assert img[7, 7] == pytest.approx(2*np.pi/0.01 - 8*np.sqrt(2), rel=1e-13)
    assert np.unravel_index(np.argmax(img), img.shape) == (7, 7)


def test_star_pixels_of_an_even_grid():
    """Even n: each of the four pixels touching the star holds pi / (2 theta_in) - sqrt2."""
    img, _ = pixel_integrals(16, 1.0, 0.01)
    assert np.allclose(img[7:9, 7:9], np.pi/(2*0.01) - np.sqrt(2), rtol=1e-13, atol=0)
    assert img[7:9, 7:9].min() > np.delete(img.ravel(), [7*16+7, 7*16+8, 8*16+7, 8*16+8]).max()


@pytest.mark.parametrize("n", [15, 16, 17, 18])
def test_symmetric_about_the_star(n):
    """The image is symmetric about the star at (n-1)/2: transposed and flipped both ways it is the same."""
    img, _ = pixel_integrals(n, 1.0, 0.2)
    for other in (img.T, img[::-1], img[:, ::-1], np.rot90(img)):
        assert np.allclose(other, img, rtol=1e-13, atol=0)


@pytest.mark.parametrize("n", SIZES)
def test_any_pixel_unit(n):
    """In any angle unit: pixel_integrals(n, p, theta p) = pixel_integrals(n, 1, theta) / p (in the inverse unit)."""
    p = 0.01*ARCSEC
    a, fa = pixel_integrals(n, p, 0.2*p)
    b, fb = pixel_integrals(n, 1.0, 0.2)
    assert fa == fb and np.allclose(a*p, b, rtol=1e-12, atol=0)


@pytest.mark.parametrize("n", [11, 12])
def test_pixels_match_supersampling(n):
    """Pixels at least 2 px from the star equal a fine midpoint-rule integral (to 1e-5)."""
    img, _ = pixel_integrals(n, 1.0, 0.2)
    c = int(np.floor((n-1)/2))
    for row, col in [(c+2, c), (c+2, c+3), (c-3, c+1), (c+5, c-4), (c, c-2)]:
        assert img[row, col] == pytest.approx(supersampled(n, 0.2, row, col), rel=1e-5), (row, col)


@pytest.mark.parametrize("n", SIZES)
def test_pixel_mean_exceeds_the_centre_value_by_three_eighths_over_r2(n):
    """Away from the star a pixel's mean of r^-3 is its centre value x (1 + 3/(8 r^2)) (second-order Taylor)."""
    img, _ = pixel_integrals(n, 1.0, 0.2)
    r = distances(n)
    far = r >= 3
    excess = (img[far]*r[far]**3 - 1)*r[far]**2
    assert np.all(excess > 0.375) and np.all(excess < 0.375*1.1)


def test_fallback_flag_switches_at_half_a_pixel():
    """The numerical fallback is used when theta_in >= half a pixel, and only then."""
    assert pixel_integrals(15, 1.0, 0.4999)[1] is False
    assert pixel_integrals(15, 1.0, 0.5)[1] is True
    assert pixel_integrals(15, 2.0, 0.9999)[1] is False


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("theta_in", [0.7, 1.7])
def test_fallback_keeps_the_total_and_the_far_pixels(n, theta_in):
    """With the fallback the total is within 1e-4 of the analytic one, pixels inside theta_in are 0, far pixels exact."""
    img, fallback = pixel_integrals(n, 1.0, theta_in)
    exact_far, _ = pixel_integrals(n, 1.0, 0.2)   # pixels not touching the disc do not depend on theta_in
    r = distances(n)
    assert fallback is True
    assert img.sum() == pytest.approx(total_integral(n, 1.0, theta_in), rel=1e-4)
    assert np.all(img[r + np.sqrt(0.5) < theta_in] == 0)
    far = r >= theta_in + 1
    assert np.array_equal(img[far], exact_far[far])


# ---------------------------------------------------------------- profile_factor

def test_profile_factor_is_one_for_isotropic_constant_speed():
    """Constant speed and isotropic scattering: the factor is 1 at any impact parameter."""
    for b in (1.0, 3.0, 100.0):
        assert profile_factor(b, 0.9, V, 10.0, None, "isotropic") == pytest.approx(1.0, rel=1e-12)


def test_profile_factor_thomson_is_fifteen_sixteenths():
    """The Thomson phase function (3/16 pi)(1 + cos^2) gives exactly 15/16 of the isotropic light at constant speed."""
    assert profile_factor(2.0, 0.9, V, 10.0, None, "thomson") == pytest.approx(15/16, rel=1e-12)


def test_profile_factor_beta_law():
    """Beta law: denser (brighter) near r_in, tending to 1 far away; v0 = v is constant speed."""
    near, mid, far = (profile_factor(b, 1.0, V, 10.0, 1.0, "isotropic") for b in (1.01, 10.0, 1000.0))
    assert near > mid > far > 1 and near > 5 and far == pytest.approx(1.0, abs=1e-3)
    assert profile_factor(1.5, 1.0, V, V, 1.0, "isotropic") == pytest.approx(1.0, rel=1e-12)


# ---------------------------------------------------------------- the Wind module

def test_wind_declares_its_family_requirements_and_symmetry():
    """Wind: kind "wind", a contrast ("relative"), needs the star's distance and radius, symmetric by construction."""
    w = Wind(mdot=MDOT)
    assert w.kind == "wind" and w.family_of() == "relative"
    assert w.requires == ("star.distance", "star.radius") and w.symmetric() is True


def test_wind_defaults():
    """The defaults are the notebooks' wind: 400 km/s, mu_e 1, from the stellar surface, isotropic, constant speed."""
    s = Wind(mdot=MDOT).settings
    assert s == dict(mdot=MDOT, v=400.0, mu_e=1.0, r_in=1.0, velocity_law="constant", beta=1.0, v0=10.0,
                     phase_function="isotropic")


def test_mdot_is_required():
    """mdot has no default: Wind().missing() names it, and rendering a scene with it stops with a ParamError."""
    assert [p.name for p in Wind().missing()] == ["mdot"]
    scene = Scene(wavelength=LAMBDA, diameter=D, sampling=4, size=15, star=EPS_ERI)
    scene.add(Wind())
    with pytest.raises(P.ParamError, match="mdot"):
        scene.render(verbose=False)


@pytest.mark.parametrize("settings, words", [
    (dict(mdot=0), "must be > 0"),
    (dict(mdot=-1e-13), "must be > 0"),
    (dict(mdot=1e-13, v=0), "must be > 0"),
    (dict(mdot=1e-13, r_in=0.5), "must be >= 1"),
    (dict(mdot=1e-13, vv=3), "did you mean v"),
    (dict(mdot=1e-13, velocity_law="betta"), "not allowed"),
    (dict(mdot=1e-13, phase_function="rayleigh"), "not allowed"),
])
def test_invalid_settings_are_param_errors(settings, words):
    """Out-of-range values, unknown names (with a suggestion) and unknown choices: ParamError at creation."""
    with pytest.raises(P.ParamError, match=words):
        Wind(**settings)


def test_mdot_with_a_unit_in_text():
    """mdot may be text with a unit ("Msun" accepted), converted to M_sun/yr."""
    assert Wind(mdot="6e-13 Msun/yr").settings["mdot"] == pytest.approx(6e-13, rel=1e-12)
    g_per_s = float((1e-13*u.M_sun/u.yr).to_value(u.g/u.s))
    assert Wind(mdot="%.12g g/s" % g_per_s).settings["mdot"] == pytest.approx(1e-13, rel=1e-10)


@pytest.mark.parametrize("n", SIZES)
def test_render_is_k_w_times_the_pixel_integrals(n):
    """render(): contrast_pixel on the scene grid (pixel_scale None) = K_w x pixel_integrals in rad."""
    g = wind_grid(n)
    img = Wind(mdot=MDOT, v=V).render(g, EPS_ERI)
    expected = wind_constant(MDOT, V, 3.22)*pixel_integrals(n, g.pixel_scale*ARCSEC, THETA_STAR)[0]
    assert img.unit == "contrast_pixel" and img.pixel_scale is None and "wind" in img.origin
    assert img.data.shape == (n, n) and img.data.dtype == np.float64
    assert np.allclose(img.data, expected, rtol=1e-13, atol=0)


@pytest.mark.parametrize("n", SIZES)
def test_far_from_the_star_the_contrast_is_k_w_over_theta_cubed(n):
    """Per steradian, a pixel away from the star (>= 3 px) holds K_w / theta^3 at its centre, to (1 + 0.42/r^2)."""
    g = wind_grid(n)
    per_sr = Wind(mdot=MDOT, v=V).render(g, EPS_ERI).data/(g.pixel_scale*ARCSEC)**2
    r = distances(n)
    far = r >= 3
    ratio = per_sr[far]/(wind_constant(MDOT, V, 3.22)/(r[far]*g.pixel_scale*ARCSEC)**3)
    assert np.all(ratio > 1) and np.all(ratio < 1 + 0.42/r[far]**2)


@pytest.mark.parametrize("n", SIZES)
def test_even_grid_gets_a_note(n):
    """An even grid size gets a note (no pixel on the star); an odd one does not."""
    notes = Wind(mdot=MDOT).render(wind_grid(n), EPS_ERI).notes
    assert any("even grid size" in note for note in notes) == (n % 2 == 0)


def test_star_larger_than_half_a_pixel_gets_the_fallback_note():
    """A star wider than half a pixel (1 mas pixels for eps Eri): the numerical fallback, with a note."""
    img = Wind(mdot=MDOT).render(wind_grid(15, pixel_scale=0.001), EPS_ERI)
    assert any("integrated numerically" in note for note in img.notes)
    assert not any("integrated numerically" in note for note in Wind(mdot=MDOT).render(wind_grid(15), EPS_ERI).notes)


def test_r_in_moves_the_inner_edge():
    """r_in = 2: the wind starts at twice the stellar radius (the same pixel integrals with 2 theta_star)."""
    g = wind_grid(15, pixel_scale=0.005)
    img = Wind(mdot=MDOT, r_in=2.0).render(g, EPS_ERI).data
    expected = wind_constant(MDOT, 400.0, 3.22)*pixel_integrals(15, g.pixel_scale*ARCSEC, 2*THETA_STAR)[0]
    assert np.allclose(img, expected, rtol=1e-13, atol=0)


@pytest.mark.parametrize("n", SIZES)
def test_mu_e_and_the_thomson_phase_function_scale_the_image(n):
    """mu_e divides the image by mu_e; the Thomson phase function multiplies it by 15/16 (constant speed)."""
    g = wind_grid(n)
    base = Wind(mdot=MDOT).render(g, EPS_ERI).data
    assert np.allclose(Wind(mdot=MDOT, mu_e=1.17).render(g, EPS_ERI).data, base/1.17, rtol=1e-12, atol=0)
    assert np.allclose(Wind(mdot=MDOT, phase_function="thomson").render(g, EPS_ERI).data, base*15/16, rtol=1e-12, atol=0)


def test_beta_law_is_brighter_near_the_star():
    """A beta-law wind is denser near the star: brighter than constant speed, most at the star."""
    g = wind_grid(15, pixel_scale=0.0025)   # 2.4 stellar radii per pixel
    ratio = (Wind(mdot=MDOT, velocity_law="beta").render(g, EPS_ERI).data /
             Wind(mdot=MDOT).render(g, EPS_ERI).data)
    assert np.all(ratio > 1) and ratio[7, 7] == ratio.max() and ratio[7, 7] > 2
    assert ratio[7, 8] > ratio[7, 10] > ratio[7, 14]


def test_beta_law_is_the_same_far_away():
    """Hundreds of stellar radii from the star a beta-law wind is the constant-speed one to 1%."""
    g = wind_grid(15, pixel_scale=0.5)
    ratio = (Wind(mdot=MDOT, velocity_law="beta").render(g, EPS_ERI).data /
             Wind(mdot=MDOT).render(g, EPS_ERI).data)
    off_star = distances(15) >= 1
    assert np.all(np.abs(ratio[off_star]-1) < 0.01)


# ---------------------------------------------------------------- in a scene

@pytest.mark.parametrize("n", SIZES)
def test_scene_without_the_star_flux_is_a_contrast(n):
    """A wind and a star without temperature: the scene is in contrast_pixel and holds render()'s image unchanged."""
    scene = Scene(wavelength=LAMBDA, diameter=D, sampling=4, size=n, star=EPS_ERI)
    scene.add(Wind(mdot=MDOT, v=V))
    scene.render(verbose=False)
    assert scene.internal_unit == "contrast_pixel" and scene.families == {"wind": "relative"} and scene.star_flux is None
    assert np.allclose(scene.images["wind"], Wind(mdot=MDOT, v=V).render(scene.grid, scene.star).data, rtol=1e-14, atol=0)


def test_scene_with_the_star_flux_is_a_brightness():
    """With the star's temperature the scene is in Jy/arcsec^2: contrast x the star's flux / the pixel area."""
    contrast = Scene(wavelength=LAMBDA, diameter=D, sampling=4, size=15, star=EPS_ERI)
    contrast.add(Wind(mdot=MDOT))
    contrast.render(verbose=False)
    jy = Scene(wavelength=LAMBDA, diameter=D, sampling=4, size=15, star=Star(distance=3.22, radius=0.73, temperature=5080))
    jy.add(Wind(mdot=MDOT))
    jy.render(verbose=False)
    assert jy.internal_unit == "Jy/arcsec^2"
    assert np.allclose(jy.images["wind"], contrast.images["wind"]*jy.star_flux/jy.grid.pixel_area, rtol=1e-12, atol=0)


def test_scene_needs_the_star_radius():
    """Without [star] radius the wind cannot start: a ParamError naming radius."""
    scene = Scene(wavelength=LAMBDA, diameter=D, sampling=4, size=15, star=Star(distance=3.22))
    scene.add(Wind(mdot=MDOT))
    with pytest.raises(P.ParamError, match="radius"):
        scene.render(verbose=False)


def test_strict_grid_accepts_a_wind():
    """[grid] resample = false: a wind is drawn on the grid, so it never stops the run."""
    scene = Scene(wavelength=LAMBDA, diameter=D, sampling=4, size=15, resample=False, star=EPS_ERI)
    scene.add(Wind(mdot=MDOT))
    assert set(scene.render(verbose=False)) == {"wind"}
    assert not any("resampled" in note for note in scene.notes)


def write_flat_fits(path, n, pixel_scale):
    """Write an (n, n) flat contrast image of 1e-6 per pixel with PIXSCALE pixel_scale arcsec to path."""
    header = fits.Header([("BUNIT", "contrast"), ("PIXSCALE", pixel_scale)])
    fits.PrimaryHDU(np.full((n, n), 1e-6), header).writeto(str(path), overwrite=True)


@pytest.mark.parametrize("scale, fails", [(1.0, False), (1.3, True)])
def test_strict_grid_names_only_the_mismatching_image(scale, fails, in_tmp):
    """resample = false with a wind and a FITS disk: a matching disk runs; a mismatching one stops the run with a
    ParamError naming the file, not the wind (the wind is drawn on the grid)."""
    pixel_scale = LAMBDA/D*P.RAD_TO_ARCSEC/4
    write_flat_fits(in_tmp/"disk.fits", 15, scale*pixel_scale)
    scene = Scene(wavelength=LAMBDA, diameter=D, sampling=4, size=15, resample=False, star=EPS_ERI)
    scene.add(Wind(mdot=MDOT))
    scene.add(FitsImage(file=str(in_tmp/"disk.fits")))
    if not fails:
        assert set(scene.render(verbose=False)) == {"wind", "fits_image"}
        assert scene.grid.pixel_scale == pytest.approx(pixel_scale, rel=1e-12)
        return
    with pytest.raises(P.ParamError) as err:
        scene.render(verbose=False)
    details = str(err.value).split("\n", 1)[1]
    assert "disk.fits" in str(err.value) and "wind" not in details


def test_observe_puts_a_wind_on_the_circular_method(in_tmp):
    """Known symmetric by its module: auto picks the circular method, and the coronagraph removes most of the light."""
    scene = Scene(wavelength=0.8e-6, diameter=2.4, sampling=4, size=15, star=EPS_ERI,
                  compute={"psf_cache_dir": str(in_tmp/"cache"), "num_cores": 2})
    scene.add(Wind(mdot=MDOT, v=V))
    res = observe(scene, ["vortex:6"], verbose=False)
    method, mode, why, error = res.methods["vortex_c6", "wind"]
    assert method == "circular" and mode is None and error == 0.0 and "known by the wind module" in why
    assert 0 < res.final().sum() < 0.5*scene.total.sum()
