"""
Unit tests for stage F4 (tests/scene_stage4_cme.py): the CME module.

Covers toycoronagraph/sources/cme.py: cme_mass_fraction (the geometry ported from astro_toolbox's
cme_imaging: shape, odd size, normalisation, mirror symmetry about the star, direction +x, the
settings, small grids), cme_contrast (N_e sigma_T / (4 pi r^2), the star masked), model_arguments,
load_cme_imaging errors, contrast_of / star_photons / Cme.from_object with a stand-in for a
cme_imaging.Star (its grid and star provided, the factor-pi and NaN notes, the distance conflict
rule), and the Cme module (defaults, parameters and errors, render on odd and even grids, explicit
dim / extent, in a Scene: the flux, resampling, [grid] resample = false, a parameter file, a wind + CME
scene FITS round trip (also next to a planet), the method chosen: spatial for the CME), and the CmeBlob template of
tests/extra_sources.py (its total, its symmetry).
The comparisons with astro_toolbox's cme_imaging itself run only if it is found at TOOLBOX.
Grids of 41 and 42 px (the default CME needs at least 39 px, see the small-grid test).
"""
import io
import os
import sys
import contextlib
import warnings
import numpy as np
import pytest
from astropy import units as u
from extra_sources import CmeBlob
from toycoronagraph import params as P
from toycoronagraph import selection
from toycoronagraph.scene import Scene, Star, Grid, resample
from toycoronagraph.sources.cme import (Cme, cme_mass_fraction, cme_contrast, model_arguments, load_cme_imaging,
                                        contrast_of, star_photons)
from toycoronagraph.sources.images import SceneFits
from toycoronagraph.sources.thomson import Wind
from toycoronagraph.sources.planets import Planet
from toycoronagraph.observe import coronagraph_for

TOOLBOX = "/Users/yu/Downloads/astro_toolbox-main"
HAVE_TOOLBOX = os.path.exists(os.path.join(TOOLBOX, "cme_imaging.py"))
needs_toolbox = pytest.mark.skipif(not HAVE_TOOLBOX, reason="astro_toolbox not found at %s" % TOOLBOX)

SIZES = [41, 42]                       # odd (the model box is the grid) and even (the model is one pixel larger)
WAVELEN, DIAMETER = 0.8e-6, 2.4
PIXSCALE = WAVELEN/DIAMETER*P.RAD_TO_ARCSEC/16    # arcsec (the stage-6 sampling)
DISTANCE, RADIUS = 3.22, 0.73
EPS_ERI = Star(distance=DISTANCE, radius=RADIUS)                       # contrast scenes
EPS_ERI_T = Star(distance=DISTANCE, radius=RADIUS, temperature=5080)   # Jy scenes (cme_imaging needs T)
# CGS constants written out (CODATA 2018 / IAU): Thomson cross-section, atomic mass unit, au, solar radius
SIGMA_T, M_U, AU, R_SUN = 6.6524587321e-25, 1.66053906660e-24, 1.495978707e13, 6.957e10


def grid_scene(n, star=EPS_ERI, **settings):
    """A scene of n px of PIXSCALE arcsec at 0.8 um, 2.4 m, around eps Eri (contrast unless star has T)."""
    return Scene(wavelength=WAVELEN, diameter=DIAMETER, size=n, pixel_scale=PIXSCALE, star=star, **settings)


def moments(fraction):
    """Mass-weighted (mean x, mean y, rms y, mean radius) of a fraction map, in pixels from the star at dim // 2."""
    c = fraction.shape[0]//2
    y, x = np.indices(fraction.shape) - c
    return ((fraction*x).sum(), (fraction*y).sum(), np.sqrt((fraction*y**2).sum()), (fraction*np.hypot(x, y)).sum())


def expected_contrast(fraction, mass_g, pixel_au, radius_au, mu_e=1.0):
    """N_e sigma_T / (4 pi r^2) per pixel from the constants written out (0 for centres inside the star)."""
    n = fraction.shape[0]
    offsets = (np.arange(n) - n//2)*pixel_au*AU
    r = np.hypot(*np.meshgrid(offsets, offsets))
    electrons = fraction*mass_g/(mu_e*M_U)
    out = np.zeros_like(fraction)
    outside = r > radius_au*AU
    out[outside] = electrons[outside]*SIGMA_T/(4*np.pi*r[outside]**2)
    return out


class FakeCmeImagingStar:
    """
    A stand-in for a cme_imaging.Star after particle_to_photon_flux() (the attributes Cme reads).

    The star emits Q photons/s; its flux grid is Q / (4 pi r^2) (NaN on the star pixel, as in
    cme_imaging); grid_cme_photons is chosen so that the CME's contrast per pixel is `pattern`
    (NaN on the star pixel). With pi_short, luminosity / photon_energy is Q / pi, as in cme_imaging.

    Args:
        pattern (np.ndarray): (dim, dim) contrast per pixel, dim odd, the star at dim // 2.
        pix_au (float): Pixel size in au.
        distance_pc (float): Distance in pc.
        pi_short (bool): Make luminosity / photon_energy pi times smaller than Q.
        computed (bool): False: grid_cme_photons is None until particle_to_photon_flux() is called.
    """
    Q = 1e45  # photons per second

    def __init__(self, pattern, pix_au=0.02, distance_pc=DISTANCE, pi_short=True, computed=True):
        dim = pattern.shape[0]
        self.distance, self.radius, self.temp = distance_pc*u.pc, RADIUS*u.R_sun, 5080*u.K
        self.pix_res, self.linear_extent = pix_au*u.au, dim*pix_au*u.au
        offsets = (np.arange(dim) - dim//2)*pix_au
        r = np.hypot(*np.meshgrid(offsets, offsets))
        self.grid_distances = r*u.au
        with np.errstate(divide="ignore"):
            self.grid_stellar_photon_fluxes = np.where(r > 0, self.Q/(4*np.pi*np.where(r > 0, r, 1)**2), np.nan)/u.s/u.au**2
        self._photons = np.where(r > 0, pattern, np.nan)*self.Q/pix_au**2/u.s/u.au**2
        self.grid_cme_photons = self._photons if computed else None
        self.photon_energy = 1.0*u.erg
        self.luminosity = self.Q*u.erg/u.s/(np.pi if pi_short else 1.0)
        self.calls = 0

    def particle_to_photon_flux(self):
        """Fill grid_cme_photons (counted in self.calls)."""
        self.calls += 1
        self.grid_cme_photons = self._photons


def pattern(dim=21, seed=1):
    """A reproducible positive (dim, dim) contrast pattern (pixel values 1e-10 .. 2e-10)."""
    return 1e-10*(1 + np.random.default_rng(seed).random((dim, dim)))


@pytest.fixture
def cme_imaging(monkeypatch):
    """astro_toolbox's cme_imaging module (sys.path restored after the test; skipped without the toolbox)."""
    if not HAVE_TOOLBOX:
        pytest.skip("astro_toolbox not found at %s" % TOOLBOX)
    monkeypatch.setattr(sys, "path", list(sys.path))
    return load_cme_imaging(TOOLBOX)


# ---------------------------------------------------------------- cme_mass_fraction: the ported geometry

@pytest.mark.parametrize("dim, dR_factor, horizontal_height_factor", [(101, 0.6, 0.0), (61, 0.9, 0.0), (101, 0.45, 1.0)])
def test_thick_shell_gives_a_fraction(dim, dR_factor, horizontal_height_factor):
    """A thick shell, whose edges skip columns (cme_imaging raises TypeError there with NumPy 2), is a valid fraction."""
    f = cme_mass_fraction(dim, dR_factor=dR_factor, horizontal_height_factor=horizontal_height_factor)
    assert np.all(np.isfinite(f)) and f.min() >= 0 and f.sum() == pytest.approx(1, rel=1e-12)


@pytest.mark.parametrize("dim", [41, 45, 61])
def test_mass_fraction_shape_dtype_and_normalisation(dim):
    """An odd dim gives a (dim, dim) float64 array, finite, non-negative, summing to 1."""
    f = cme_mass_fraction(dim)
    assert f.shape == (dim, dim) and f.dtype == np.float64
    assert np.all(np.isfinite(f)) and f.min() >= 0 and f.sum() == pytest.approx(1, rel=1e-12)


@pytest.mark.parametrize("dim", [40, 44])
def test_even_dim_is_made_odd(dim):
    """An even dim is made odd (+1), as in cme_imaging: the same array as for dim + 1."""
    f = cme_mass_fraction(dim)
    assert f.shape == (dim+1, dim+1) and np.array_equal(f, cme_mass_fraction(dim+1))


@pytest.mark.parametrize("settings", [{}, dict(dR_factor=0.35, vertical_height_factor=0.6, horizontal_height_factor=0.5),
                                      dict(smear=True, smear_sig=2.0)])
def test_mirror_symmetric_about_the_star_row(settings):
    """The upper half is mirrored to the lower half about the star's row (dim // 2): flipping the rows changes nothing."""
    f = cme_mass_fraction(45, **settings)
    assert np.allclose(f[::-1], f, rtol=1e-14, atol=0)
    assert moments(f)[1] == pytest.approx(0, abs=1e-12)


def test_points_along_plus_x():
    """The CME points along +x: no mass left of the star's column, a positive mean x (rows = y, columns = x)."""
    f = cme_mass_fraction(45)
    assert not f[:, :45//2].any() and moments(f)[0] > 0
    smeared = cme_mass_fraction(45, smear=True)
    assert moments(smeared)[0] > 0 and smeared[:, :22].sum() < 0.1


def test_horizontal_height_factor_moves_the_shell_outwards():
    """horizontal_height_factor > 0 moves the mass away from the star (larger mean radius)."""
    assert moments(cme_mass_fraction(61, horizontal_height_factor=0.5))[3] > moments(cme_mass_fraction(61))[3]


def test_vertical_height_factor_stretches_the_shell_vertically():
    """vertical_height_factor > 1 spreads the mass in y (larger rms y), < 1 squeezes it."""
    rms = lambda v: moments(cme_mass_fraction(61, vertical_height_factor=v))[2]
    assert rms(0.5) < rms(1.0) < rms(1.5)


def test_smear_keeps_the_total_and_spreads_the_mass():
    """smear = true: a Gaussian smoothing; still summing to 1, on more pixels, with a lower peak."""
    sharp, smooth = cme_mass_fraction(45), cme_mass_fraction(45, smear=True, smear_sig=2.0)
    assert smooth.sum() == pytest.approx(1, rel=1e-12)
    assert (smooth > 0).sum() > (sharp > 0).sum() and smooth.max() < sharp.max()


@pytest.mark.parametrize("dR_factor", [0, 1, 1.2, -0.1])
def test_dR_factor_must_be_between_0_and_1(dR_factor):
    """dR_factor outside (0, 1): a ParamError saying so (cme_imaging only asserts dR_factor < 1)."""
    with pytest.raises(P.ParamError, match="between 0 and 1"):
        cme_mass_fraction(45, dR_factor=dR_factor)


@pytest.mark.parametrize("dim", [21, 31])
def test_small_grid_gives_a_fraction_or_a_clear_error(dim):
    """A box too small for the shell's layers: either a valid fraction (finite, summing to 1) or a ParamError."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            f = cme_mass_fraction(dim)
        except P.ParamError:
            return
    assert np.all(np.isfinite(f)) and f.sum() == pytest.approx(1, rel=1e-12)


def test_cme_on_a_small_grid_is_finite_or_a_clear_error():
    """A CME on a 31 px scene: either finite values or a ParamError, never a NaN image."""
    scene = grid_scene(31)
    scene.add(Cme())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            images = scene.render(verbose=False)
        except P.ParamError:
            return
    assert np.all(np.isfinite(images["cme"]))


# ---------------------------------------------------------------- cme_contrast

def test_cme_contrast_is_n_e_sigma_t_over_4_pi_r2():
    """Each pixel: (fraction x mass / (mu_e m_u)) sigma_T / (4 pi r^2), r the pixel centre's distance from the star."""
    fraction = pattern(11, seed=3)
    fraction /= fraction.sum()
    got = cme_contrast(fraction, 2e19, 0.05, 0.004)
    assert got.shape == (11, 11) and got.dtype == np.float64
    assert np.allclose(got, expected_contrast(fraction, 2e19, 0.05, 0.004), rtol=1e-8, atol=0)


def test_cme_contrast_masks_the_star_without_warnings():
    """Pixels whose centre is inside the star are 0 (the star pixel too: no division by zero, no NaN, no warning)."""
    fraction = np.full((11, 11), 1/121)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got = cme_contrast(fraction, 1e19, 0.05, 0.11)   # the star covers the pixels within 2.2 px of its centre
    r = np.hypot(*np.meshgrid(np.arange(11)-5, np.arange(11)-5))
    assert np.all(np.isfinite(got)) and np.all(got[r*0.05 <= 0.11] == 0) and np.all(got[r*0.05 > 0.11] > 0)


def test_cme_contrast_scales_with_mass_and_mu_e():
    """The contrast is proportional to the mass and to 1 / mu_e."""
    fraction = cme_mass_fraction(41)
    base = cme_contrast(fraction, 1e19, 0.01, 0.003)
    assert np.allclose(cme_contrast(fraction, 3e19, 0.01, 0.003), 3*base, rtol=1e-14, atol=0)
    assert np.allclose(cme_contrast(fraction, 1e19, 0.01, 0.003, mu_e=1.17), base/1.17, rtol=1e-14, atol=0)


# ---------------------------------------------------------------- helpers for cme_imaging itself

def test_model_arguments_reads_quantities_from_text():
    """Text with a number and a unit becomes a Quantity; numbers, words and booleans are passed as they are."""
    args = model_arguments({"mass": "1e19 g", "speed": "500 km/s", "dim": 61, "halo": True, "model": "halo"})
    assert args["mass"] == 1e19*u.g and args["speed"] == 500*u.km/u.s
    assert args["dim"] == 61 and args["halo"] is True and args["model"] == "halo"


def test_load_cme_imaging_without_the_file_is_a_file_not_found_error(in_tmp, monkeypatch):
    """toolbox_path without cme_imaging.py: FileNotFoundError naming the file and the path."""
    monkeypatch.setattr(sys, "path", list(sys.path))
    with pytest.raises(FileNotFoundError, match="cme_imaging.py"):
        load_cme_imaging(str(in_tmp))


def test_contrast_of_a_cme_imaging_star():
    """contrast_of: grid_cme_photons x pix_res^2 / the star's photons from its flux grid (NaN kept on the star)."""
    p = pattern()
    model = FakeCmeImagingStar(p)
    assert star_photons(model) == pytest.approx(FakeCmeImagingStar.Q, rel=1e-12)
    got = contrast_of(model)
    assert got.shape == (21, 21) and np.isnan(got[10, 10]) and np.isnan(got).sum() == 1
    ok = np.isfinite(got)
    assert np.allclose(got[ok], p[ok], rtol=1e-12, atol=0)


# ---------------------------------------------------------------- Cme.from_object

def test_from_object_provides_the_star_and_the_grid():
    """from_object: no requirements; inspect() gives the object's distance, pixel scale, radius, temperature and size."""
    source = Cme.from_object(FakeCmeImagingStar(pattern(), pix_au=0.02, distance_pc=3.0))
    assert source.requires == () and Cme.requires == ("star.distance", "star.radius")
    found = {key: value for key, value, origin in source.inspect()}
    assert all(origin == "the cme_imaging Star" for _, _, origin in source.inspect())
    assert found["star.distance"] == pytest.approx(3.0) and found["grid.pixel_scale"] == pytest.approx(0.02/3.0)
    assert found["star.radius"] == pytest.approx(RADIUS) and found["star.temperature"] == pytest.approx(5080)
    assert found["grid.size"] == 21 and Cme().inspect() == []


def test_from_object_calls_particle_to_photon_flux_when_needed():
    """from_object computes the photon grids of a model that has none yet, and leaves a computed one alone."""
    fresh, done = FakeCmeImagingStar(pattern(), computed=False), FakeCmeImagingStar(pattern())
    Cme.from_object(fresh)
    Cme.from_object(done)
    assert fresh.calls == 1 and fresh.grid_cme_photons is not None and done.calls == 0


@pytest.mark.parametrize("pi_short", [True, False])
def test_from_object_renders_its_contrast_with_notes(pi_short):
    """render(): the object's contrast on its own pixel scale, NaN set to 0 with a note; a note when its
    luminosity lacks the factor pi."""
    p = pattern()
    model = FakeCmeImagingStar(p, pix_au=0.02, distance_pc=3.0, pi_short=pi_short)
    img = Cme.from_object(model).render(Grid(21, 0.02/3.0, WAVELEN, DIAMETER), EPS_ERI)
    assert img.unit == "contrast_pixel" and img.pixel_scale == pytest.approx(0.02/3.0, rel=1e-12)
    assert img.data[10, 10] == 0 and np.allclose(np.delete(img.data.ravel(), 220), np.delete(p.ravel(), 220), rtol=1e-12)
    assert any("1 pixel(s) inside the star" in note for note in img.notes)
    assert any("factor pi" in note for note in img.notes) == pi_short


def test_from_object_in_a_scene_takes_its_grid():
    """A scene with only Cme.from_object(...) takes its size, pixel scale and star from the object."""
    scene = Scene(wavelength=WAVELEN, diameter=DIAMETER)
    scene.add(Cme.from_object(FakeCmeImagingStar(pattern(), pix_au=0.02, distance_pc=3.0)))
    scene.render(verbose=False)
    assert scene.grid.size == 21 and scene.grid.pixel_scale == pytest.approx(0.02/3.0, rel=1e-12)
    assert scene.params.origin("grid.size") == "the cme_imaging Star" and scene.star.distance == pytest.approx(3.0)
    contrast = scene.images["cme"]*scene.grid.pixel_area/scene.star_flux     # Jy/arcsec^2 back to contrast
    assert contrast[10, 10] == 0 and contrast.sum() == pytest.approx(pattern().sum()-pattern()[10, 10], rel=1e-10)


@pytest.mark.parametrize("distance, on_conflict, outcome", [(3.2, "error", "error"), (3.2, "warn", "warning"),
                                                            (3.21, "error", "ok")])
def test_distance_conflict_with_the_object(distance, on_conflict, outcome):
    """The object's distance vs [star] distance = 3.22: beyond 0.5% an error (or a note with on_conflict = "warn")."""
    scene = Scene(wavelength=WAVELEN, diameter=DIAMETER, star=Star(distance=DISTANCE), checks={"on_conflict": on_conflict})
    scene.add(Cme.from_object(FakeCmeImagingStar(pattern(), distance_pc=distance)))
    if outcome == "error":
        with pytest.raises(P.ParamError, match="distance"):
            scene.render(verbose=False)
        return
    scene.render(verbose=False)
    assert any("star.distance" in note for note in scene.notes) == (outcome == "warning")


# ---------------------------------------------------------------- the Cme module

def test_cme_declares_its_family_requirements_and_defaults():
    """Cme: a contrast ("relative") needing the star's distance and radius; the cme_imaging defaults; symmetry measured."""
    cme = Cme()
    assert cme.kind == "cme" and cme.family_of() == "relative" and cme.requires == ("star.distance", "star.radius")
    assert cme.settings == dict(mass=1e19, dR_factor=0.2, vertical_height_factor=1.0, horizontal_height_factor=0.0,
                                smear=False, smear_sig=3.0, mu_e=1.0, extent="auto", dim="auto", toolbox_path="",
                                cme_args={})
    assert cme.symmetric() is None and cme.missing() == []


def test_toolbox_path_also_requires_the_temperature_without_changing_the_class():
    """With toolbox_path the instance also needs [star] temperature; the class and other instances do not."""
    assert Cme(toolbox_path="somewhere").requires == ("star.distance", "star.radius", "star.temperature")
    assert Cme.requires == Cme().requires == ("star.distance", "star.radius")


def test_mass_accepts_text_with_units():
    """mass may be text with a unit, converted to g."""
    assert Cme(mass="2e19 g").settings["mass"] == pytest.approx(2e19) and Cme(mass="1 kg").settings["mass"] == 1000


@pytest.mark.parametrize("settings, words", [
    (dict(mass=0), "must be > 0"),
    (dict(dR_factor=0), "must be > 0"),
    (dict(dR_factor=1.5), "must be <= 1"),
    (dict(dim=2), "must be >= 3"),
    (dict(dim=41.5), "not an integer"),
    (dict(extent=0), "must be > 0"),
    (dict(smear="yes"), "not true/false"),
    (dict(masss=1e19), "did you mean mass"),
])
def test_invalid_settings_are_param_errors(settings, words):
    """Out-of-range or mistyped values and unknown names (with a suggestion): ParamError at creation."""
    with pytest.raises(P.ParamError, match=words):
        Cme(**settings)


def test_dR_factor_of_one_stops_at_render():
    """dR_factor = 1 passes the parameter check (max 1) but the geometry refuses it with a clear ParamError."""
    cme = Cme(dR_factor=1)
    with pytest.raises(P.ParamError, match="between 0 and 1"):
        cme.render(Grid(41, PIXSCALE, WAVELEN, DIAMETER), EPS_ERI)


def test_render_on_an_odd_grid_is_the_grid():
    """dim and extent "auto" on an odd grid: the model box is the grid (same size and pixel scale), N_e sigma_T/(4 pi r^2)."""
    img = Cme(mass=2e19).render(Grid(41, PIXSCALE, WAVELEN, DIAMETER), EPS_ERI)
    pixel_au = PIXSCALE*DISTANCE
    assert img.unit == "contrast_pixel" and img.data.shape == (41, 41) and img.data.dtype == np.float64
    assert img.pixel_scale == pytest.approx(PIXSCALE, rel=1e-12) and "odd number" not in img.origin
    expected = expected_contrast(cme_mass_fraction(41), 2e19, pixel_au, RADIUS*R_SUN/AU)
    assert np.allclose(img.data, expected, rtol=1e-8, atol=0)


def test_render_on_an_even_grid_is_one_pixel_larger():
    """On an even grid the model box is made odd (one pixel more, same pixel scale) and its origin says so."""
    img = Cme().render(Grid(42, PIXSCALE, WAVELEN, DIAMETER), EPS_ERI)
    assert img.data.shape == (43, 43) and img.pixel_scale == pytest.approx(PIXSCALE, rel=1e-12)
    assert "odd number of pixels" in img.origin


@pytest.mark.parametrize("dim", [41, 40])
def test_explicit_dim_and_extent(dim):
    """dim and extent set: a box of dim pixels (made odd) across extent au, so extent / dim / d arcsec per pixel."""
    img = Cme(dim=dim, extent=0.5).render(Grid(45, PIXSCALE, WAVELEN, DIAMETER), EPS_ERI)
    odd = dim + (dim % 2 == 0)
    assert img.data.shape == (odd, odd) and img.pixel_scale == pytest.approx(0.5/odd/DISTANCE, rel=1e-12)


@pytest.mark.parametrize("n", SIZES)
def test_cme_is_measured_not_symmetric(n):
    """The CME's symmetry is not known by its module (None) and measures as not symmetric."""
    scene = grid_scene(n)
    scene.add(Cme())
    scene.render(verbose=False)
    symmetric, why = selection.image_symmetry(scene.images["cme"])
    assert Cme().symmetric() is None and symmetric is False and "not symmetric" in why


# ---------------------------------------------------------------- in a scene

def test_scene_on_an_odd_grid_keeps_the_model_and_its_flux():
    """Odd grid: the scene's image is the model's (no resampling), summing to sum N_e sigma_T / (4 pi r^2)."""
    scene = grid_scene(41)
    scene.add(Cme(mass=2e19))
    scene.render(verbose=False)
    assert scene.internal_unit == "contrast_pixel" and not any("resampled" in note for note in scene.notes)
    expected = expected_contrast(cme_mass_fraction(41), 2e19, PIXSCALE*DISTANCE, RADIUS*R_SUN/AU)
    assert scene.images["cme"].sum() == pytest.approx(expected.sum(), rel=1e-8)


@pytest.mark.parametrize("n", SIZES)
def test_twice_the_mass_twice_the_light(n):
    """The image is proportional to the mass (also after resampling on an even grid)."""
    images = []
    for mass in (1e19, 2e19):
        scene = grid_scene(n)
        scene.add(Cme(mass=mass))
        images.append(scene.render(verbose=False)["cme"])
    assert np.allclose(images[1], 2*images[0], rtol=1e-12, atol=0)


def test_scene_in_jy_is_the_contrast_times_the_star_flux():
    """With the star's temperature the scene is in Jy/arcsec^2: contrast x the star's flux / the pixel area."""
    contrast, jy = grid_scene(41), grid_scene(41, star=EPS_ERI_T)
    for scene in (contrast, jy):
        scene.add(Cme())
        scene.render(verbose=False)
    assert jy.internal_unit == "Jy/arcsec^2"
    assert np.allclose(jy.images["cme"], contrast.images["cme"]*jy.star_flux/jy.grid.pixel_area, rtol=1e-12, atol=0)


def test_scene_on_an_even_grid_resamples_with_the_flux_conserved():
    """Even grid: the 43 px model is resampled onto the 42 px grid (flux conserved), with a note naming the CME."""
    scene = grid_scene(42)
    scene.add(Cme())
    scene.render(verbose=False)
    model = Cme().render(scene.grid, scene.star)
    assert any(note.startswith("Note: resampled") and "CME (" in note for note in scene.notes)
    assert np.allclose(scene.images["cme"], resample(model.data, model.pixel_scale, PIXSCALE, 42), rtol=1e-12, atol=0)
    assert scene.images["cme"].sum() <= model.data.sum()   # only the half-pixel border outside the grid is lost


@pytest.mark.parametrize("n, dim, words", [(41, "auto", None), (41, 61, "61 px"), (42, "auto", "odd grid size")])
def test_strict_grid(n, dim, words):
    """[grid] resample = false: "auto" runs on an odd grid; another dim, or an even grid, stops naming the CME."""
    scene = grid_scene(n, resample=False)
    scene.add(Cme(**({} if dim == "auto" else {"dim": dim})))
    if words is None:
        assert set(scene.render(verbose=False)) == {"cme"}
        return
    with pytest.raises(P.ParamError) as err:
        scene.render(verbose=False)
    assert "CME (" in str(err.value) and words in str(err.value)


@pytest.mark.parametrize("n", SIZES)
def test_scene_fits_round_trip_of_a_wind_and_a_cme(n, in_tmp):
    """A wind + CME scene written with Scene.write reads back (kind scene_fits) as the same components on the same grid."""
    scene = grid_scene(n)
    scene.add(Wind(mdot=6e-13, v=400))
    scene.add(Cme())
    scene.render(verbose=False)
    path = str(in_tmp/"wind_cme.fits")
    scene.write(path)
    back = Scene(star=EPS_ERI, diameter=DIAMETER)
    back.add(SceneFits(file=path))
    back.render(verbose=False)
    assert list(back.images) == ["wind", "cme"] and back.grid == scene.grid
    assert back.internal_unit == scene.internal_unit == "contrast_pixel"
    for name in scene.images:
        assert np.allclose(back.images[name], scene.images[name], rtol=1e-12, atol=0), name


def test_scene_fits_read_back_combines_with_a_planet(in_tmp):
    """A wind + CME scene FITS read back next to a planet: the images unchanged, the planet a point source."""
    scene = grid_scene(41)
    scene.add(Wind(mdot=6e-13, v=400))
    scene.add(Cme())
    scene.render(verbose=False)
    scene.write(str(in_tmp/"wind_cme.fits"))
    back = Scene(star=EPS_ERI, diameter=DIAMETER)
    back.add(SceneFits(file=str(in_tmp/"wind_cme.fits")))
    back.add(Planet(mode="polar", pos=[0.3*41*PIXSCALE, 35], brightness=1e-6, brightness_unit="contrast"))
    back.render(verbose=False)
    assert list(back.images) == ["wind", "cme"] and list(back.points) == ["planet"]
    for name in scene.images:
        assert np.array_equal(back.images[name], scene.images[name]), name


def test_cme_blob_template_total_is_the_point_source_value():
    """The CmeBlob template (tests/extra_sources.py): a small blob well inside the field sums to N_e sigma_T / (4 pi r^2)
    at its centre (electrons = mass / m_p), to 1%."""
    scene = grid_scene(61)                      # half width 30 px = 0.42 au
    scene.add(CmeBlob(mass=1e19, r=0.3, pa=30, width=0.015))   # 21.7 px from the star, sigma 1.1 px
    scene.render(verbose=False)
    point = 1e19/1.67262192369e-24*SIGMA_T/(4*np.pi*(0.3*AU)**2)
    assert scene.images["cme_blob"].sum() == pytest.approx(point, rel=0.01)


@pytest.mark.parametrize("n", SIZES)
def test_cme_goes_to_the_spatial_method_a_centred_blob_and_a_wind_to_the_circular_one(n, in_tmp):
    """auto: the measured, not symmetric CME needs the spatial method; a centred blob and the wind are known symmetric
    by their modules (circular). Only the method is chosen here (no PSF is computed)."""
    scene = grid_scene(n, compute={"psf_cache_dir": str(in_tmp/"cache"), "num_cores": 2})
    sources = {"wind": Wind(mdot=6e-13, v=400), "cme": Cme(mass=1e20), "cme_blob": CmeBlob(r=0, width=0.05, mass=1e19)}
    for source in sources.values():
        scene.add(source)
    scene.render(verbose=False)
    coro = coronagraph_for(scene.grid, "vortex:6", scene.params)[1]
    methods = {}
    for name, source in sources.items():
        known = source.symmetric()
        symmetry = None if known is None else (known, "known")
        methods[name] = coro.choose(scene.images[name], symmetry, P.AUTO, P.AUTO)[0]
    assert methods == {"wind": "circular", "cme": "spatial", "cme_blob": "circular"}
    assert CmeBlob(r=0.1).symmetric() is False and not os.path.exists(str(in_tmp/"cache"))


def test_parameter_file_equals_python(in_tmp):
    """A cme from a parameter file (mass as text with a unit) is the same image as in Python."""
    path = in_tmp/"cme.toml"
    path.write_text('[observation]\nwavelength = %r\ndiameter = %r\n[grid]\nsize = 41\npixel_scale = %r\n'
                    '[star]\ndistance = %r\nradius = %r\n[[source]]\nkind = "cme"\nmass = "2e19 g"\n'
                    % (WAVELEN, DIAMETER, PIXSCALE, DISTANCE, RADIUS))
    from_file = Scene.from_file(str(path))
    from_file.render(verbose=False)
    scene = grid_scene(41)
    scene.add(Cme(mass=2e19))
    scene.render(verbose=False)
    assert np.allclose(from_file.images["cme"], scene.images["cme"], rtol=1e-12, atol=0)


# ---------------------------------------------------------------- astro_toolbox's cme_imaging itself (if found)

@needs_toolbox
@pytest.mark.parametrize("dim, settings", [
    (41, {}), (44, {}), (45, dict(dR_factor=0.35, vertical_height_factor=0.6, horizontal_height_factor=0.5)),
    (127, dict(dR_factor=0.1, vertical_height_factor=1.5)), (41, dict(smear=True, smear_sig=2.0)),
])
def test_ported_geometry_is_cme_imaging_bit_for_bit(dim, settings, cme_imaging):
    """cme_mass_fraction() = cme_imaging's CME.build_cme_normalized_grid exactly, for several sizes and settings."""
    c = dict(settings)
    smear, sig = c.pop("smear", False), c.pop("smear_sig", 3.0)
    with contextlib.redirect_stdout(io.StringIO()):
        cme = cme_imaging.CME(dim=dim)
        cme.build_cme_normalized_grid(smear=smear, smear_sig=sig, **c)
    assert np.array_equal(cme_mass_fraction(dim, smear=smear, smear_sig=sig, **c), cme.grid_norm)


@needs_toolbox
def test_toolbox_path_equals_the_built_in_cme(cme_imaging):
    """Cme(toolbox_path, cme_args) built by cme_imaging = the built-in CME on the same grid (to round-off)."""
    built_in, external = grid_scene(41, star=EPS_ERI_T), grid_scene(41, star=EPS_ERI_T)
    built_in.add(Cme(mass=2e19))
    external.add(Cme(toolbox_path=TOOLBOX, cme_args={"mass": "2e19 g"}))
    built_in.render(verbose=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        external.render(verbose=False)
    assert np.abs(external.images["cme"]-built_in.images["cme"]).max() < 1e-12*built_in.images["cme"].max()
    assert any("factor pi" in note for note in external.notes)


@needs_toolbox
def test_cme_imaging_on_a_small_grid_gives_a_clear_error(cme_imaging):
    """cme_imaging's model box too small to hold the shell (all NaN): a ParamError, not a silent zero image."""
    scene = grid_scene(31, star=EPS_ERI_T)
    scene.add(Cme(toolbox_path=TOOLBOX, cme_args={"mass": "2e19 g"}))
    with pytest.raises(P.ParamError, match="NaN everywhere"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scene.render(verbose=False)


@needs_toolbox
def test_cme_args_without_a_unit_get_a_hint(cme_imaging):
    """A cme_args mass without its unit: a ParamError showing how to write it."""
    scene = grid_scene(41, star=EPS_ERI_T)
    scene.add(Cme(toolbox_path=TOOLBOX, cme_args={"mass": 2e19}))
    with pytest.raises(P.ParamError, match='mass = "1e19 g"'):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scene.render(verbose=False)
