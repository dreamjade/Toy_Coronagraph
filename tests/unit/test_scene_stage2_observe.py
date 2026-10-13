"""
Unit tests for stage F2 (tests/scene_stage2_observe.py): a scene through one or more coronagraphs.

Covers toycoronagraph/observe.py: observe() (combine "separate" / "together", the method chosen per
component, components=[...], iwa_ignore, several coronagraphs, planets with their exact PSFs placed rows-up),
Coronagraph / coronagraph_for (shared objects and caches), and Result (final, in_unit, contributions,
contrast, save (folder names, FITS headers), load: also for points-only scenes, which need SCNSIZE, and for
contrast scenes and star_flux), [grid] resample = false (stops before any PSF), and Target vs observe().
Tiny scenes: 16 and 17 px (even and odd) at 4 px per lambda/D, with the template modules
of tests/extra_sources.py (Ring: symmetric, Jy; CmeBlob: not symmetric, contrast) and planets.

The PSFs are computed once per module in a pytest temporary folder (num_cores = 2); every test runs
in its own temporary folder (in_tmp), so run folders never land in the project.
"""
import os
import json
import numpy as np
import pytest
from astropy.io import fits
from extra_sources import Ring, CmeBlob
from toycoronagraph import params as P
from toycoronagraph import units as U
from toycoronagraph.scene import Scene, Star, Grid, point_labels
from toycoronagraph.sources.planets import Planet
from toycoronagraph.sources.images import FitsImage
from toycoronagraph.observe import observe, Result, Coronagraph, coronagraph_for, _digest
from toycoronagraph.optics import Optics
from toycoronagraph.psf import exact_point_sources

SIZES = [16, 17]                       # even (star between four pixels) and odd (star on a pixel)
SPATIAL_N = 16                         # the spatial method's PSF model is slow to build: tests needing it use one size
WAVELEN, DIAMETER, SAMPLING = 0.8e-6, 2.4, 4
PIXSCALE = WAVELEN/DIAMETER*P.RAD_TO_ARCSEC/SAMPLING     # arcsec: 4 px per lambda/D
DISTANCE = 3.22
STAR = Star(distance=DISTANCE, radius=0.73, temperature=5080)   # flux known: scenes in Jy/arcsec^2
CONTRAST_STAR = Star(distance=DISTANCE, radius=0.73)            # no flux: scenes in contrast_pixel
PLANET_JY = 2e-4


@pytest.fixture(scope="module")
def cache(tmp_path_factory):
    """One PSF cache folder for this module (in pytest's temporary area), so each PSF is computed once."""
    return str(tmp_path_factory.mktemp("psf_cache"))


@pytest.fixture(autouse=True)
def _work_in_tmp(in_tmp):
    """Run every test in its own temporary folder (run folders and any stray file land there)."""
    return in_tmp


def make_scene(n, cache, star=STAR, ring=True, blob=True, planet=True, ring_scale=1.0, blob_method="circular",
               planet_pa=35.0, planet_brightness=PLANET_JY, planet_unit="Jy", **settings):
    """
    A tiny scene on n px at 4 px per lambda/D.

    Args:
        n (int): Grid size.
        cache (str): The PSF cache folder ([compute] psf_cache_dir).
        star (Star): The star (STAR: Jy scene; CONTRAST_STAR: contrast scene if every component is relative).
        ring, blob, planet (bool): Include a Ring (Jy, symmetric), a CmeBlob (contrast, at PA 40 deg) and a
            Planet (at 0.35 n pixels from the star).
        ring_scale (float): The ring component's scale.
        blob_method (str): The blob component's own method: "circular" by default, so that most tests need no
            spatial PSF model (slow to build); "auto" puts it on the spatial method (it is not symmetric).
        planet_pa (float): The planet's position angle in deg (from +x towards +y).
        planet_brightness (float): The planet's brightness, in planet_unit.
        planet_unit (str): "Jy", or "contrast" (a fraction of the star's flux).
        settings: More Scene settings.

    Returns:
        Scene: Not rendered.
    """
    settings = dict(dict(wavelength=WAVELEN, diameter=DIAMETER, size=n, pixel_scale=PIXSCALE, star=star,
                         compute={"psf_cache_dir": cache, "num_cores": 2}), **settings)
    scene = Scene(**settings)
    if ring:
        scene.add(Ring(radius=0.25*n*PIXSCALE, width=0.05*n*PIXSCALE, peak=0.02), scale=ring_scale)
    if blob:
        scene.add(CmeBlob(r=0.2*n*PIXSCALE*DISTANCE, pa=40, width=0.06*n*PIXSCALE*DISTANCE, mass=1e21),
                  method=blob_method)
    if planet:
        scene.add(Planet(mode="polar", pos=[0.35*n*PIXSCALE, planet_pa], brightness=planet_brightness,
                         brightness_unit=planet_unit))
    return scene


def rel_max(a, b):
    """Largest absolute difference of two images, relative to the peak of |b|."""
    return np.abs(a-b).max()/np.abs(b).max()


# ---------------------------------------------------------------- observe(): images, methods, combine

@pytest.mark.parametrize("n", SIZES)
def test_final_images_are_n_by_n_float64_in_the_scene_unit(n, cache):
    """Every final image is an (N, N) float64 array in the scene's internal unit; final() is their sum."""
    scene = make_scene(n, cache)
    res = observe(scene, ["vortex:2"], verbose=False)
    assert res.unit == scene.internal_unit == "Jy/arcsec^2"
    assert list(res.finals) == ["vortex_c2"] and list(res.finals["vortex_c2"]) == ["ring", "cme_blob", "planet"]
    for image in res.finals["vortex_c2"].values():
        assert image.shape == (n, n) and image.dtype == np.float64
    assert np.allclose(res.final(), sum(res.finals["vortex_c2"].values()), rtol=1e-14, atol=0)


@pytest.mark.parametrize("n", SIZES)
def test_observe_renders_the_scene_when_needed(n, cache):
    """observe() renders a scene that is not rendered yet."""
    scene = make_scene(n, cache, blob=False, planet=False)
    assert scene.images is None
    observe(scene, ["vortex:2"], verbose=False)
    assert set(scene.images) == {"ring"}


@pytest.mark.parametrize("n", SIZES)
def test_no_coronagraphs_given_means_a_charge_2_vortex(n, cache):
    """coronagraphs=None and none in the parameters: one charge-2 vortex."""
    res = observe(make_scene(n, cache, blob=False, planet=False), verbose=False)
    assert list(res.finals) == ["vortex_c2"]


@pytest.mark.parametrize("n, method", [(SPATIAL_N, "spatial"), (16, "circular"), (17, "circular")])
def test_together_equals_separate_with_the_same_method(n, method, cache):
    """The coronagraph is linear: with one method for everything, summing before or after it is the same."""
    scene = make_scene(n, cache)
    separate = observe(scene, ["vortex:2"], method=method, verbose=False)
    together = observe(scene, ["vortex:2"], method=method, combine="together", verbose=False)
    assert set(separate.finals["vortex_c2"]) == {"ring", "cme_blob", "planet"}
    assert set(together.finals["vortex_c2"]) == {"total", "planet"}
    assert rel_max(together.final(), separate.final()) < 1e-12
    assert np.array_equal(together.final(component="planet"), separate.final(component="planet"))


def test_auto_separate_gives_each_component_its_method(cache):
    """auto + separate: the symmetric ring on the circular method, the blob on the spatial one, planets exact."""
    res = observe(make_scene(SPATIAL_N, cache, blob_method="auto"), ["vortex:2"], verbose=False)
    ring, blob, planet = (res.methods["vortex_c2", c] for c in ("ring", "cme_blob", "planet"))
    assert ring[0] == "circular" and ring[1] is None and ring[3] == 0.0
    assert "known by the ring module" in ring[2]
    assert blob[0] == "spatial" and blob[1] in ("rotate", "direct") and blob[3] > 0
    assert planet[0] == "exact"


def test_together_is_spatial_when_a_component_is_not_symmetric(cache):
    """auto + together: one component known not to be symmetric makes the sum spatial (component methods ignored)."""
    res = observe(make_scene(SPATIAL_N, cache, planet=False), ["vortex:2"], combine="together", verbose=False)
    method, _, why, _ = res.methods["vortex_c2", "total"]
    assert method == "spatial" and "a component is not" in why


@pytest.mark.parametrize("n", SIZES)
def test_together_is_circular_when_every_component_is_symmetric(n, cache):
    """auto + together: the sum of symmetric components is symmetric, so circular."""
    scene = make_scene(n, cache, blob=False, planet=False)
    scene.add(Ring(radius=0.1*n*PIXSCALE, width=0.05*n*PIXSCALE, peak=0.01))
    res = observe(scene, ["vortex:2"], combine="together", verbose=False)
    method, _, why, _ = res.methods["vortex_c2", "total"]
    assert method == "circular" and "every component is" in why


@pytest.mark.parametrize("combine", ["both", "Together"])
def test_unknown_combine_is_a_param_error(combine, cache):
    """combine must be "separate" or "together"; anything else is a ParamError saying so."""
    with pytest.raises(P.ParamError, match="separate.*together"):
        observe(make_scene(16, cache, blob=False, planet=False), ["vortex:2"], combine=combine, verbose=False)


def test_combine_none_uses_the_method_section(cache):
    """combine=None takes [method] combine (default "separate"); an explicit combine overrides it."""
    scene = make_scene(16, cache, blob=False, planet=False)
    assert observe(scene, ["vortex:2"], verbose=False).combine == "separate"
    scene = make_scene(16, cache, blob=False, planet=False, method={"combine": "together"})
    assert observe(scene, ["vortex:2"], verbose=False).combine == "together"
    assert observe(scene, ["vortex:2"], combine="separate", verbose=False).combine == "separate"


def test_unknown_component_is_named_with_the_available_ones(cache):
    """components=[...] with a name not in the scene: a ParamError naming it and listing the scene's."""
    with pytest.raises(P.ParamError) as err:
        observe(make_scene(16, cache), ["vortex:2"], components=["comet"], verbose=False)
    assert "comet" in str(err.value) and "ring" in str(err.value) and "planet" in str(err.value)


@pytest.mark.parametrize("n", SIZES)
def test_components_alone_add_up_to_the_combined_run(n, cache):
    """Each component alone (components=[name]) equals its part of the combined run; the parts add up to it."""
    scene = make_scene(n, cache)
    full = observe(scene, ["vortex:2"], verbose=False)
    total = 0
    for name in full.finals["vortex_c2"]:
        alone = observe(scene, ["vortex:2"], components=[name], verbose=False)
        assert list(alone.finals["vortex_c2"]) == [name]
        assert np.array_equal(alone.final(), full.final(component=name)), name
        total = total + alone.final()
    assert np.allclose(total, full.final(), rtol=1e-12, atol=0)


@pytest.mark.parametrize("n", SIZES)
def test_two_coronagraphs_in_one_call_equal_two_calls(n, cache):
    """observe(scene, [a, b]) gives each coronagraph the same final images as observe(scene, [a]) and [b]."""
    scene = make_scene(n, cache, blob=False)
    both = observe(scene, ["vortex:2", "vortex:6"], verbose=False)
    assert list(both.finals) == ["vortex_c2", "vortex_c6"]
    for spec, name in [("vortex:2", "vortex_c2"), ("vortex:6", "vortex_c6")]:
        one = observe(scene, [spec], verbose=False)
        assert np.array_equal(one.final(name), both.final(name))
    assert not np.allclose(both.final("vortex_c2"), both.final("vortex_c6"))


def test_empty_run_gives_a_zero_image(cache):
    """A run with no component (components=[]) has an (N, N) zero final image."""
    res = observe(make_scene(16, cache, blob=False, planet=False), ["vortex:2"], components=[], verbose=False)
    assert res.finals["vortex_c2"] == {}
    final = res.final()
    assert isinstance(final, np.ndarray) and final.shape == (16, 16) and not final.any()


def test_empty_run_plots(cache, in_tmp):
    """A run with no component still plots (its total is zero), with log profiles too."""
    res = observe(make_scene(16, cache, blob=False, planet=False), ["vortex:2"], components=[], verbose=False)
    res.plot(path=str(in_tmp/"finals.png"), show=False)
    res.plot_profiles(log=True, path=str(in_tmp/"profiles.png"), show=False)
    assert (in_tmp/"finals.png").exists() and (in_tmp/"profiles.png").exists()


# ---------------------------------------------------------------- Coronagraph objects and caches

def test_digest_depends_on_the_values_only():
    """_digest: 40 hex characters, the same for equal values in float32 or float64, different otherwise."""
    a = np.arange(6.0).reshape(2, 3)
    assert len(_digest(a)) == 40 and _digest(a) == _digest(a.astype(np.float32)) == _digest(a.tolist())
    assert _digest(a) != _digest(a+1e-12)


def test_coronagraph_for_reads_strings_and_tables(cache):
    """"vortex:6", ":6" and {"type": "vortex", "charge": 6} name the same coronagraph object; "vortex" is charge 2."""
    scene = make_scene(16, cache, blob=False, planet=False)
    scene.build()
    table, coro = coronagraph_for(scene.grid, "vortex:6", scene.params)
    assert table.name == coro.name == "vortex_c6" and coro.tag == "vortex_c6_N16_range2_lyot0.95"
    assert coronagraph_for(scene.grid, ":6", scene.params)[1] is coro
    assert coronagraph_for(scene.grid, {"type": "vortex", "charge": 6}, scene.params)[1] is coro
    table2, coro2 = coronagraph_for(scene.grid, "vortex", scene.params)
    assert table2.name == "vortex_c2" and coro2 is not coro
    assert coro.optics.img_pixel == 16 and coro.cache_dir == cache and coro.num_cores == 2


def test_coronagraph_get_is_one_object_per_optics_and_settings(cache):
    """Coronagraph.get returns the same object for the same optics, name and settings, a new one otherwise."""
    optics = Optics.vortex(2, 16, 2.0, 0.95)
    a = Coronagraph.get(optics, "vortex_c2", cache_dir=cache, num_cores=2)
    assert Coronagraph.get(optics, "vortex_c2", num_cores=2, cache_dir=cache) is a
    assert Coronagraph.get(optics, "vortex_c2", cache_dir=cache, num_cores=1) is not a
    assert a.rot_number == 16  # 0 = N rotations


def test_kernel_size_choice_is_cached_on_disk(tmp_path, monkeypatch):
    """The kernel size chosen for an image and the model's error report are saved in the PSF cache folder
    and read back by a new Coronagraph (a new session) without choosing or checking again."""
    scene = make_scene(SPATIAL_N, str(tmp_path), planet=False)
    weight = sum(scene.render(verbose=False).values())
    coro = coronagraph_for(scene.grid, "vortex:2", scene.params)[1]
    K = coro.spatial_model("rotate", weight)["kernel_size"]
    files = list(tmp_path.glob("kernel_vortex_c2_*.json"))
    assert len(files) == 1 and json.loads(files[0].read_text()) == {"kernel_size": K}
    Coronagraph._instances.clear()  # a new session
    assert len(list(tmp_path.glob("error_vortex_c2_*_rotate_K%d_t*.npz" % K))) == 1
    monkeypatch.setattr("toycoronagraph.observe.choose_kernel_size", lambda *a, **k: pytest.fail("chosen again"))
    monkeypatch.setattr("toycoronagraph.psf.spatial_error_report", lambda *a, **k: pytest.fail("checked again"))
    again = coronagraph_for(scene.grid, "vortex:2", scene.params)[1]
    model = again.spatial_model("rotate", weight)
    assert again is not coro and model["kernel_size"] == K
    assert np.array_equal(model["err_median"], coro.spatial_model("rotate", weight)["err_median"], equal_nan=True)


def _same(a, b):
    """Two model entries are equal (arrays: same values, NaNs equal, and dtype; other values: ==)."""
    a, b = np.asarray(a), np.asarray(b)
    return a.dtype == b.dtype and np.array_equal(a, b, equal_nan=a.dtype.kind in "fc")


def _spatial_scene(tmp_path):
    """A scene on SPATIAL_N px, its total image (the kernel-size weight) and its charge-2 vortex."""
    scene = make_scene(SPATIAL_N, str(tmp_path), planet=False)
    weight = sum(scene.render(verbose=False).values())
    return scene, weight, coronagraph_for(scene.grid, "vortex:2", scene.params)[1]


def test_spatial_model_is_cached_on_disk_with_what_it_was_made_from(tmp_path, monkeypatch):
    """The built model is saved (model_<tag>_rotate_K<K>_<16 hex>.npz, its description inside: the optics in full,
    mode, kernel size, error target, radial PSF file, code) and a new session reads it back instead of building
    it: every entry identical, the final image identical, nothing built."""
    from toycoronagraph.spatial import spatial_psf
    scene, weight, coro = _spatial_scene(tmp_path)
    model = coro.spatial_model("rotate", weight)
    K = model["kernel_size"]
    files = list(tmp_path.glob("model_vortex_c2_N%d_*_rotate_K%d_*.npz" % (SPATIAL_N, K)))
    assert len(files) == 1 and not list(tmp_path.glob("*.part"))
    with np.load(files[0]) as f:
        stored = json.loads(f["description"].item())
    assert stored["optics"] == repr(coro.optics) and (stored["mode"], stored["kernel_size"]) == ("rotate", K)
    assert stored["error_target"] == coro.error_target and stored["radial_psfs"][0] == os.path.basename(coro.radial_psfs_file())
    assert len(stored["code"]) == 16
    Coronagraph._instances.clear()  # a new session
    monkeypatch.setattr("toycoronagraph.observe.spatial_psf_calculation", lambda *a, **k: pytest.fail("built again"))
    again = coronagraph_for(scene.grid, "vortex:2", scene.params)[1]
    back = again.spatial_model("rotate", weight)
    assert again is not coro and not [what for what, _ in again.computed if what == "spatial model"]
    assert set(back) == set(model) and all(_same(back[k], model[k]) for k in model)
    assert np.array_equal(spatial_psf(weight, back, report=False), spatial_psf(weight, model, report=False))


@pytest.mark.parametrize("change", ["error target", "radial PSFs recomputed", "code", "another description in the file"])
def test_spatial_model_is_built_again_when_what_it_was_made_from_differs(change, tmp_path, monkeypatch):
    """A cached model is not used if anything it was made from differs: the error target, the radial PSF file
    (recomputed: a new modification time), the code that builds it, or the description stored in the file
    (it is checked, not only the file name); the model is then built again and saved."""
    scene, weight, coro = _spatial_scene(tmp_path)
    K = coro.spatial_model("rotate", weight)["kernel_size"]
    Coronagraph._instances.clear()
    settings = {}
    if change == "error target":
        settings = dict(method={"error_target": 0.03})
    elif change == "radial PSFs recomputed":
        path = coro.radial_psfs_file()
        os.utime(path, ns=(os.stat(path).st_atime_ns, os.stat(path).st_mtime_ns+10**9))
    elif change == "code":
        monkeypatch.setattr("toycoronagraph.observe._code_digest", lambda *modules: "0"*16)
    else:
        [path] = tmp_path.glob("model_*.npz")
        with np.load(path) as f:
            entries = {k: f[k] for k in f.files if k != "description"}
        np.savez(path, description=json.dumps({"optics": "another"}), **entries)
    scene = make_scene(SPATIAL_N, str(tmp_path), planet=False, **settings)
    scene.render(verbose=False)
    again = coronagraph_for(scene.grid, "vortex:2", scene.params)[1]
    model = again.spatial_model("rotate", weight)
    assert [what for what, _ in again.computed if what == "spatial model"] == ["spatial model"]
    if change != "error target":
        assert model["kernel_size"] == K
    if change == "another description in the file":  # rewritten with the right one
        [path] = tmp_path.glob("model_*.npz")
        with np.load(path) as f:
            assert json.loads(f["description"].item())["optics"] == repr(again.optics)


def test_a_changed_component_is_the_only_one_computed_again(cache):
    """After changing the ring's scale only its final image is computed again (the blob's comes from the cache)."""
    scene = make_scene(SPATIAL_N, cache, planet=False)
    scene.build()
    coro = coronagraph_for(scene.grid, "vortex:2", scene.params)[1]   # shared with earlier tests of this module
    start = len(coro.computed)
    res = observe(scene, ["vortex:2"], method="spatial", verbose=False)
    assert res.coronagraphs["vortex_c2"] is coro
    built = [what for what, _ in coro.computed[start:]]
    assert built.count("spatial model") <= 1   # at most one spatial model for both components
    before = len(coro.computed)
    observe(make_scene(SPATIAL_N, cache, planet=False, ring_scale=1.01), ["vortex:2"], method="spatial", verbose=False)
    assert [what for what, _ in coro.computed[before:]] == ["final image"]


# ---------------------------------------------------------------- planets: exact PSFs, rows up, contrast()

@pytest.mark.parametrize("n", SIZES)
def test_planet_gets_its_exact_psf_at_its_pixel(n, cache):
    """A planet's final image is its exact PSF at Grid.pixel(r, theta), its flux per pixel area in Jy/arcsec^2."""
    scene = make_scene(n, cache, ring=False, blob=False)
    res = observe(scene, ["vortex:2"], verbose=False)
    g = scene.grid
    optics = Optics.vortex(2, n, g.visual_range/2, 0.95)
    exact = exact_point_sources(optics, [g.pixel(0.35*n*PIXSCALE, np.radians(35))], [PLANET_JY/g.pixel_area])
    assert np.allclose(res.final(component="planet"), exact, rtol=1e-12, atol=0)


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("pa", [0, 90, 180, 270])
def test_planets_are_placed_rows_up(n, pa, cache):
    """Rows are +y (as in FITS): a planet at PA 90 deg peaks above the star, at 0 deg to its right (star at (N-1)/2)."""
    res = observe(make_scene(n, cache, ring=False, blob=False, planet_pa=pa), ["vortex:2"], verbose=False)
    image = res.final(component="planet")
    row, col = np.unravel_index(np.argmax(image), image.shape)
    c, r, t = (n-1)/2, 0.35*n, np.radians(pa)
    dy, dx = row-c, col-c                                 # pixels from the star; rows = +y
    assert np.hypot(dx-r*np.cos(t), dy-r*np.sin(t)) < 0.5*SAMPLING   # within 0.5 lambda/D of the planet
    assert dx*np.cos(t) + dy*np.sin(t) > 0.25*n                       # on the planet's side (not mirrored)


@pytest.mark.parametrize("n", SIZES)
def test_contrast_gives_each_components_light_in_the_planet_core(n, cache):
    """contrast(): each component summed over the planet's core (pixels above half its peak), times the pixel area."""
    scene = make_scene(n, cache)
    res = observe(scene, ["vortex:2"], verbose=False)
    shares = res.contrast("planet")
    finals = res.finals["vortex_c2"]
    assert set(shares) == {"ring", "cme_blob", "planet"} and all(isinstance(v, float) for v in shares.values())
    core = finals["planet"] > 0.5*finals["planet"].max()
    for name, value in shares.items():
        assert np.isclose(value, finals[name][core].sum()*scene.grid.pixel_area, rtol=1e-12, atol=0)
    # the core holds part of the planet's PLANET_JY Jy (the throughput and the PSF wings take the rest)
    assert 0.1*PLANET_JY < shares["planet"] < PLANET_JY and shares["ring"] > 0 and shares["cme_blob"] > 0


@pytest.mark.parametrize("n", SIZES)
def test_contrast_scene_planet_is_a_fraction_of_the_star(n, cache):
    """In a contrast scene a planet's brightness is a fraction of the star's flux, not multiplied by the pixel area."""
    scene = make_scene(n, cache, star=CONTRAST_STAR, ring=False, blob=False, planet_brightness=1e-6,
                       planet_unit="contrast")
    res = observe(scene, ["vortex:2"], verbose=False)
    assert res.unit == "contrast_pixel" and res.star_flux is None
    assert point_labels(res.points, res.unit, res.grid.pixel_area)["planet"][1] == ["planet 1e-06 of the star"]
    assert 0.1e-6 < res.contrast("planet")["planet"] < 1e-6


# ---------------------------------------------------------------- units: contrast results, star_flux

@pytest.mark.parametrize("n", SIZES)
def test_contrast_result_converts_to_a_brightness_with_the_star_flux(n, cache):
    """A contrast result x the star's flux (a Star or Jy) = the same scene computed in Jy/arcsec^2, and back."""
    contrast = observe(make_scene(n, cache, star=CONTRAST_STAR, ring=False, planet=False), ["vortex:2"],
                       method="circular", verbose=False)
    jy = observe(make_scene(n, cache, ring=False, planet=False), ["vortex:2"], method="circular", verbose=False)
    assert contrast.unit == "contrast_pixel" and jy.unit == "Jy/arcsec^2"
    assert jy.star_flux == pytest.approx(STAR.flux_at(WAVELEN), rel=1e-12)
    for flux in (STAR, jy.star_flux):
        assert np.allclose(contrast.final(unit="Jy/arcsec^2", star_flux=flux), jy.final(), rtol=1e-10, atol=0)
    assert np.allclose(jy.final(unit="contrast_pixel"), contrast.final(), rtol=1e-10, atol=0)
    assert np.allclose(jy.final(unit="Jy/pixel"), jy.final()*jy.grid.pixel_area, rtol=1e-14, atol=0)


def test_contrast_result_without_the_star_flux_cannot_become_a_brightness(cache):
    """Converting a contrast result to Jy without a star flux: a UnitError naming the star's flux."""
    res = observe(make_scene(16, cache, star=CONTRAST_STAR, ring=False, planet=False), ["vortex:2"], verbose=False)
    with pytest.raises(U.UnitError, match="star's flux"):
        res.final(unit="Jy/arcsec^2")
    with pytest.raises(U.UnitError, match="star's flux"):
        res.final(unit="Jy/arcsec^2", star_flux=Star(distance=DISTANCE))


def test_in_unit_none_or_auto_returns_the_same_array():
    """in_unit(image, None / "auto" / the result's unit) returns the image itself, unconverted."""
    res = Result(Grid(16, PIXSCALE, WAVELEN, DIAMETER), "Jy/arcsec^2", STAR, 1.0, "separate", {})
    image = np.ones((16, 16))
    for unit in (None, "auto", "Jy/arcsec^2"):
        assert res.in_unit(image, unit) is image


# ---------------------------------------------------------------- Result on hand-made finals

def hand_made_result():
    """A Result with two coronagraphs and hand-made 4 x 4 finals (zeros in the top row of the total)."""
    res = Result(Grid(4, 0.01, WAVELEN, DIAMETER), "Jy/arcsec^2", STAR, 2.0, "separate", {})
    a = np.zeros((4, 4))
    a[1:, :] = 1.0
    b = np.zeros((4, 4))
    b[2:, :] = 3.0
    res.finals = {"first": {"a": a, "b": b}, "second": {"a": 2*a}}
    return res


def test_final_defaults_to_the_first_coronagraph_and_the_sum():
    """final(): the first coronagraph and the sum of its components by default; one component when named."""
    res = hand_made_result()
    a, b = res.finals["first"]["a"], res.finals["first"]["b"]
    assert np.array_equal(res.final(), a+b) and np.array_equal(res.final("first"), a+b)
    assert np.array_equal(res.final(component="b"), b) and np.array_equal(res.final("second"), 2*a)


def test_final_unknown_names_raise_key_error():
    """final() with an unknown coronagraph or component: KeyError."""
    res = hand_made_result()
    with pytest.raises(KeyError):
        res.final("third")
    with pytest.raises(KeyError):
        res.final(component="c")


def test_contributions_sum_to_one_and_are_zero_where_the_total_is():
    """contributions(): fractions summing to 1 where the total is not 0, and 0 (no NaN) where it is."""
    res = hand_made_result()
    shares = res.contributions()
    total = res.final()
    assert set(shares) == {"a", "b"}
    assert np.allclose(sum(shares.values())[total != 0], 1, rtol=0, atol=1e-15)
    assert all(np.all(s[total == 0] == 0) and np.all(np.isfinite(s)) for s in shares.values())
    assert np.allclose(shares["b"][3], 0.75)


@pytest.mark.parametrize("n", SIZES)
def test_contributions_of_a_run_sum_to_one(n, cache):
    """A real run's contributions add up to 1 wherever the total final image is not 0."""
    res = observe(make_scene(n, cache), ["vortex:2"], verbose=False)
    total = res.final()
    assert np.allclose(sum(res.contributions().values())[total != 0], 1, rtol=0, atol=1e-12)


# ---------------------------------------------------------------- the run folder: save and load

@pytest.mark.parametrize("n", SIZES)
def test_save_writes_the_run_folder(n, cache, in_tmp):
    """save(): params_used.toml, scene.fits, one FITS per coronagraph (FINAL_TOTAL = the components' sum, star at
    CRPIX = (N+1)/2) and the quicklook plots."""
    res = observe(make_scene(n, cache), ["vortex:2"], verbose=False)
    folder = res.save(str(in_tmp/"run"))
    assert folder == res.folder == str(in_tmp/"run")
    assert sorted(os.listdir(folder)) == ["params_used.toml", "quicklook", "scene.fits", "vortex_c2.fits"]
    assert sorted(os.listdir(os.path.join(folder, "quicklook"))) == ["contributions.png", "finals.png"]
    with fits.open(os.path.join(folder, "vortex_c2.fits")) as h:
        assert [hdu.name for hdu in h[1:]] == ["FINAL_TOTAL", "FINAL_RING", "FINAL_CME_BLOB", "FINAL_PLANET"]
        assert np.allclose(h["FINAL_TOTAL"].data, sum(hdu.data for hdu in h[2:]), rtol=1e-12, atol=0)
        p0, ring = h[0].header, h["FINAL_RING"].header
        assert p0["CORONAG"] == "vortex_c2" and p0["COMBINE"] == "separate" and p0["BUNIT"] == "Jy/arcsec^2"
        assert p0["STARFLUX"] == pytest.approx(res.star_flux, rel=1e-12) and p0["CHARGE"] == 2
        assert ring["CRPIX1"] == ring["CRPIX2"] == (n+1)/2 and ring["METHOD"] == "circular" and ring["MODE"] == ""
        # make_scene puts the blob on the circular method (its own setting); planets are exact
        assert h["FINAL_CME_BLOB"].header["METHOD"] == "circular" and h["FINAL_PLANET"].header["METHOD"] == "exact"
        assert ring["ERRTYP"] == 0.0 and len(ring["WHY"]) <= 68
        parhash = p0["PARHASH"]
    with fits.open(os.path.join(folder, "scene.fits")) as h:
        assert h[0].header["PARHASH"] == parhash and h[0].header["SCNSIZE"] == n
        assert [hdu.name for hdu in h[1:]] == ["RING", "CME_BLOB", "POINTS"]
        assert h["POINTS"].data["FLUX"][0] == pytest.approx(PLANET_JY, rel=1e-12)   # Jy per source


def test_save_never_overwrites_unless_asked(cache, in_tmp):
    """A second save to the same folder goes to <folder>_2, then _3; overwrite=True writes into it."""
    res = observe(make_scene(16, cache, ring=False, blob=False), ["vortex:2"], verbose=False)
    base = str(in_tmp/"run")
    assert [res.save(base) for _ in range(3)] == [base, base+"_2", base+"_3"]
    assert res.save(base, overwrite=True) == base


def test_save_without_a_folder_uses_the_run_settings(cache, in_tmp):
    """save() with no folder: [run] output_dir / "scene" (a scene built in Python has no file name)."""
    res = observe(make_scene(16, cache, blob=False, planet=False), ["vortex:2"], verbose=False)
    assert res.save() == os.path.join("results", "scene")
    assert (in_tmp/"results"/"scene"/"vortex_c2.fits").exists()


def test_save_without_a_folder_is_named_after_the_parameter_file(cache, in_tmp):
    """[run] name = "auto" and a scene read from my_run.toml: save() writes to [run] output_dir / "my_run"."""
    res = observe(make_scene(16, cache, blob=False, planet=False), ["vortex:2"], verbose=False)
    first = res.save(str(in_tmp/"first"))
    os.replace(os.path.join(first, "params_used.toml"), str(in_tmp/"my_run.toml"))
    again = observe(Scene.from_file(str(in_tmp/"my_run.toml")), verbose=False)
    assert again.save(overwrite=True) == os.path.join("results", "my_run")
    assert (in_tmp/"results"/"my_run"/"vortex_c2.fits").exists()


def test_spatial_component_is_saved_with_its_mode_and_error(cache, in_tmp):
    """A component on the spatial method: METHOD "spatial", its MODE and estimated error in the FITS header, and back
    in Result.load()'s methods."""
    res = observe(make_scene(SPATIAL_N, cache, ring=False, planet=False, blob_method="auto"), ["vortex:2"],
                  verbose=False)
    method, mode, _, error = res.methods["vortex_c2", "cme_blob"]
    folder = res.save(str(in_tmp/"run"))
    with fits.open(os.path.join(folder, "vortex_c2.fits")) as h:
        header = h["FINAL_CME_BLOB"].header
        assert header["METHOD"] == method == "spatial" and header["MODE"] == mode
        assert header["ERRTYP"] == pytest.approx(error, rel=1e-12) and error > 0
    loaded = Result.load(folder)
    assert loaded.methods["vortex_c2", "cme_blob"][:2] == ("spatial", mode)
    assert loaded.methods["vortex_c2", "cme_blob"][3] == pytest.approx(error, rel=1e-12)


@pytest.mark.parametrize("n", SIZES)
def test_iwa_ignore_removes_light_from_the_images_only(n, cache):
    """iwa_ignore=True leaves out the target inside the inner working angle: less light in the ring's final image
    (circular method), the planet's exact PSF unchanged."""
    scene = make_scene(n, cache, blob=False)
    scene.add(Ring(radius=0.05*n*PIXSCALE, width=0.05*n*PIXSCALE, peak=0.05))   # light near the star
    full = observe(scene, ["vortex:2"], verbose=False)
    cut = observe(scene, ["vortex:2"], iwa_ignore=True, verbose=False)
    assert np.array_equal(cut.final(component="planet"), full.final(component="planet"))
    assert cut.final(component="ring_2").sum() < full.final(component="ring_2").sum()


@pytest.mark.parametrize("n", SIZES)
def test_load_gives_back_the_saved_result(n, cache, in_tmp):
    """Result.load(): the same grid, unit, star flux, final images, methods and planet positions (two coronagraphs)."""
    # the blob on the circular method (its own setting): no spatial model needed for charge 6
    res = observe(make_scene(n, cache, blob_method="circular"), ["vortex:2", "vortex:6"], verbose=False)
    loaded = Result.load(res.save(str(in_tmp/"run")))
    assert loaded.grid == res.grid and loaded.unit == res.unit and loaded.combine == "separate"
    assert loaded.star_flux == pytest.approx(res.star_flux, rel=1e-12) and loaded.star.distance == DISTANCE
    assert list(loaded.finals) == ["vortex_c2", "vortex_c6"] and loaded.coronagraphs == {}
    for c in res.finals:
        assert list(loaded.finals[c]) == list(res.finals[c])
        for name in res.finals[c]:
            assert np.allclose(loaded.final(c, name), res.final(c, name), rtol=1e-12, atol=0)
            assert loaded.methods[c, name][:2] == res.methods[c, name][:2]
    assert np.allclose(loaded.points["planet"][0], res.points["planet"][0], rtol=1e-12)
    assert np.allclose(loaded.points["planet"][1], res.points["planet"][1], rtol=1e-12)
    assert set(loaded.inputs) == {"ring", "cme_blob"}


@pytest.mark.parametrize("n", SIZES)
def test_points_only_scene_round_trip_uses_scnsize(n, cache, in_tmp):
    """A scene of planets only: scene.fits has no image, so load() takes the grid size from SCNSIZE."""
    scene = make_scene(n, cache, ring=False, blob=False, planet_pa=90)
    res = observe(scene, ["vortex:2"], verbose=False)
    assert res.inputs == {} and list(res.finals["vortex_c2"]) == ["planet"]
    folder = res.save(str(in_tmp/"run"))
    with fits.open(os.path.join(folder, "scene.fits")) as h:
        assert [hdu.name for hdu in h[1:]] == ["POINTS"] and h[0].header["SCNSIZE"] == n
    loaded = Result.load(folder)
    assert loaded.grid == res.grid and loaded.inputs == {}
    assert np.allclose(loaded.final(), res.final(), rtol=1e-12, atol=0)
    assert np.isclose(loaded.contrast("planet")["planet"], res.contrast("planet")["planet"], rtol=1e-12)
    loaded.plot(path=str(in_tmp/"loaded.png"), show=False)   # the planet's input panel is blank, marked


@pytest.mark.parametrize("n", SIZES)
def test_contrast_scene_round_trip_keeps_the_fraction(n, cache, in_tmp):
    """A contrast scene's planet: the saved table holds its fraction of the star; load() gives back the same contrast."""
    scene = make_scene(n, cache, star=CONTRAST_STAR, ring=False, planet_brightness=1e-6, planet_unit="contrast")
    res = observe(scene, ["vortex:2"], verbose=False)
    folder = res.save(str(in_tmp/"run"))
    with fits.open(os.path.join(folder, "scene.fits")) as h:
        assert h["POINTS"].data["FLUX"][0] == pytest.approx(1e-6, rel=1e-12) and h[0].header["STARFLUX"] == 0
    loaded = Result.load(folder)
    assert loaded.unit == "contrast_pixel" and loaded.star_flux is None
    assert np.isclose(loaded.contrast("planet")["planet"], res.contrast("planet")["planet"], rtol=1e-12)
    assert np.allclose(loaded.points["planet"][1], res.points["planet"][1], rtol=1e-12)


@pytest.mark.parametrize("n", SIZES)
def test_params_used_reproduces_the_run(n, cache, in_tmp):
    """The saved params_used.toml, read back with Scene.from_file, gives the same final images."""
    res = observe(make_scene(n, cache), ["vortex:2"], verbose=False)
    folder = res.save(str(in_tmp/"run"))
    again = observe(Scene.from_file(os.path.join(folder, "params_used.toml")), verbose=False)
    assert list(again.finals) == ["vortex_c2"]
    assert np.allclose(again.final(), res.final(), rtol=1e-12, atol=0)


def test_load_without_coronagraph_files_is_none(cache, in_tmp):
    """Result.load() of a folder without any coronagraph FITS file returns None."""
    res = observe(make_scene(16, cache, blob=False, planet=False), ["vortex:2"], verbose=False)
    folder = res.save(str(in_tmp/"run"))
    os.remove(os.path.join(folder, "vortex_c2.fits"))
    assert Result.load(folder) is None


# ---------------------------------------------------------------- strict grid, Target

def write_ring_fits(path, n, cache):
    """Write the n px ring of make_scene to path as a Jy/arcsec^2 FITS image (BUNIT, PIXSCALE, WAVELEN)."""
    scene = make_scene(n, cache, blob=False, planet=False)
    scene.render(verbose=False)
    header = fits.Header([("BUNIT", "Jy/arcsec^2"), ("PIXSCALE", PIXSCALE), ("WAVELEN", WAVELEN)])
    fits.PrimaryHDU(scene.images["ring"], header).writeto(str(path), overwrite=True)


@pytest.mark.parametrize("n", SIZES)
def test_strict_grid_runs_a_matching_image(n, cache, in_tmp):
    """[grid] resample = false: a FITS image on the grid's pixel scale, a module image and a planet all run."""
    write_ring_fits(in_tmp/"ring.fits", n, cache)
    scene = Scene(wavelength=WAVELEN, diameter=DIAMETER, resample=False, star=STAR,
                  compute={"psf_cache_dir": cache, "num_cores": 2})
    scene.add(FitsImage(file=str(in_tmp/"ring.fits")))
    scene.add(CmeBlob(r=0.2*n*PIXSCALE*DISTANCE, pa=40, width=0.06*n*PIXSCALE*DISTANCE, mass=1e21), method="circular")
    scene.add(Planet(mode="polar", pos=[0.35*n*PIXSCALE, 35.0], brightness=PLANET_JY))
    res = observe(scene, ["vortex:2"], verbose=False)
    assert set(res.finals["vortex_c2"]) == {"fits_image", "cme_blob", "planet"} and res.grid.size == n


def test_strict_grid_mismatch_stops_before_any_psf(cache, in_tmp):
    """[grid] resample = false and an image on another pixel scale: a ParamError naming the file, no PSF computed."""
    write_ring_fits(in_tmp/"ring.fits", 16, cache)
    empty = in_tmp/"empty_cache"
    scene = Scene(wavelength=WAVELEN, diameter=DIAMETER, resample=False, pixel_scale=1.3*PIXSCALE, star=STAR,
                  compute={"psf_cache_dir": str(empty), "num_cores": 2})
    scene.add(FitsImage(file=str(in_tmp/"ring.fits")))
    with pytest.raises(P.ParamError, match="ring.fits"):
        observe(scene, ["vortex:2"], verbose=False)
    assert not empty.exists() or not os.listdir(str(empty))


@pytest.mark.filterwarnings("ignore:FigureCanvasAgg is non-interactive")   # Target.plot_final calls plt.show()
@pytest.mark.parametrize("n", SIZES)
def test_target_and_observe_agree(n, cache, in_tmp):
    """Target (charge 2, circular, no planets) and observe() read the same parameter file to the same final image."""
    import matplotlib.pyplot as plt
    import toycoronagraph.main as toy
    write_ring_fits(in_tmp/"ring.fits", n, cache)
    (in_tmp/"ring.toml").write_text('[observation]\ndiameter = %r\n[[source]]\nkind = "fits_image"\nfile = "ring"\n'
                                    '[compute]\npsf_cache_dir = "%s"\nnum_cores = 2\n' % (DIAMETER, cache))
    target = toy.Target("ring.toml")
    image = target.plot_final(charge=2, add_planet=False, method="circular", plot_dpi=20)
    plt.close("all")
    res = observe(Scene.from_file("ring.toml"), ["vortex:2"], method="circular", verbose=False)
    assert image.shape == (n, n) and np.allclose(image, res.final("vortex_c2"), rtol=1e-12, atol=0)


def test_together_run_round_trip(cache, in_tmp):
    """A combine = "together" run saved and loaded keeps its "total" final image (and the planet's)."""
    res = observe(make_scene(16, cache), ["vortex:2"], combine="together", verbose=False)
    loaded = Result.load(res.save(str(in_tmp/"run")))
    assert list(loaded.finals["vortex_c2"]) == ["total", "planet"]
    assert np.allclose(loaded.final(), res.final(), rtol=1e-12, atol=0)


def test_together_run_round_trip_with_a_component_named_images(cache, in_tmp):
    """A point source named "images" in a together run does not collide with the summed images' extension."""
    scene = make_scene(16, cache, planet=False)
    scene.add(Planet(mode="polar", pos=[0.35*16*PIXSCALE, 35.0], brightness=PLANET_JY), name="images")
    res = observe(scene, ["vortex:2"], combine="together", verbose=False)
    loaded = Result.load(res.save(str(in_tmp/"run")))
    assert sorted(loaded.finals["vortex_c2"]) == ["images", "total"]
    for name in ("images", "total"):
        assert np.allclose(loaded.finals["vortex_c2"][name], res.finals["vortex_c2"][name], rtol=1e-12, atol=0)
