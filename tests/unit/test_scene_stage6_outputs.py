"""
Unit tests for stage F6 (tests/scene_stage6_outputs.py): the outputs of a scene run.

Covers the contrast_peak unit (Grid.peak_norm, the backends' unocculted(): vortex and Lazuli,
Coronagraph.peak_norm, Result.peak_norms / in_unit / final, Scene.in_unit, PEAKNORM in save / load, a
backend without unocculted()), the circular method's rings beyond the radial PSFs (out to the
corners), the inner working angle (Coronagraph.iwa, Result.iwa, its line in plot_profiles, IWA in the files),
the plot settings ([display] defaults, values given to a plot, keywords passed to matplotlib),
the Result figures (plot selections and the one-row layout, plot_profiles panels and
log range, plot_contributions: rows, shares, the dominant map), psf.cir_core_throughput_plot in
lambda/D, the lossless video writer (toycoronagraph.video) and observe.video.
Tiny scenes: 16 and 17 px at 4 px per lambda/D, with the template modules of tests/extra_sources.py.
"""
import json
import shutil
import subprocess
import numpy as np
import pytest
import matplotlib.pyplot as plt
from astropy.io import fits
from extra_sources import Ring, CmeBlob
from toycoronagraph import params as P
from toycoronagraph import units as U
from toycoronagraph import psf
from toycoronagraph.optics import Optics, register_backend
from toycoronagraph.scene import Scene, Star, Grid
from toycoronagraph.sources.planets import Planet
from toycoronagraph.observe import observe, video, Result, Coronagraph
from toycoronagraph.video import write_video

SIZES = [16, 17]
sizes = pytest.mark.parametrize("n", SIZES)
WAVELEN, DIAMETER, SAMPLING = 0.8e-6, 2.4, 4
PIXSCALE = WAVELEN/DIAMETER*P.RAD_TO_ARCSEC/SAMPLING     # arcsec: 4 px per lambda/D
STAR = Star(distance=3.22, radius=0.73, temperature=5080)   # Jy scenes
CONTRAST_STAR = Star(distance=3.22, radius=0.73)            # contrast scenes
LYOT = 0.95
needs_ffmpeg = pytest.mark.skipif(shutil.which("ffmpeg") is None or shutil.which("ffprobe") is None,
                                  reason="ffmpeg / ffprobe not found")


@pytest.fixture(scope="module")
def cache(tmp_path_factory):
    """One PSF cache folder for this module, so each PSF is computed once."""
    return str(tmp_path_factory.mktemp("psf_cache"))


@pytest.fixture(autouse=True)
def _work_in_tmp(in_tmp):
    """Run every test in its own temporary folder."""
    return in_tmp


def make_scene(n, cache, star=STAR, planet=True, **settings):
    """A ring (symmetric, Jy), a blob (contrast; circular method, so no spatial model is built) and a planet."""
    scene = Scene(wavelength=WAVELEN, diameter=DIAMETER, size=n, pixel_scale=PIXSCALE, star=star,
                  compute={"psf_cache_dir": cache, "num_cores": 2}, **settings)
    scene.add(Ring(radius=0.25*n*PIXSCALE, width=0.05*n*PIXSCALE, peak=0.02))
    scene.add(CmeBlob(r=0.2*n*PIXSCALE*3.22, pa=40, width=0.06*n*PIXSCALE*3.22, mass=1e21), method="circular")
    if planet:
        scene.add(Planet(mode="polar", pos=[0.35*n*PIXSCALE, 120], brightness=1e-6, brightness_unit="contrast"),
                  name="planet")
    return scene


def frames_in(path):
    """(number of frames, codec) of a video file (ffprobe)."""
    out = subprocess.run(["ffprobe", "-v", "error", "-count_frames", "-select_streams", "v:0", "-show_entries",
                          "stream=nb_read_frames,codec_name", "-of", "json", str(path)], capture_output=True, text=True,
                         check=True)
    stream = json.loads(out.stdout)["streams"][0]
    return int(stream["nb_read_frames"]), stream["codec_name"]


# ---------------------------------------------------------------- the unocculted peak

def test_grid_peak_norm_is_the_ideal_aperture_peak():
    """Grid.peak_norm = pi/4 (pixel / lambda/D)^2: an Airy peak per pixel, as a fraction of the light."""
    g = Grid(17, PIXSCALE, WAVELEN, DIAMETER)
    assert g.peak_norm == pytest.approx(np.pi/4/SAMPLING**2, rel=1e-12)


@sizes
def test_vortex_unocculted_peak_is_lyot_to_the_fourth(n, lyot=LYOT):
    """The star without the vortex: peak lyot^4 x the ideal aperture's (to 1%), throughput 1; Coronagraph.peak_norm
    is its value at the star, whatever the parity (the image's largest pixel only for an odd size)."""
    optics = Optics.vortex(2, n, n/SAMPLING/2, lyot)
    star = optics.backend().unocculted()
    ideal = np.pi/4*optics.pixel_scale**2
    peak = Coronagraph(optics, "v").peak_norm
    assert peak == star.block((0.0, 0.0), (1, 1))[0, 0]
    assert peak/ideal == pytest.approx(lyot**4, rel=0.01)
    assert star.throughput == 1.0
    if n % 2:  # the star on a pixel: its peak is the image's largest pixel
        assert star.frame().max() == pytest.approx(peak, rel=1e-9)
    else:      # between four pixels: lower
        assert star.frame().max() < 0.95*peak


def test_vortex_peak_norm_does_not_depend_on_the_parity():
    """The peak at the star is the same for 16 and 17 px (to 0.5%: only the pupil sampling differs)."""
    a, b = (Coronagraph(Optics.vortex(6, n, n/SAMPLING/2), "v").peak_norm for n in SIZES)
    assert a == pytest.approx(b, rel=0.005)


@pytest.mark.parametrize("charge, expected", [(0, 1.0), (2, LYOT**4)])
def test_lazuli_unocculted_peak(charge, expected):
    """Lazuli: the telescope alone peaks at the ideal aperture's value, with a charge (without the vortex, the Lyot
    stop kept) at lyot^4 of it (to 2%: 64 pupil samples); the planet far out of the telescope alone (charge = 0) reaches that peak."""
    pytest.importorskip("prysm")
    pytest.importorskip("config_um")
    optics = Optics.lazuli(charge, 17, 17/SAMPLING/2, [632e-9], [1.0], 632e-9, 3.0, 64, LYOT)
    backend = optics.backend()
    peak = Coronagraph(optics, "l").peak_norm
    assert peak/(np.pi/4*optics.pixel_scale**2) == pytest.approx(expected, rel=0.02)
    assert backend.unocculted().throughput == pytest.approx(1.0, rel=1e-12)
    if charge == 0:
        assert backend.source(1.5, 0.0).block((1.5, 0.0), (1, 1))[0, 0] == pytest.approx(peak, rel=1e-9)


@register_backend("no_unocculted_test")
class NoUnoccultedBackend:
    """A backend without unocculted(): its coronagraph has no contrast_peak reference."""
    rotationally_symmetric, bandlimit = True, 1.0

    def __init__(self, optics):
        self.optics = optics


def test_backend_without_unocculted_has_no_peak_norm():
    """peak_norm is None, and contrast_peak of its images stops with a clear error."""
    coro = Coronagraph(Optics("no_unocculted_test", 16, 2.0), "plain")
    assert coro.peak_norm is None
    res = Result(Grid(16, PIXSCALE, WAVELEN, DIAMETER), "contrast_pixel", CONTRAST_STAR, None, "separate", {})
    res.finals["plain"], res.peak_norms["plain"] = {"x": np.ones((16, 16))}, None
    with pytest.raises(U.UnitError, match="contrast_peak needs the peak"):
        res.final("plain", unit="contrast_peak")


# ---------------------------------------------------------------- contrast_peak in results

@sizes
def test_contrast_peak_uses_each_coronagraphs_peak(n, cache):
    """final(unit="contrast_peak") = contrast_pixel / that coronagraph's peak (two Lyot stops, two peaks); inputs use
    the ideal aperture's, in Result and in Scene alike."""
    scene = make_scene(n, cache)
    res = observe(scene, ["vortex:2", {"type": "vortex", "charge": 2, "lyot_mask_size": 0.8, "name": "small_lyot"}],
                  verbose=False)
    peaks = res.peak_norms
    assert peaks["small_lyot"]/peaks["vortex_c2"] == pytest.approx((0.8/LYOT)**4, rel=0.03)
    for name in res.finals:
        assert np.allclose(res.final(name, unit="contrast_peak"), res.final(name, unit="contrast_pixel")/peaks[name],
                           rtol=1e-12, atol=0)
        assert np.allclose(res.final(name, "ring", unit="contrast_peak"),
                           res.in_unit(res.finals[name]["ring"], "contrast_peak", coronagraph=name), rtol=1e-12, atol=0)
    ring = res.inputs["ring"]
    expected = res.in_unit(ring, "contrast_pixel")/scene.grid.peak_norm
    assert np.allclose(res.in_unit(ring, "contrast_peak"), expected, rtol=1e-12, atol=0)
    assert np.allclose(scene.in_unit(ring, "contrast_peak"), expected, rtol=1e-12, atol=0)


@sizes
def test_peaknorm_saved_and_loaded(n, cache):
    """PEAKNORM in each coronagraph's primary header; load() restores peak_norms and the contrast_peak images; an
    unknown peak is written 0 and loaded as None."""
    scene = make_scene(n, cache)
    res = observe(scene, ["vortex:2", "vortex:4"], verbose=False)
    folder = res.save("run", overwrite=True)
    for name, peak in res.peak_norms.items():
        assert fits.getheader("%s/%s.fits" % (folder, name))["PEAKNORM"] == pytest.approx(peak, rel=1e-14)
    again = Result.load(folder)
    assert again.peak_norms == pytest.approx(res.peak_norms, rel=1e-14)
    assert np.allclose(again.final("vortex_c4", unit="contrast_peak"), res.final("vortex_c4", unit="contrast_peak"),
                       rtol=1e-6, atol=0)
    res.peak_norms["vortex_c4"] = None
    folder = res.save("run_unknown", overwrite=True)
    assert fits.getheader("%s/vortex_c4.fits" % folder)["PEAKNORM"] == 0
    assert Result.load(folder).peak_norms["vortex_c4"] is None


# ---------------------------------------------------------------- the inner working angle

@pytest.mark.parametrize("charge, expected", [(2, 0.863), (4, 1.605), (6, 2.379)])
def test_vortex_iwa(charge, expected):
    """Coronagraph.iwa: where a source's throughput reaches 0.5 (lambda/D); the vortex's known values, the same
    for any image (a small field is searched beyond its edge)."""
    for n, half_width in [(17, 2.0), (65, 8.0), (97, 6.0)]:
        coro = Coronagraph(Optics.vortex(charge, n, half_width, LYOT), "v")
        assert coro.iwa == pytest.approx(expected, abs=0.02)
        throughput = lambda x: coro.optics.backend().source(x, 0.0).throughput
        assert throughput(coro.iwa-0.01) < 0.5 < throughput(coro.iwa+0.01)


def test_no_iwa_without_a_coronagraph():
    """A system without a mask (the vortex optics with charge = 0: the telescope and its Lyot stop) has no IWA."""
    assert Coronagraph(Optics.vortex(0, 33, 4.0, LYOT), "telescope").iwa is None


@sizes
def test_iwa_in_results_plots_and_files(n, cache):
    """Result.iwa: Coronagraph.iwa in arcsec, computed when asked; plot_profiles marks it (a dotted line in the curve's
    colour on one panel, grey on each coronagraph's panel; none with iwa=False); IWA in the files, loaded back."""
    res = observe(make_scene(n, cache, planet=False), ["vortex:2", "vortex:4"], verbose=False)
    assert res.iwas == {}
    for name, coro in res.coronagraphs.items():
        assert res.iwa(name) == pytest.approx(coro.iwa*res.grid.lambda_d, rel=1e-12)
    assert res.iwa("vortex_c2") < res.iwa("vortex_c4")
    dotted = lambda ax: {l.get_xdata()[0]: l.get_color() for l in ax.lines if l.get_linestyle() == ":"}
    fig = res.plot_profiles(axis="arcsec", show=False)
    colours = {l.get_label(): l.get_color() for l in fig.axes[0].lines}
    assert dotted(fig.axes[0]) == pytest.approx({res.iwa(c): colours[c] for c in ("vortex_c2", "vortex_c4")})
    fig = res.plot_profiles(components=True, axis="lambda/D", show=False)
    for ax, name in zip(fig.axes, ["vortex_c2", "vortex_c4"]):
        assert list(dotted(ax)) == pytest.approx([res.iwa(name)/res.grid.lambda_d])
    fig = res.plot_profiles(components=True, iwa=False, show=False)
    assert not any(dotted(ax) for ax in fig.axes)
    plt.close("all")
    folder = res.save("run", overwrite=True)
    assert fits.getheader("%s/vortex_c4.fits" % folder)["IWA"] == pytest.approx(res.iwa("vortex_c4"), rel=1e-12)
    assert Result.load(folder).iwas == pytest.approx(res.iwas, rel=1e-12)


def test_unknown_iwa_is_none():
    """A Result without the Coronagraph (e.g. loaded from files written before IWA) has no IWA, and plots without it."""
    res = Result(Grid(16, PIXSCALE, WAVELEN, DIAMETER), "contrast_pixel", CONTRAST_STAR, None, "separate", {})
    assert res.iwa("vortex_c2") is None


# ---------------------------------------------------------------- the circular method out to the corners

@sizes
def test_circular_rings_reach_the_corners(n, cache):
    """A uniform image through the circular method: no step at the inscribed circle (the last radial PSF's source):
    the rings beyond use the last PSF shifted, so the corners keep their light (without them they lose about half)."""
    coro = Coronagraph.get(Optics.vortex(2, n, n/SAMPLING/2), "v", num_cores=2, cache_dir=cache)
    final = psf.cir_psf(np.ones((n, n)), [], [], PIXSCALE, False, False, n, n, coro.radial_psfs_file())
    c = (n-1)/2
    r = np.hypot(*(np.indices((n, n))-c))
    last = psf.radial_source_positions(n)[-1]
    inside = final[(r > last-1.5) & (r <= last)].mean()
    beyond = final[(r > last) & (r <= last+1.5)].mean()
    assert beyond/inside > 0.9
    assert np.all(final[r > last] > 0.5*inside)


@sizes
def test_circular_method_stays_linear(n, cache):
    """With the corner rings the method is still linear in the image."""
    coro = Coronagraph.get(Optics.vortex(2, n, n/SAMPLING/2), "v", num_cores=2, cache_dir=cache)
    r = np.hypot(*(np.indices((n, n))-(n-1)/2))
    a, b = np.exp(-r/3), (r+1)**-2.0
    run = lambda img: psf.cir_psf(img, [], [], PIXSCALE, False, False, n, n, coro.radial_psfs_file())
    assert np.allclose(run(2*a-3*b), 2*run(a)-3*run(b), rtol=1e-10, atol=1e-12*np.abs(run(a)).max())


# ---------------------------------------------------------------- figures

def titles(fig):
    return [ax.get_title() for ax in fig.axes if ax.get_title()]


@sizes
def test_plot_selections_and_one_row(n, cache):
    """plot(): coronagraphs and components pick rows and columns, input=False drops the input row, and a single
    column is laid out as one row; contrast_peak titles say what each row is relative to."""
    res = observe(make_scene(n, cache), ["vortex:2", "vortex:4"], verbose=False)
    fig = res.plot(unit="contrast_peak", show=False)
    assert len(titles(fig)) == 3*4 and "ideal aperture" in fig._suptitle.get_text()
    fig = res.plot(components=["total"], show=False)
    assert titles(fig) == ["input: total", "vortex_c2: total", "vortex_c4: total"]
    assert len({round(ax.get_position().y0, 6) for ax in fig.axes if ax.get_title()}) == 1
    fig = res.plot(coronagraphs="vortex_c4", components=["ring", "planet"], input=False, show=False)
    assert titles(fig) == ["vortex_c4: ring", "vortex_c4: planet"]
    plt.close("all")


@sizes
def test_plot_profiles_panels(n, cache):
    """plot_profiles: one panel without components; one per coronagraph with them (legend on the first only), each
    component's input dashed and final solid in the same colour; log axes span at most log_range decades."""
    res = observe(make_scene(n, cache), ["vortex:2", "vortex:4"], verbose=False)
    fig = res.plot_profiles(unit="contrast_peak", show=False)
    curves = [l.get_label() for l in fig.axes[0].lines if not l.get_label().startswith("_")]   # not the IWA lines
    assert len(fig.axes) == 1 and curves == ["input", "vortex_c2", "vortex_c4"]
    fig = res.plot_profiles(components=True, show=False)
    panels = [ax for ax in fig.axes if ax.get_visible()]
    assert [ax.get_title() for ax in panels] == ["vortex_c2", "vortex_c4"]
    assert panels[0].get_legend() is not None and panels[1].get_legend() is None
    lines = {l.get_label(): l for l in panels[0].lines}
    assert lines["ring (input)"].get_linestyle() == "--" and lines["ring"].get_linestyle() == "-"
    assert lines["ring (input)"].get_color() == lines["ring"].get_color()
    bottom, top = panels[0].get_ylim()
    assert top/bottom == pytest.approx(2e8, rel=1e-6)
    plt.close("all")


@sizes
def test_plot_contributions_maps(n, cache):
    """plot_contributions: rows input + coronagraphs, a column per component and the dominant map; each row's shares
    sum to 100%; the dominant map is the largest fraction; the point source's input panel is blank; one colour bar.
    A single name, and input=False, select rows."""
    res = observe(make_scene(n, cache), ["vortex:2", "vortex:4"], verbose=False)
    fig = res.plot_contributions(show=False)
    names = ["ring", "cme_blob", "planet"]
    t = titles(fig)
    assert [x.split(":")[0] for x in t] == [r for r in ("input", "vortex_c2", "vortex_c4") for _ in range(4)]
    assert len(fig.axes) == 3*4+1  # and the colour bar
    panels = [ax for ax in fig.axes if ax.get_title()]
    for i in range(3):
        shares = [float(x.rsplit(", ", 1)[1].rstrip("%")) for x in t[4*i:4*i+3]]
        assert sum(shares) == pytest.approx(100, abs=0.1)
    assert not panels[2].images  # the planet's input: blank, marked
    for i, name in enumerate(["vortex_c2", "vortex_c4"], start=1):
        fractions, total = res.contributions(name), res.final(name)
        dominant = np.asarray(panels[4*i+3].images[0].get_array())
        expected = np.argmax([fractions[c] for c in names], axis=0)
        assert np.array_equal(dominant[total > 0], expected[total > 0])
    assert titles(res.plot_contributions("vortex_c4", input=False, show=False))[0] == "vortex_c4: ring, %s" % t[8].rsplit(", ", 1)[1]
    plt.close("all")


def test_core_throughput_plot_in_lambda_d(in_tmp, cache):
    """cir_core_throughput_plot(pixel_scale=...): distances and the IWA in lambda/D; without it, in pixels."""
    coro = Coronagraph.get(Optics.vortex(2, 17, 17/SAMPLING/2), "v", num_cores=2, cache_dir=cache)
    psfs = np.load(coro.radial_psfs_file())
    psf.cir_core_throughput_plot(psfs, plot_dpi=50, pixel_scale=1/SAMPLING)
    ax = plt.gcf().axes[0]
    assert np.allclose(ax.lines[0].get_xdata(), np.arange(len(psfs))/SAMPLING)
    iwa = psf.cir_iwa(psfs)
    assert ax.get_legend().get_texts()[0].get_text() == "IWA = %.3g lambda/D" % (iwa/SAMPLING)
    assert "lambda/D" in ax.get_xlabel()
    psf.cir_core_throughput_plot(psfs, plot_dpi=50)
    assert plt.gcf().axes[0].get_legend().get_texts()[0].get_text() == "IWA = %d pixel" % iwa
    plt.close("all")


# ---------------------------------------------------------------- videos

@needs_ffmpeg
def test_write_video_is_lossless_vp9(in_tmp):
    """write_video: one VP9 frame per call of update, in a .webm file."""
    fig, ax = plt.subplots(figsize=(2, 2))
    im = ax.imshow(np.zeros((8, 8)), vmin=0, vmax=4)
    path = write_video(fig, lambda i: [im.set_data(np.full((8, 8), i))], 5, str(in_tmp/"v.webm"), fps=5, dpi=40)
    plt.close(fig)
    assert frames_in(path) == (5, "vp9")


@needs_ffmpeg
@sizes
def test_observe_video(n, cache):
    """observe.video: one frame per setting; the others' final images reused (only the changing ring's final images are
    computed); difference=True gives the ring alone, growing outwards; its settings come back; the default file name."""
    scene = make_scene(n, cache)
    coro = observe(scene, ["vortex:2"], verbose=False).coronagraphs["vortex_c2"]
    before = len(coro.computed)
    frames = [{"radius": f*n*PIXSCALE} for f in (0.25, 0.3, 0.35)]
    path, full = video(scene, "vortex:2", "ring", frames, unit="contrast_pixel", dpi=40)
    assert path == "ring_vortex_c2.webm" and frames_in(path) == (3, "vp9")
    assert [what for what, _ in coro.computed[before:]] == ["final image"]*2  # 0.25 n px: the scene as observed
    _, alone = video(scene, "vortex:2", "ring", frames, path="alone.webm", unit="contrast_pixel", dpi=40, difference=True)
    others = observe(scene, ["vortex:2"], components=["cme_blob", "planet"], verbose=False)
    others = others.final("vortex_c2", unit="contrast_pixel")
    for f, a in zip(full, alone):
        assert np.allclose(f-a, others, rtol=1e-12, atol=1e-12*np.abs(f).max())
    grid = scene.grid
    r = np.hypot(*grid.coordinates())
    mean_r = [(a*r).sum()/a.sum() for a in alone]
    assert mean_r[0] < mean_r[1] < mean_r[2]
    assert scene.component_of["ring"].source.settings["radius"] == pytest.approx(0.25*n*PIXSCALE)


# ---------------------------------------------------------------- notes

def test_unused_star_notes_wait_for_files_that_decide_their_unit(cache):
    """A scene_fits image's unit comes from its file, so the star's keys are not called unused before it is read
    (here its Jy images and a planet in contrast need the star's flux); a scene of contrast modules only still gets
    the note that the temperature is not used."""
    from toycoronagraph.sources.images import SceneFits
    jy = make_scene(16, cache, planet=False)
    jy.render(verbose=False)
    jy.write("jy_scene.fits")
    scene = Scene(star=STAR, compute={"psf_cache_dir": cache}, diameter=DIAMETER)
    scene.add(SceneFits(file="jy_scene.fits"))
    scene.add(Planet(mode="polar", pos=[0.3*16*PIXSCALE, 30], brightness=1e-6, brightness_unit="contrast"), name="planet")
    scene.render(verbose=False)
    assert scene.internal_unit == "Jy/arcsec^2" and not [n for n in scene.notes if "is not used" in n]
    contrast = Scene(wavelength=WAVELEN, diameter=DIAMETER, size=16, pixel_scale=PIXSCALE, star=STAR)
    contrast.add(Planet(mode="polar", pos=[0.3*16*PIXSCALE, 30], brightness=1e-6, brightness_unit="contrast"))
    contrast.render(verbose=False)
    assert any("temperature is not used" in n for n in contrast.notes)


# ---------------------------------------------------------------- comparing two sets of PSFs

@sizes
def test_plot_compare_two_coronagraphs(n, cache):
    """plot_compare: both final images on one colour scale; the relative difference (second - first) / first where
    the first is within log_range decades of its peak (NaN, grey, elsewhere), its median and 99th percentile in the
    title; the azimuthal means of the input and of both; one component or the total."""
    res = observe(make_scene(n, cache, planet=False), ["vortex:2", "vortex:4"], verbose=False)
    fig = res.plot_compare("vortex_c2", "vortex_c4", component="ring", unit="contrast_peak", log_range=4, show=False)
    first, second, diff, profiles = fig.axes[:4]
    assert first.get_title() == "vortex_c2: ring" and second.get_title() == "vortex_c4: ring"
    assert first.images[0].get_clim() == second.images[0].get_clim()
    a, b = (res.final(c, "ring", unit="contrast_peak") for c in ("vortex_c2", "vortex_c4"))
    kept = a > a.max()*1e-4
    shown = np.asarray(np.ma.filled(diff.images[0].get_array(), np.nan), dtype=float)
    assert np.allclose(shown[kept], ((b-a)/a)[kept]) and np.all(np.isnan(shown[~kept]))
    size = np.abs((b-a)/a)[kept]
    assert diff.get_title() == "relative difference: median %.2g%%, 99%% within %.2g%%" % (
        100*np.median(size), 100*np.percentile(size, 99))
    assert [l.get_label() for l in profiles.lines] == ["input", "vortex_c2", "vortex_c4"]
    fig = res.plot_compare("vortex_c2", "vortex_c4", show=False)  # the total
    assert fig.axes[0].get_title() == "vortex_c2: total"
    plt.close("all")


# ---------------------------------------------------------------- plot settings

def test_display_setting_order(cache):
    """A value given to the plot wins; else [display] of the parameters; else the declared default."""
    from toycoronagraph.scene import display_setting
    scene = make_scene(16, cache, display={"log_range": 3}, planet=False)
    scene.render(verbose=False)
    assert display_setting(scene.params, "log_range", 7) == 7
    assert display_setting(scene.params, "log_range", None) == 3
    assert display_setting(scene.params, "dpi", None) == 150
    assert display_setting(None, "log", None) is True and display_setting(None, "axis", None) == "arcsec"


def colour_span(ax):
    """Decades shown by an image panel's colour scale (log scale: vmax - vmin of log10 values)."""
    low, high = ax.images[0].get_clim()
    return high-low


def test_plots_follow_display_and_pass_keywords_to_matplotlib(cache):
    """[display] log_range sets the images' colour span (profiles: PROFILE_DECADES more); a value given to the plot
    wins; other keywords reach matplotlib: imshow for images and share maps, plot for every curve."""
    from toycoronagraph.scene import PROFILE_DECADES
    scene = make_scene(17, cache, display={"log_range": 3, "unit": "contrast_pixel"}, planet=False)
    res = observe(scene, ["vortex:2"], verbose=False)
    fig = scene.plot(show=False, cmap="magma")
    assert colour_span(fig.axes[0]) == pytest.approx(3) and fig.axes[0].images[0].get_cmap().name.startswith("magma")
    fig = res.plot(show=False, log_range=6, interpolation="nearest")
    assert colour_span(fig.axes[0]) == pytest.approx(6) and fig.axes[0].images[0].get_interpolation() == "nearest"
    fig = res.plot_profiles(components=True, show=False, iwa=False, lw=3)
    bottom, top = fig.axes[0].get_ylim()
    assert np.log10(top/bottom) == pytest.approx(3+PROFILE_DECADES+np.log10(2))   # the top is 2x the highest point
    assert {line.get_linewidth() for line in fig.axes[0].lines} == {3}
    fig = res.plot_contributions(show=False, cmap="magma")
    assert fig.axes[0].images[0].get_cmap().name.startswith("magma")
    fig = res.plot_compare("vortex_c2", "vortex_c2", show=False, cmap="magma")
    assert fig.axes[0].images[0].get_cmap().name.startswith("magma") and colour_span(fig.axes[0]) == pytest.approx(3)
    plt.close("all")


def test_saved_resolution_follows_display_dpi(cache):
    """The saved file's resolution: [display] dpi, or the dpi given to the plot."""
    scene = make_scene(17, cache, display={"dpi": 40}, planet=False)
    scene.plot(path="low.png", show=False)
    scene.plot(path="high.png", show=False, dpi=80)
    low, high = (plt.imread(name).shape[1] for name in ("low.png", "high.png"))
    assert high/low == pytest.approx(2, rel=0.05)


# ---------------------------------------------------------------- the sampling warning

def sampling_warnings(*args):
    """The SamplingWarning messages of optics_for(*args), shown whatever the warning filters."""
    import warnings
    from toycoronagraph.optics import SamplingWarning
    from toycoronagraph.observe import optics_for
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        optics_for(*args)
    return [str(w.message) for w in record if issubclass(w.category, SamplingWarning)]


@pytest.mark.parametrize("sampling, warned", [(2.36, True), (4.9, True), (5.0, False), (8.0, False)])
def test_vortex_warns_below_five_pixels_per_lambda_d(sampling, warned):
    """The vortex model samples its pupil from the image's pixels: below MIN_SAMPLING (5) pixels per lambda/D a
    SamplingWarning names the sampling and the largest pixel scale that is enough."""
    from toycoronagraph.optics import MIN_SAMPLING
    values = {p.name: p.default for p in P.CORONAGRAPHS["vortex"][1]}
    messages = sampling_warnings("vortex", values, 41, 41/sampling/2, WAVELEN, DIAMETER)
    assert len(messages) == warned
    if warned:
        assert "%.3g pixels per lambda/D" % sampling in messages[0] and "sampling = %d" % MIN_SAMPLING in messages[0]
        assert "pixel_scale <= %.4g arcsec" % (WAVELEN/DIAMETER*P.RAD_TO_ARCSEC/MIN_SAMPLING) in messages[0]


@pytest.mark.parametrize("samples, warned", [(256, True), (512, False)])
def test_lazuli_warns_when_its_pupil_is_too_coarse_for_the_image(samples, warned):
    """Lazuli's pupil has its own sampling: fewer than 8/3 x MIN_SAMPLING samples per lambda/D of the image's
    half-width (a 69 lambda/D wide image: 256 samples) gives a SamplingWarning naming the samples needed."""
    pytest.importorskip("prysm")
    pytest.importorskip("config_um")
    values = {p.name: p.default for p in P.CORONAGRAPHS["lazuli"][1]}
    values.update(charge=6, pupil_samples=samples)
    messages = sampling_warnings("lazuli", values, 345, 34.5, 632e-9, 3.0)
    assert len(messages) == warned
    if warned:
        assert "set pupil_samples >= 460" in messages[0]


@pytest.mark.parametrize("missing", ["prysm", "config_um", "config_um_wcc"])
def test_lazuli_without_its_packages_says_how_to_install_them(missing, monkeypatch):
    """type = "lazuli" without prysm or a UASAL config package: an ImportError naming the package and its pip command."""
    import sys
    import toycoronagraph
    from toycoronagraph import optics
    from toycoronagraph.observe import optics_for, LAZULI_PACKAGES
    # toycoronagraph.lazuli is imported afresh (earlier tests may have imported it), with the package missing; the
    # backend registry a fresh import would change is restored afterwards
    monkeypatch.delitem(sys.modules, "toycoronagraph.lazuli", raising=False)
    monkeypatch.delattr(toycoronagraph, "lazuli", raising=False)
    monkeypatch.setattr(optics, "_BACKENDS", dict(optics._BACKENDS))
    for name in [m for m in sys.modules if m == missing or m.startswith(missing + ".")] + [missing]:
        monkeypatch.setitem(sys.modules, name, None)                           # as if it were not installed
    values = {p.name: p.default for p in P.CORONAGRAPHS["lazuli"][1]}
    with pytest.raises(ImportError, match="needs %s, which is not installed: pip install %s"
                       % (missing, LAZULI_PACKAGES[missing].replace(".", r"\.").replace("+", r"\+"))):
        optics_for("lazuli", values, 17, 2.0, 632e-9, 3.0)
