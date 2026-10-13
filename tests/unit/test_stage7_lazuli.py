"""
Unit tests for stage 7: the Lazuli telescope backend ([[coronagraph]] type = "lazuli", computed with prysm).

Modules covered:
    toycoronagraph.lazuli    config parameters, detector_pixel_scale, spectral_weights, read_opd,
                             pupil_samples, optics_from_table, LazuliBackend / LazuliSource
    toycoronagraph.optics    Optics.lazuli, Optics.settings / from_settings with tuple parameters
    toycoronagraph.params    the "lazuli" coronagraph type (unit strings, "auto", cross-key checks, names);
                             params_docs.json_schema for its keys, params_docs.new_config lazuli templates
    toycoronagraph.observe   optics_for / OPTICS_BUILDERS, coronagraph_for, Coronagraph.tag, observe(),
                             the saved FITS header; psf._same_settings (PSF cache checks)
    toycoronagraph.main / toycoronagraph.config   Target and the older layout with type = "lazuli"

Small images (N = 15-22 pixels at 4 pixels per lambda/D, odd and even where the centre convention
matters; 2P+1 pixels for the one-period light sums) and 32-128 pupil samples keep every test short;
the brute-force images, the comparisons with hcipy at 512 pupil samples and the figures stay in
tests/stage7_lazuli.py.
"""
import os
import re
import sys
import subprocess
import textwrap
import numpy as np
import pytest
from astropy.io import fits

pytest.importorskip("prysm")
pytest.importorskip("config_um")
pytest.importorskip("config_um_wcc")

from extra_sources import Ring  # noqa: E402  (tests/ is on sys.path, see conftest.py)
from toycoronagraph import lazuli as L  # noqa: E402
from toycoronagraph import params as P  # noqa: E402
from toycoronagraph import config as C  # noqa: E402
from toycoronagraph import psf  # noqa: E402
from toycoronagraph.optics import Optics  # noqa: E402
from toycoronagraph.scene import Scene, Grid  # noqa: E402
from toycoronagraph.sources.images import FitsImage  # noqa: E402
from toycoronagraph.observe import (observe, optics_for, coronagraph_for, Coronagraph,  # noqa: E402
                                    OPTICS_BUILDERS)
from toycoronagraph.tool import arcsec_to_pixel  # noqa: E402

pytestmark = pytest.mark.filterwarnings("ignore:FigureCanvasAgg is non-interactive")

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WAVELEN = 632e-9     # m: the scene's wavelength
LAZULI_D = 3.0       # m: Lazuli's pupil (beam size on M1); as scene diameter, scene lambda/D = Lazuli's
SAMPLING = 4         # image pixels per scene lambda/D
PUPIL = 64           # pupil samples across D (the configuration's 512 is slower)
PIXSCALE = WAVELEN/LAZULI_D*P.RAD_TO_ARCSEC/SAMPLING  # arcsec per pixel at SAMPLING


# ---------------------------------------------------------------- helpers

def lazuli_optics(charge=0, n=21, wavelengths=(WAVELEN,), weights=(1.0,), diameter=LAZULI_D, pupil=PUPIL, opd="",
                  psf_range=None):
    """
    An Optics.lazuli on a small image (scene wavelength WAVELEN, Lyot stop 0.95).

    Args:
        charge (int): 0 (the telescope alone) or an even vortex charge.
        n (int): Image size in pixels.
        wavelengths, weights (tuple of float): PSF wavelengths in m and their weights.
        diameter (float): The scene's diameter in m (LAZULI_D: scene lambda/D = Lazuli's lambda/D).
        pupil (int): Pupil samples across D (the OPD map's size when opd is set).
        opd (str): OPD map file, "" for none.
        psf_range (float): Half-width of the image in scene lambda/D; None: SAMPLING pixels per lambda/D.

    Returns:
        Optics: kind 'lazuli'.
    """
    return Optics.lazuli(charge, n, psf_range or n/SAMPLING/2, wavelengths, weights, WAVELEN, diameter, pupil,
                         0.95, opd)


def table_values(table):
    """
    A [[coronagraph]] table checked by params.make, with the defaults of its type filled in.

    Args:
        table (dict): The table as written in a parameter file (with "type").

    Returns:
        dict: Every key of the type's declaration (as optics_for wants them).
    """
    checked = P.make({}, [], [table]).coronagraphs[0]
    return dict({p.name: p.default for p in P.CORONAGRAPHS[checked.type][1]}, **checked.values)


def problems(table):
    """
    The error message of params.make for one [[coronagraph]] table.

    Args:
        table (dict): The table.

    Returns:
        str: The ParamError's text (the test fails if there is no error).
    """
    with pytest.raises(P.ParamError) as err:
        P.make({}, [], [table])
    return str(err.value)


def astigmatism(n, rms_nm):
    """
    An OPD map rho^2 cos(2 theta) (rho = 1 at the pupil edge) of a given RMS.

    Args:
        n (int): Samples on one edge.
        rms_nm (float): RMS over the pupil, in nm.

    Returns:
        np.ndarray: (n, n) float64 OPD in nm (rows = y, columns = x).
    """
    line = (np.arange(n)-(n-1)/2)/(n/2)
    x, y = np.meshgrid(line, line)
    return rms_nm*np.sqrt(6)*(x**2-y**2)


def ring_image(n):
    """
    A symmetric Gaussian ring centred on pixel ((n-1)/2, (n-1)/2), in Jy/arcsec^2.

    Args:
        n (int): Image size.

    Returns:
        np.ndarray: (n, n) float64.
    """
    y, x = np.mgrid[:n, :n]-(n-1)/2
    return np.exp(-0.5*((np.hypot(x, y)-0.25*n)/1.0)**2)


def write_image(path, img, pixscale=PIXSCALE):
    """
    Write an image as a FITS file in Jy/arcsec^2 with PIXSCALE and WAVELEN (WAVELEN) in its header.

    Args:
        path (str): The file.
        img (np.ndarray): (n, n) image.
        pixscale (float): Pixel scale in arcsec.
    """
    fits.PrimaryHDU(img, fits.Header([("BUNIT", "Jy/arcsec^2"), ("PIXSCALE", pixscale), ("WAVELEN", WAVELEN)])
                    ).writeto(path, overwrite=True)


def legacy_toml(path, image="ring", coronagraph='type = "lazuli"\npupil_samples = 32\n'):
    """
    Write an older-layout parameter file ([target], [telescope], [coronagraph], [compute]) for Target.

    Args:
        path (str): The .toml file.
        image (str): fits_filename (relative to the working directory).
        coronagraph (str): The lines of the [coronagraph] section.
    """
    with open(path, "w") as f:
        f.write('[target]\nfits_filename = "%s"\n[telescope]\ndiameter = %g\n[coronagraph]\n%s'
                '[compute]\npsf_cache_dir = "cache"\nnum_cores = 1\n' % (image, LAZULI_D, coronagraph))


def ring_scene(n, **settings):
    """
    A Scene with a Ring on n x n pixels at SAMPLING pixels per lambda/D (D = LAZULI_D), PSF cache in ./cache.

    Args:
        n (int): Image size.
        settings: Other Scene settings.

    Returns:
        Scene: Not rendered yet (one core, so no worker processes).
    """
    scene = Scene(**dict(dict(wavelength=WAVELEN, diameter=LAZULI_D, size=n, pixel_scale=PIXSCALE,
                              compute={"psf_cache_dir": "cache", "num_cores": 1}), **settings))
    scene.add(Ring(radius=0.25*n*PIXSCALE, width=0.04*n*PIXSCALE, peak=1.0))
    return scene


@pytest.fixture
def fresh_coronagraphs():
    """Forget the Coronagraph objects shared between runs (Coronagraph._instances), before and after a test."""
    Coronagraph._instances.clear()
    yield
    Coronagraph._instances.clear()


# ---------------------------------------------------------------- a) parameters from the configs

def test_lazuli_parameters_match_the_configs():
    """D = 3.0 m (beam size on M1), f_eff = 42.047 m, 512 pupil samples and 3.76 um detector pixels."""
    t = L.lazuli_parameters()
    assert (t.diameter, t.focal_length, t.pupil_samples, t.pixel_size) == (3.0, 42.047, 512, 3.76e-6)
    assert isinstance(t.pupil_samples, int) and isinstance(t.diameter, float)


def test_lazuli_parameters_are_read_once():
    """lazuli_parameters is cached: the same object on every call."""
    assert L.lazuli_parameters() is L.lazuli_parameters()


def test_lazuli_curves_are_wavelength_tables_in_metres():
    """The six throughput curves (M1, M2-M4, r' filter, QE) are matching 1-D float arrays, wavelengths in m."""
    curves = L.lazuli_parameters().curves
    assert set(curves) == {"m1", "m2", "m3", "m4", "filter", "qe"}
    for name, (lam, values) in curves.items():
        assert lam.ndim == 1 and lam.shape == values.shape and lam.dtype == np.float64, name
        assert lam.min() < 600e-9 < lam.max() < 1e-4, name  # covers the r' band, in m (not nm or um)


def test_detector_pixel_scale_is_pixel_over_focal_length():
    """detector_pixel_scale = pixel_size / f_eff in arcsec, about 0.018445."""
    expected = 3.76e-6/42.047*P.RAD_TO_ARCSEC
    assert L.detector_pixel_scale() == pytest.approx(expected, rel=1e-12)
    assert L.detector_pixel_scale() == pytest.approx(0.018445, abs=1e-6)


# ---------------------------------------------------------------- the Lazuli team's get_psf

def test_pupil_is_a_centred_gray_edged_circle():
    """pupil(P): (P, P) coordinates in mm centred on the array and a symmetric 0-1 aperture of diameter D."""
    x, y, aperture = L.pupil(PUPIL)
    assert x.shape == y.shape == aperture.shape == (PUPIL, PUPIL)
    assert np.isclose(x[0, -1]-x[0, 0], 3000*(PUPIL-1)/PUPIL) and np.isclose(x.mean(), 0) and np.isclose(y.mean(), 0)
    assert aperture.min() >= 0 and aperture.max() == 1
    assert np.allclose(aperture, aperture[::-1]) and np.allclose(aperture, aperture.T)
    assert aperture[0, 0] == 0  # the corners are outside the circle


def test_tele_data_uses_the_backend_values():
    """tele_data: the backend's pupil, no OPD, wavelengths in um, dx = D/P mm, fnum = f_eff/D, the detector pixel."""
    td = L.tele_data([0.62, 0.65], samples=PUPIL)
    assert np.array_equal(td["Aperture"], L.pupil(PUPIL)[2]) and not td["OPD"].any()
    assert td["Wavelengths"] == [0.62, 0.65]
    assert td["dx"] == pytest.approx(3000/PUPIL) and td["fnum"] == pytest.approx(42.047/3.0)
    assert td["Pixsize"] == pytest.approx(3.76)
    assert np.allclose(td["Throughput"], L.spectral_weights([0.62e-6, 0.65e-6]))


@pytest.mark.parametrize("shift_px", [(0.0, 0.0), (2.3, -1.6)])
def test_get_psf_is_charge_0_on_the_detector_pixels(shift_px):
    """get_psf (wavelengths in um) = the charge-0 backend on the detector's pixels, once both sum to 1;
    get_psf's shift moves its grid, so the backend's source sits at minus the shift."""
    n = 21
    td = L.tele_data([WAVELEN*1e6], samples=PUPIL)
    pixel = L.detector_pixel_scale()/(WAVELEN/LAZULI_D*P.RAD_TO_ARCSEC)  # scene lambda/D
    b = lazuli_optics(n=n, psf_range=n*pixel/2).backend()
    shift = tuple(v*td["Pixsize"] for v in shift_px)
    theirs = L.get_psf(td, shift=shift, samples=n)
    ours = b.source(-shift[0]/b.um_per_lambda_d, -shift[1]/b.um_per_lambda_d).frame()
    assert theirs.shape == (n, n) and theirs.sum() == pytest.approx(1)
    assert np.abs(ours/ours.sum()-theirs).max() < 1e-10*theirs.max()


def test_get_psf_scalar_throughput_weighs_every_wavelength_alike():
    """A float Throughput with several wavelengths weighs them equally (a list of equal weights)."""
    td = L.tele_data([0.6, 0.66], samples=PUPIL)
    equal = L.get_psf(dict(td, Throughput=[1.0, 1.0]), samples=15)
    assert np.allclose(L.get_psf(dict(td, Throughput=1.0), samples=15), equal, rtol=1e-12, atol=0)


# ---------------------------------------------------------------- b) spectral weights

def test_spectral_weights_are_nonnegative_and_sum_to_one():
    """spectral_weights returns (n,) float64 weights >= 0 summing to 1."""
    lam = np.arange(400, 801, 5)*1e-9
    w = L.spectral_weights(lam)
    assert w.shape == lam.shape and w.dtype == np.float64
    assert np.all(w >= 0) and w.sum() == pytest.approx(1.0, abs=1e-12)


def test_spectral_weights_positive_across_the_r_band():
    """Inside 570-690 nm every weight is within a factor ~3 of the largest."""
    lam = np.arange(400, 801, 5)*1e-9
    w = L.spectral_weights(lam)
    band = (lam >= 570e-9) & (lam <= 690e-9)
    assert w[band].min()/w.max() > 0.3


def test_spectral_weights_vanish_outside_the_r_band():
    """Below 530 nm and above 740 nm the r' filter lets no light through."""
    lam = np.arange(400, 801, 5)*1e-9
    w = L.spectral_weights(lam)
    outside = (lam < 530e-9) | (lam > 740e-9)
    assert w[outside].max() < 1e-15


def test_spectral_weights_of_one_wavelength_is_one():
    """A single in-band wavelength gets weight 1."""
    assert L.spectral_weights([WAVELEN]) == pytest.approx([1.0], abs=1e-15)


def test_spectral_weights_out_of_band_wavelength_gets_zero():
    """One wavelength outside the band and one inside: weights (0, 1)."""
    assert L.spectral_weights([450e-9, WAVELEN]) == pytest.approx([0.0, 1.0], abs=1e-15)


@pytest.mark.parametrize("wavelengths", [[450e-9], [900e-9, 1e-6], []], ids=["blue", "infrared", "empty"])
def test_spectral_weights_without_light_raise(wavelengths):
    """No wavelength in the r' band (or none at all): ValueError naming the r' filter."""
    with pytest.raises(ValueError, match="r' filter"):
        L.spectral_weights(wavelengths)


def test_spectral_weights_follow_the_band_not_its_sampling():
    """Regression: the weights include the band each wavelength stands for, so the weighted mean wavelength of
    550-700 nm is the same (to 2 nm) sampled every 10 nm or 4x denser below 600 nm, and in any order."""
    even = np.arange(550, 701, 10)*1e-9
    uneven = np.concatenate([np.arange(550, 600, 2.5), np.arange(600, 701, 25)])*1e-9
    means = [np.sum(L.spectral_weights(lams)*lams) for lams in (even, uneven)]
    assert abs(means[1]-means[0]) < 2e-9
    assert L.spectral_weights(uneven[::-1]) == pytest.approx(L.spectral_weights(uneven)[::-1], abs=1e-15)


# ---------------------------------------------------------------- c) the "lazuli" coronagraph type in params

def test_lazuli_declaration_defaults():
    """Declared defaults: charge 0, Lyot stop 0.95, wavelengths / throughput / pupil_samples "auto", no OPD."""
    defaults = {p.name: p.default for p in P.CORONAGRAPHS["lazuli"][1]}
    assert defaults == {"charge": 0, "lyot_mask_size": 0.95, "wavelengths": "auto", "throughput": "auto", "opd": "",
                        "pupil_samples": "auto", "symmetric": "auto"}


def test_wavelength_unit_strings_are_converted_to_metres():
    """wavelengths may mix numbers in m and text with a length unit ("620 nm", "0.645 um")."""
    values = P.make({}, [], [{"type": "lazuli", "wavelengths": ["620 nm", "0.645 um", 6.5e-7]}]).coronagraphs[0].values
    assert np.allclose(values["wavelengths"], [620e-9, 645e-9, 6.5e-7], rtol=1e-12, atol=0)


def test_auto_wavelengths_and_throughput_are_kept_as_auto():
    """"auto" passes the checks unchanged (resolved later by optics_from_table)."""
    values = P.make({}, [], [{"type": "lazuli", "wavelengths": "auto", "throughput": "auto"}]).coronagraphs[0].values
    assert (values["wavelengths"], values["throughput"]) == ("auto", "auto")


def test_auto_wavelengths_use_the_scene_wavelength():
    """wavelengths = "auto": one wavelength, the scene's, with weight 1."""
    optics = L.optics_from_table(table_values({"type": "lazuli", "pupil_samples": PUPIL}), 21, 2.6, 640e-9, 2.4, "%d %d")
    assert optics.param("wavelengths") == (640e-9,) and optics.param("weights") == (1.0,)
    assert (optics.param("scene_wavelength"), optics.param("scene_diameter")) == (640e-9, 2.4)


def test_auto_throughput_is_the_spectral_weights():
    """throughput = "auto": the weights are spectral_weights(wavelengths)."""
    values = table_values({"type": "lazuli", "wavelengths": ["620 nm", "650 nm"], "pupil_samples": PUPIL})
    optics = L.optics_from_table(values, 21, 2.6, WAVELEN, LAZULI_D, "%d %d")
    assert np.allclose(optics.param("weights"), L.spectral_weights([620e-9, 650e-9]), rtol=1e-12, atol=0)


def test_auto_throughput_outside_the_band_raises():
    """A scene at 1 um with "auto" wavelengths and throughput: ValueError naming the r' filter."""
    with pytest.raises(ValueError, match="r' filter"):
        L.optics_from_table(table_values({"type": "lazuli", "pupil_samples": PUPIL}), 21, 2.6, 1e-6, LAZULI_D, "%d %d")


def test_explicit_throughput_is_normalised():
    """throughput [1, 3] gives weights (0.25, 0.75)."""
    values = table_values({"type": "lazuli", "wavelengths": [6e-7, 6.5e-7], "throughput": [1, 3], "pupil_samples": PUPIL})
    optics = L.optics_from_table(values, 21, 2.6, WAVELEN, LAZULI_D, "%d %d")
    assert optics.param("weights") == pytest.approx((0.25, 0.75), abs=1e-15)


@pytest.mark.parametrize("charge", [1, 3, 5])
def test_odd_charge_is_refused(charge):
    """An odd charge is a ParamError: "a vortex charge must be even"."""
    assert "charge = %d: a vortex charge must be even" % charge in problems({"type": "lazuli", "charge": charge})


def test_negative_charge_is_refused():
    """charge must be >= 0 (0 = the telescope alone)."""
    assert "charge = -2: must be >= 0" in problems({"type": "lazuli", "charge": -2})


@pytest.mark.parametrize("table, words", [
    ({"type": "lazuli", "wavelengths": [6e-7, 6.5e-7], "throughput": [1.0]}, "throughput has 1 values for 2 wavelength(s)"),
    ({"type": "lazuli", "throughput": [1.0, 2.0]}, "throughput has 2 values for 1 wavelength(s)"),
], ids=["two_wavelengths_one_weight", "auto_wavelength_two_weights"])
def test_throughput_length_mismatch_is_refused(table, words):
    """throughput needs one weight per wavelength ("auto" wavelengths count as one)."""
    assert words in problems(table)


def test_bad_wavelength_item_adds_no_wrong_throughput_count():
    """Two wavelengths (one with a bad unit) and two weights: only the unit is reported, not a wrong count."""
    message = problems({"type": "lazuli", "wavelengths": ["632 kg", "650 nm"], "throughput": [1.0, 1.0]})
    assert 'item "632 kg": cannot convert to m' in message
    assert "for 1 wavelength(s)" not in message


def test_all_zero_throughput_is_refused():
    """Weights that are all 0 are refused."""
    assert "throughput: the weights are all 0" in problems({"type": "lazuli", "throughput": [0.0]})


def test_empty_wavelength_list_is_refused():
    """wavelengths = [] is refused (at least one wavelength, or "auto")."""
    assert "wavelengths = []: give at least one wavelength" in problems({"type": "lazuli", "wavelengths": []})


@pytest.mark.parametrize("wavelengths, words", [
    (["632 kg"], 'item "632 kg": cannot convert to m'),
    ([6e-7, -1e-7], "item -1e-07: must be > 0"),
    ([6e-7, True], "item true: not a number"),
    ("632 nm", "not a list"),
], ids=["mass_unit", "negative", "bool", "not_a_list"])
def test_bad_wavelength_items_are_named(wavelengths, words):
    """Each wavelength is checked on its own and a bad item is named in the message."""
    assert words in problems({"type": "lazuli", "wavelengths": wavelengths})


def test_negative_weight_is_refused():
    """A negative throughput item is refused (weights are >= 0)."""
    assert "item -1.0: must be >= 0" in problems({"type": "lazuli", "throughput": [-1.0]})


def test_too_few_pupil_samples_are_refused():
    """pupil_samples must be >= 16."""
    assert "pupil_samples = 8: must be >= 16" in problems({"type": "lazuli", "pupil_samples": 8})


@pytest.mark.parametrize("key", ["wavelengths", "throughput", "opd", "pupil_samples"])
def test_lazuli_keys_are_refused_under_a_vortex(key):
    """The lazuli-only keys are unknown keys in a vortex table."""
    assert "unknown key %s" % key in problems({"type": "vortex", key: "auto" if key != "opd" else ""})


@pytest.mark.parametrize("table, name", [({"type": "lazuli"}, "lazuli_c0"), ({"type": "lazuli", "charge": 2}, "lazuli_c2"),
                                         ({"type": "vortex"}, "vortex_c2"), ({"type": "lazuli", "name": "mine"}, "mine")])
def test_automatic_coronagraph_names(table, name):
    """Automatic names are <type>_c<charge> with the type's default charge (lazuli: 0)."""
    assert P.make({}, [], [table]).coronagraphs[0].name == name


@pytest.mark.parametrize("sources, size", [(["planet"], None), (["fits_image"], 511)], ids=["planet", "example_image"])
def test_lazuli_templates_build_their_optics(in_tmp, fresh_coronagraphs, sources, size):
    """Regression: new_config with a lazuli coronagraph writes an r'-band wavelength and, for the bundled 512-pixel
    example image, [grid] size = 511, so the lazuli optics are built from the file as written (no PSF computed)."""
    from toycoronagraph.params_docs import new_config
    text = new_config(sources, ["lazuli:0"], path="lazuli.toml")
    assert "wavelength = 6.33e-07" in text and (size is None or "\nsize = %d " % size in text)
    scene = Scene.from_file("lazuli.toml")
    scene.build()
    _, coro = coronagraph_for(scene.grid, scene.params.coronagraphs[0], scene.params)
    assert coro.optics.img_pixel % 2 == 1 and coro.optics.param("wavelengths") == (6.33e-7,)


def test_json_schema_accepts_lazuli_tables_and_keeps_vortex_keys():
    """Regression: the JSON schema merges the keys shared by both types and accepts text items with a unit."""
    jsonschema = pytest.importorskip("jsonschema")
    from toycoronagraph.params_docs import json_schema
    schema = json_schema()
    shared = schema["properties"]["coronagraph"]["items"]["properties"]
    for key in ("charge", "lyot_mask_size", "symmetric"):  # declared by both types: neither overwrites the other
        assert len(shared[key]["anyOf"]) == 2, key
    jsonschema.validate({"coronagraph": [{"type": "lazuli", "charge": 0, "wavelengths": ["632 nm", 6.5e-7],
                                          "throughput": "auto", "opd": "", "pupil_samples": 64},
                                         {"type": "vortex", "charge": 6, "lyot_mask_size": 0.9}]}, schema)
    for bad in [{"type": "vortex", "wavelengths": [6e-7]}, {"type": "lazuli", "charge": -2},
                {"type": "lazuli", "wavelengths": [True]}, {"type": "vortex", "charge": 0},
                {"type": "vortex", "charge": 1}, {"charge": 0}]:  # (no type: a vortex)
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate({"coronagraph": [bad]}, schema)


# ---------------------------------------------------------------- d) Optics.lazuli and its settings

def test_optics_lazuli_stores_sorted_tuples_and_normalised_weights():
    """Optics.lazuli: kind/module, sorted params, tuple wavelengths and weights summing to 1; hashable."""
    o = lazuli_optics(2, wavelengths=[6e-7, 6.5e-7], weights=[2, 6])
    assert (o.kind, o.module) == ("lazuli", "toycoronagraph.lazuli")
    names = [k for k, _ in o.params]
    assert names == sorted(names) and set(names) == {"charge", "lyot_mask_size", "opd", "opd_sha1", "pupil_samples",
                                                     "scene_diameter", "scene_wavelength", "wavelengths", "weights"}
    assert o.param("wavelengths") == (6e-7, 6.5e-7) and o.param("weights") == pytest.approx((0.25, 0.75))
    assert isinstance(o.param("weights"), tuple) and hash(o) == hash(lazuli_optics(2, wavelengths=(6e-7, 6.5e-7),
                                                                                   weights=(1, 3)))


@pytest.mark.parametrize("wavelengths, weights, with_opd", [((WAVELEN,), (1.0,), False), ((6e-7, 6.5e-7), (0.3, 0.7), False),
                                                             ((WAVELEN,), (1.0,), True)],
                         ids=["one_wavelength", "two_wavelengths", "opd_map"])
def test_optics_settings_round_trip(in_tmp, wavelengths, weights, with_opd):
    """Optics.from_settings(settings()) gives the same Optics, also after a .npz save (arrays back to tuples)."""
    if with_opd:
        np.save("opd.npy", astigmatism(32, 10.0))
    o = lazuli_optics(2, wavelengths=wavelengths, weights=weights, pupil=32, opd="opd.npy" if with_opd else "")
    assert Optics.from_settings(o.settings()) == o
    np.savez("settings.npz", **o.settings())
    with np.load("settings.npz") as saved:
        back = Optics.from_settings(dict(saved))
    assert back == o and isinstance(back.param("wavelengths"), tuple) and isinstance(back.param("charge"), int)
    assert isinstance(back.param("opd"), str) and isinstance(back.param("opd_sha1"), str)


def test_optics_lazuli_opd_stored_as_absolute_path_with_content_digest(in_tmp):
    """The OPD map is stored as an absolute path plus a SHA-1 of its content, so the PSF cache tag follows the content."""
    np.save("opd.npy", astigmatism(32, 10.0))
    first = lazuli_optics(pupil=32, opd="opd.npy")
    assert first.param("opd") == os.path.abspath("opd.npy") and re.fullmatch("[0-9a-f]{16}", first.param("opd_sha1"))
    np.save("opd.npy", astigmatism(32, 20.0))
    second = lazuli_optics(pupil=32, opd="opd.npy")
    assert second.param("opd_sha1") != first.param("opd_sha1") and second != first
    assert Coronagraph(first, "a").tag != Coronagraph(second, "a").tag
    assert lazuli_optics().param("opd_sha1") == ""


# ---------------------------------------------------------------- e) OPD maps and pupil samples

def test_read_opd_fits_in_microns_equals_npy_in_nm(in_tmp):
    """A FITS map with BUNIT = um and a .npy map in nm give the same (n, n) float64 OPD in nm; NaN -> 0."""
    opd = astigmatism(16, 20.0)
    opd[0, 0] = np.nan
    fits.PrimaryHDU(opd*1e-3, fits.Header([("BUNIT", "um")])).writeto("opd_um.fits")
    np.save("opd_nm.npy", opd)
    from_fits, from_npy = L.read_opd("opd_um.fits"), L.read_opd("opd_nm.npy")
    assert from_fits.shape == (16, 16) and from_fits.dtype == np.float64
    assert np.allclose(from_fits, from_npy, rtol=1e-12, atol=1e-12) and from_npy[0, 0] == 0


def test_read_opd_fits_without_bunit_is_in_nm(in_tmp):
    """A FITS map without BUNIT is read in nm."""
    opd = astigmatism(8, 5.0)
    fits.PrimaryHDU(opd).writeto("opd.fits")
    assert np.array_equal(L.read_opd("opd.fits"), opd)


@pytest.mark.parametrize("shape", [(4, 5), (4,), (2, 4, 4)])
def test_read_opd_refuses_maps_that_are_not_square(in_tmp, shape):
    """A map that is not a square 2-D array is a ValueError."""
    np.save("opd.npy", np.zeros(shape))
    with pytest.raises(ValueError, match="square 2-D"):
        L.read_opd("opd.npy")


def test_read_opd_refuses_a_bunit_that_is_not_a_length(in_tmp):
    """BUNIT = Jy cannot be converted to nm: ValueError."""
    fits.PrimaryHDU(np.zeros((4, 4)), fits.Header([("BUNIT", "Jy")])).writeto("opd.fits")
    with pytest.raises(ValueError):
        L.read_opd("opd.fits")


def test_pupil_samples_from_the_map_the_request_or_the_config(in_tmp):
    """pupil_samples: the OPD map's size if any, else the requested number, else the config's 512."""
    np.save("opd.npy", np.zeros((32, 32)))
    assert L.pupil_samples("") == 512 and L.pupil_samples("", 64) == 64
    assert L.pupil_samples("opd.npy") == 32 and L.pupil_samples("opd.npy", 32) == 32
    with pytest.raises(ValueError, match="pupil_samples = 64 but the OPD map"):
        L.pupil_samples("opd.npy", 64)


# ---------------------------------------------------------------- f) image geometry

def test_frame_and_block_shapes_and_dtypes():
    """frame() is (N, N) float64, block() (ny, nx) float64; total and throughput exist when the source is built."""
    src = lazuli_optics().backend().source(0.5, -0.3)
    assert src.frame().shape == (21, 21) and src.frame().dtype == np.float64
    assert src.block((0.2, 0.1), (3, 6)).shape == (3, 6) and src.block((0.2, 0.1), (3, 6)).dtype == np.float64
    assert isinstance(float(src.total), float) and isinstance(float(src.throughput), float)


@pytest.mark.parametrize("n", [15, 20, 21])
@pytest.mark.parametrize("d_col, d_row", [(3, 0), (0, 3), (-2, 4)])
def test_source_orientation_x_columns_y_rows(n, d_col, d_row):
    """A source at +x lands at larger columns, +y at larger rows (star at (N-1)/2), as Grid.pixel / arcsec_to_pixel."""
    o = lazuli_optics(n=n)
    s, c = o.pixel_scale, o.centre
    half = 0.5 if n % 2 == 0 else 0.0  # even N: the star is between pixels, so offsets of k + 1/2 land on one
    dx, dy = d_col+half, d_row+half    # in pixels
    frame = o.backend().source(dx*s, dy*s).frame()
    peak = np.unravel_index(np.argmax(frame), frame.shape)
    assert peak == (c+dy, c+dx)
    r, angle = np.hypot(dx, dy)*PIXSCALE, np.arctan2(dy, dx)
    grid = Grid(size=n, pixel_scale=PIXSCALE, wavelength=WAVELEN, diameter=LAZULI_D)
    assert np.allclose(grid.pixel(r, angle), peak) and np.allclose(arcsec_to_pixel((r, angle), PIXSCALE, n), peak)


@pytest.mark.parametrize("n", [15, 21])
def test_on_axis_frame_is_centred_on_the_middle_pixel(n):
    """The star's image peaks at ((N-1)/2, (N-1)/2) and is symmetric under flips and transposition."""
    frame = lazuli_optics(n=n).backend().source(0, 0).frame()
    c = (n-1)//2
    assert np.unravel_index(np.argmax(frame), frame.shape) == (c, c)
    for other in (frame[::-1], frame[:, ::-1], frame.T):
        assert np.abs(frame-other).max() < 1e-12*frame.max()


def test_even_image_size_puts_the_star_between_four_pixels():
    """Optics.lazuli itself on an even N follows the shared convention: the star at (N-1)/2, between pixels."""
    frame = lazuli_optics(n=20).backend().source(0, 0).frame()
    assert np.abs(frame-frame[::-1]).max() < 1e-12*frame.max() and np.abs(frame-frame[:, ::-1]).max() < 1e-12*frame.max()
    assert np.allclose(frame[9:11, 9:11], frame.max(), rtol=1e-12, atol=0)


def test_source_along_y_is_the_transpose_of_source_along_x():
    """Swapping x and y of a source transposes its image (rows = y, columns = x)."""
    b = lazuli_optics().backend()
    along_x, along_y = b.source(1.3, 0.0).frame(), b.source(0.0, 1.3).frame()
    assert np.abs(along_x.T-along_y).max() < 1e-12*along_x.max()


@pytest.mark.parametrize("n", [20, 21])
@pytest.mark.parametrize("row0, col0, shape", [(3, 2, (8, 12)), (5, 11, (9, 6)), (7, 6, (6, 8)), (10, 10, (1, 4)),
                                               (0, 0, None)], ids=["8x12", "9x6", "6x8", "1x4", "whole"])
def test_block_matches_the_frame_for_even_and_odd_shapes(n, row0, col0, shape):
    """block() on any (even or odd) shape equals the same pixels of frame(), for an even or odd image size."""
    o = lazuli_optics(2, n=n, pupil=32)
    s, c = o.pixel_scale, o.centre
    shape = shape or (n, n)
    src = o.backend().source(2.3*s, -1.6*s)
    frame = src.frame()
    center = ((col0+(shape[1]-1)/2-c)*s, (row0+(shape[0]-1)/2-c)*s)
    block = src.block(center, shape)
    assert np.abs(block-frame[row0:row0+shape[0], col0:col0+shape[1]]).max() < 1e-12*frame.max()


@pytest.mark.parametrize("shape", [(2, 2), (4, 4), (4, 6), (3, 4)])
def test_even_block_centred_on_the_star_is_symmetric(shape):
    """Regression (half-pixel offset): pixel centres of a block are symmetric about its centre for even shapes."""
    block = lazuli_optics().backend().source(0, 0).block((0.0, 0.0), shape)
    assert np.abs(block-block[::-1]).max() < 1e-12*block.max()
    assert np.abs(block-block[:, ::-1]).max() < 1e-12*block.max()


@pytest.mark.parametrize("shape", [(2, 4), (3, 2), (1, 1)])
def test_block_pixel_centres_sit_at_half_pixels_from_its_centre(shape):
    """Pixel (i, j) of block(center, (ny, nx)) is centred at center + ((j-(nx-1)/2) s, (i-(ny-1)/2) s), even shapes too."""
    o = lazuli_optics(2, pupil=32)
    s = o.pixel_scale
    src = o.backend().source(0.7, -0.3)
    (x0, y0), (ny, nx) = (0.1, 0.2), shape
    block = src.block((x0, y0), shape)
    single = np.array([[src.block((x0+(j-(nx-1)/2)*s, y0+(i-(ny-1)/2)*s), (1, 1))[0, 0] for j in range(nx)]
                       for i in range(ny)])
    assert np.abs(block-single).max() < 1e-12*block.max()


# ---------------------------------------------------------------- g) units and normalisation

@pytest.mark.parametrize("pupil", [33, 64, 65])
def test_charge0_total_and_throughput_are_one(pupil):
    """Without a coronagraph a source's total light (Parseval) and throughput are 1 (even or odd pupil sampling)."""
    src = lazuli_optics(pupil=pupil).backend().source(0.3, -0.2)
    assert abs(src.total-1) < 1e-12 and abs(src.throughput-1) < 1e-12


@pytest.mark.parametrize("diameter", [LAZULI_D, 2.4])
def test_charge0_peak_is_pi_p2_over_4(diameter):
    """The peak pixel is pi p^2/4 (p: pixel in Lazuli lambda/D) up to the gray edge (< 1/P); exact with it."""
    o = lazuli_optics(diameter=diameter)
    b = o.backend()
    peak = b.source(0, 0).frame().max()
    p = o.pixel_scale*LAZULI_D/diameter
    assert abs(peak/(np.pi*p**2/4)-1) < 1/PUPIL
    exact = (p/PUPIL)**2*b.aperture.sum()**2/b.entering  # pi p^2/4 for a pupil without the gray edge
    assert abs(peak/exact-1) < 1e-9


def test_one_period_of_the_focal_field_holds_all_the_light():
    """Over 2P+1 pixels of P/(2P+1) lambda/D (one period of the sampled focal field) the image sums to 1."""
    frame = lazuli_optics(n=2*PUPIL+1, psf_range=PUPIL/2).backend().source(0.3, -0.2).frame()
    assert abs(frame.sum()-1) < 1e-10


def test_total_with_a_vortex_is_the_light_over_one_focal_period():
    """Charge 2: total (Parseval on the exit field) equals the image summed over one focal-field period, and is < 1."""
    pupil = 32
    src = lazuli_optics(2, n=2*pupil+1, pupil=pupil, psf_range=pupil/2).backend().source(1.1, -0.6)
    assert abs(src.frame().sum()/src.total-1) < 1e-10 and 0 < src.total < 1


def test_light_in_a_growing_field_approaches_one():
    """The light missing from a +-r lambda/D field is positive, small and shrinks as r grows."""
    missing = []
    for half in (5, 10):
        n = int(2*half/0.25) | 1
        missing.append(1-lazuli_optics(n=n, psf_range=n*0.25/2).backend().source(0, 0).frame().sum())
    assert 0 < missing[1] < missing[0] < 0.05


def test_peak_scales_as_inverse_wavelength_squared():
    """The unocculted peak at lambda_1 over that at lambda_2 is (lambda_2/lambda_1)^2."""
    p600, p680 = [lazuli_optics(wavelengths=(lam,)).backend().source(0, 0).frame().max() for lam in (600e-9, 680e-9)]
    assert abs(p680/p600/(600/680)**2-1) < 1e-9


# ---------------------------------------------------------------- h) the vortex coronagraph

def test_charge2_on_axis_leak_is_tiny():
    """A charge-2 vortex leaves < 1e-4 of the star's light and a peak < 1e-6 of the unocculted one."""
    star = lazuli_optics(2, pupil=128).backend().source(0, 0)
    unocculted = lazuli_optics(0, pupil=128).backend().source(0, 0).frame().max()
    assert star.throughput < 1e-4 and star.frame().max() < 1e-6*unocculted


def test_charge2_throughput_rises_with_offset():
    """The charge-2 throughput increases with the source offset, from ~0.23 at 0.5 lambda/D to > 0.9 at 4."""
    b = lazuli_optics(2).backend()
    t = [b.source(r, 0).throughput for r in (0.5, 1, 2, 4)]
    assert np.all(np.diff(t) > 0) and 0.15 < t[0] < 0.3 and 0.9 < t[-1] <= 1


def test_charge2_throughput_agrees_with_the_hcipy_vortex():
    """Two independent codes: lazuli charge 2 (128 pupil samples) vs optics.VortexBackend, within 0.02."""
    ours, hcipy = lazuli_optics(2, pupil=128).backend(), Optics.vortex(2, 21, 21/SAMPLING/2).backend()
    for r in (0.5, 1, 2, 4):
        assert abs(ours.source(r, 0).throughput-float(hcipy.source(r, 0).throughput)) < 0.02, r


def test_lyot_stop_passes_its_area_fraction():
    """With a charge the reference light is the Lyot stop's share (~0.95^2); without, the whole aperture."""
    vortex, clear = lazuli_optics(2).backend(), lazuli_optics(0).backend()
    assert abs(vortex.lyot_light/vortex.entering-0.95**2) < 2/PUPIL
    assert clear.lyot_stop is None and clear.vortex is None and clear.lyot_light == clear.entering


# ---------------------------------------------------------------- i) symmetry flags and band limit

def test_rotationally_symmetric_only_without_an_opd_map(in_tmp):
    """rotationally_symmetric is True for the clear pupil and False with an OPD map."""
    np.save("opd.npy", astigmatism(32, 20.0))
    assert lazuli_optics(pupil=32).backend().rotationally_symmetric is True
    assert lazuli_optics(pupil=32, opd="opd.npy").backend().rotationally_symmetric is False


@pytest.mark.parametrize("diameter, wavelengths", [(LAZULI_D, (WAVELEN,)), (2.4, (WAVELEN,)), (2.4, (680e-9, 600e-9))])
def test_bandlimit_formula(diameter, wavelengths):
    """bandlimit = (D / scene diameter) (scene wavelength / shortest wavelength), cycles per scene lambda/D."""
    b = lazuli_optics(wavelengths=wavelengths, weights=(1.0,)*len(wavelengths), diameter=diameter).backend()
    assert b.bandlimit == pytest.approx(LAZULI_D/diameter*WAVELEN/min(wavelengths), rel=1e-12)


def test_method_choice_follows_the_symmetry_flags(in_tmp):
    """Clear pupil: circular, or spatial 'rotate'; with an OPD map: spatial 'direct' and 'rotate' is refused."""
    np.save("opd.npy", astigmatism(32, 20.0))
    image = np.ones((21, 21))
    target = (True, "target symmetric")
    clear = Coronagraph(lazuli_optics(pupil=32), "clear", num_cores=1)
    aberrated = Coronagraph(lazuli_optics(pupil=32, opd="opd.npy"), "aberrated", num_cores=1)
    assert clear.psf_symmetry() == (True, "PSFs symmetric: lazuli model")
    assert aberrated.psf_symmetry() == (False, "PSFs not symmetric: lazuli model")
    assert clear.choose(image, target)[:2] == ("circular", None)
    assert clear.choose(image, target, "spatial")[:2] == ("spatial", "rotate")
    assert aberrated.choose(image, target)[:2] == ("spatial", "direct")
    with pytest.raises(ValueError, match="rotate"):
        aberrated.choose(image, target, "spatial", "rotate")


# ---------------------------------------------------------------- j) several wavelengths

@pytest.mark.parametrize("charge", [0, 2])
def test_polychromatic_image_is_the_weighted_sum(charge):
    """With several wavelengths, frame, block, total and throughput are the weighted sums of one-wavelength ones."""
    lams, weights = (600e-9, 632e-9, 680e-9), (0.2, 0.5, 0.3)
    poly = lazuli_optics(charge, wavelengths=lams, weights=weights, pupil=32).backend().source(1.3, -0.4)
    singles = [lazuli_optics(charge, wavelengths=(lam,), pupil=32).backend().source(1.3, -0.4) for lam in lams]
    frame = poly.frame()
    assert np.abs(frame-sum(w*s.frame() for w, s in zip(weights, singles))).max() < 1e-12*frame.max()
    block = poly.block((0.25, -0.5), (4, 5))
    assert np.abs(block-sum(w*s.block((0.25, -0.5), (4, 5)) for w, s in zip(weights, singles))).max() < 1e-12*frame.max()
    assert poly.total == pytest.approx(sum(w*s.total for w, s in zip(weights, singles)), rel=1e-12)
    assert poly.throughput == pytest.approx(sum(w*s.throughput for w, s in zip(weights, singles)), rel=1e-12)


# ---------------------------------------------------------------- k) builders, odd sizes, cache names

def test_every_coronagraph_type_has_an_optics_builder():
    """OPTICS_BUILDERS has a builder for each declared [[coronagraph]] type."""
    assert set(OPTICS_BUILDERS) == set(P.CORONAGRAPHS)


def test_optics_for_lazuli_builds_lazuli_optics():
    """optics_for("lazuli", ...) gives kind 'lazuli' on the image size, half-width and scene wavelength/diameter."""
    o = optics_for("lazuli", table_values({"type": "lazuli", "charge": 2, "pupil_samples": PUPIL}), 21, 2.625, WAVELEN, 2.4)
    assert (o.kind, o.img_pixel, o.psf_range, o.param("charge"), o.param("pupil_samples")) == ("lazuli", 21, 2.625, 2, PUPIL)
    assert (o.param("scene_wavelength"), o.param("scene_diameter")) == (WAVELEN, 2.4)


@pytest.mark.parametrize("n", [20, 22])
def test_even_image_size_is_refused_for_lazuli(n):
    """An even image size is a ParamError naming the two odd sizes around it in [grid] size."""
    with pytest.raises(P.ParamError) as err:
        optics_for("lazuli", table_values({"type": "lazuli", "pupil_samples": PUPIL}), n, 2.5, WAVELEN, LAZULI_D)
    assert "odd image size (got %d)" % n in str(err.value) and "[grid] size = %d or %d" % (n-1, n+1) in str(err.value)


def test_even_image_size_is_fine_for_the_vortex():
    """The vortex builder accepts an even size (only lazuli needs an odd one)."""
    assert optics_for("vortex", table_values({"type": "vortex"}), 20, 2.5, WAVELEN, LAZULI_D).img_pixel == 20


def test_pupil_sampling_too_coarse_for_the_field_is_refused(in_tmp):
    """The focal field repeats every pupil_samples lambda/D: on 101 pixels of 0.5 lambda/D (PSF frames reaching
    62.5 lambda/D from a source) 62 samples, or a 62 x 62 OPD map, are refused naming 63; 63 samples work."""
    values = table_values({"type": "lazuli", "pupil_samples": 62})
    with pytest.raises(P.ParamError, match="repeats every 62 lambda/D.*set pupil_samples >= 63"):
        optics_for("lazuli", values, 101, 25.25, WAVELEN, LAZULI_D)
    np.save("opd.npy", np.zeros((62, 62)))
    with pytest.raises(P.ParamError, match="OPD map of at least 63 x 63 samples"):
        optics_for("lazuli", dict(values, opd="opd.npy", pupil_samples="auto"), 101, 25.25, WAVELEN, LAZULI_D)
    assert optics_for("lazuli", dict(values, pupil_samples=63), 101, 25.25, WAVELEN, LAZULI_D).param("pupil_samples") == 63


def test_pupil_sampling_period_is_in_scene_lambda_d_at_the_shortest_wavelength():
    """The period is pupil_samples Lazuli lambda/D at the shortest wavelength, in scene lambda/D: 64 samples at
    600 nm on a 2.4 m scene at 632 nm repeat every 64 (600/3.0)/(632/2.4) = 48.6 scene lambda/D."""
    values = table_values({"type": "lazuli", "wavelengths": ["600 nm", "680 nm"], "pupil_samples": 64})
    with pytest.raises(P.ParamError, match="repeats every 48.61 lambda/D \\(at 600 nm\\)"):
        optics_for("lazuli", values, 81, 20.25, WAVELEN, 2.4)


def test_coronagraph_for_reads_type_and_charge_strings(in_tmp, fresh_coronagraphs):
    """coronagraph_for(grid, "lazuli:2") gives lazuli_c2 with the config's 512 pupil samples (no PSF computed)."""
    scene = ring_scene(21)
    scene.build()
    table, coro = coronagraph_for(scene.grid, "lazuli:2", scene.params)
    assert (table.name, coro.optics.kind, coro.optics.param("charge")) == ("lazuli_c2", "lazuli", 2)
    assert coro.optics.param("pupil_samples") == 512 and coro.optics.param("wavelengths") == (WAVELEN,)
    assert coro.computed == []


@pytest.mark.parametrize("symmetric", [True, False])
def test_symmetric_setting_of_a_lazuli_table_overrides_the_model(in_tmp, fresh_coronagraphs, symmetric):
    """[[coronagraph]] symmetric = true/false reaches the Coronagraph and wins over the backend's declaration."""
    scene = ring_scene(21)
    scene.build()
    _, coro = coronagraph_for(scene.grid, {"type": "lazuli", "symmetric": symmetric, "pupil_samples": 32}, scene.params)
    sym, why = coro.psf_symmetry()
    assert coro.symmetric is symmetric and sym is symmetric and "set in the parameters" in why


def test_coronagraph_tag_keeps_the_vortex_format_and_hashes_lazuli():
    """Cache tags: vortex unchanged (vortex_c2_N21_range..._lyot0.95), lazuli <kind>_N_range_<12 hex digits>."""
    assert Coronagraph(Optics.vortex(2, 21, 2.625), "v").tag == "vortex_c2_N21_range2.625_lyot0.95"
    tag = Coronagraph(lazuli_optics(2), "l").tag
    assert re.fullmatch(r"lazuli_N21_range2\.625_[0-9a-f]{12}", tag)
    other = Coronagraph(lazuli_optics(2, wavelengths=(6e-7, 6.5e-7), weights=(1, 1)), "l").tag
    assert other != tag and other.startswith("lazuli_N21_range2.625_")


def test_same_settings_tells_close_wavelengths_apart():
    """Regression: cached settings compare tuples to a relative 1e-5 (632 vs 633 nm differ) and by shape."""
    cached = {"optics_wavelengths": np.array([632e-9]), "optics_opd": np.array(""), "psf_range": np.array(2.625)}
    assert psf._same_settings(cached, {"optics_wavelengths": (632e-9,), "optics_opd": "", "psf_range": 2.625})
    assert not psf._same_settings(cached, {"optics_wavelengths": (633e-9,)})
    assert not psf._same_settings(cached, {"optics_wavelengths": (632e-9, 650e-9)})
    assert not psf._same_settings(cached, {"optics_opd": "/elsewhere/opd.npy"})


# ---------------------------------------------------------------- l) observe()

def test_observe_refuses_an_even_grid_before_any_psf(in_tmp, fresh_coronagraphs):
    """observe() on an even grid with a vortex and a lazuli coronagraph stops before computing a PSF."""
    scene = ring_scene(22)
    with pytest.raises(P.ParamError, match=r"odd image size \(got 22\)"):
        observe(scene, ["vortex:2", {"type": "lazuli", "pupil_samples": PUPIL}], verbose=False)
    assert not os.path.exists("cache")


def test_observe_small_ring_through_lazuli(in_tmp, fresh_coronagraphs):
    """A ring through lazuli_c0 and lazuli_c2: (N, N) images, circular method, the coronagraph removes light."""
    scene = ring_scene(21)
    res = observe(scene, [{"type": "lazuli", "pupil_samples": 32}, {"type": "lazuli", "charge": 2, "pupil_samples": 32}],
                  verbose=False)
    assert set(res.finals) == {"lazuli_c0", "lazuli_c2"}
    clear, vortex = res.final("lazuli_c0"), res.final("lazuli_c2")
    assert clear.shape == (21, 21) and clear.dtype == np.float64 and np.all(np.isfinite(clear))
    assert res.methods["lazuli_c0", "ring"][0] == "circular" and res.methods["lazuli_c2", "ring"][0] == "circular"
    assert 0.8 < clear.sum()/scene.total.sum() < 1 and vortex.sum() < clear.sum()
    assert os.path.basename(res.coronagraphs["lazuli_c2"].radial_psfs_file()).startswith("radial_lazuli_N21_")


def test_saved_fits_header_writes_wavelength_lists_as_text(in_tmp, fresh_coronagraphs):
    """Result.save: the coronagraph FITS header names the optics and lists wavelengths and weights as text."""
    res = observe(ring_scene(15), [{"type": "lazuli", "wavelengths": ["620 nm", "0.645 um"], "throughput": [1, 3],
                                    "pupil_samples": 32}], verbose=False)
    folder = res.save("run", overwrite=True)
    header = fits.getheader(os.path.join(folder, "lazuli_c0.fits"))
    assert (header["OPTICS"], header["CHARGE"], header["PUPIL_SA"]) == ("lazuli", 0, 32)
    assert (header["WAVELENG"], header["WEIGHTS"]) == ("6.2e-07, 6.45e-07", "0.25, 0.75")


# ---------------------------------------------------------------- m) the older layout and Target

def test_legacy_lazuli_charge_defaults_to_zero():
    """Regression: an older-layout [coronagraph] type = "lazuli" without charge has charge 0 (vortex: 2)."""
    assert C.from_dict({"coronagraph": {"type": "lazuli"}}).charge == 0
    assert C.from_dict({"coronagraph": {"type": "vortex"}}).charge == 2


def test_legacy_lazuli_reads_its_keys_with_units():
    """The older layout reads the lazuli keys like the scene layout (unit strings converted to m)."""
    cfg = C.from_dict({"coronagraph": {"type": "lazuli", "charge": 2, "wavelengths": ["632 nm", 6.5e-7],
                                       "pupil_samples": 64}})
    assert (cfg.coronagraph_type, cfg.charge, cfg.throughput, cfg.opd, cfg.pupil_samples) == ("lazuli", 2, "auto", "", 64)
    assert np.allclose(cfg.wavelengths, [632e-9, 6.5e-7], rtol=1e-12, atol=0)


@pytest.mark.parametrize("coronagraph, words", [
    ({"type": "lazuli", "charge": 3}, "charge = 3: a vortex charge must be even"),
    ({"type": "lazuli", "charge": -2}, "charge = -2: must be >= 0"),
    ({"type": "lazuli", "wavelengths": [6e-7, 6.5e-7], "throughput": [1.0]}, "throughput has 1 values for 2 wavelength(s)"),
    ({"type": "lazuli", "wavelengths": ["632 kg"]}, "cannot convert to m"),
    ({"type": "vortex", "charge": 0}, "a vortex charge must be positive and even"),
    ({"type": "starshade"}, "must be one of"),
], ids=["odd_charge", "negative_charge", "throughput_length", "mass_unit", "vortex_charge0", "unknown_type"])
def test_legacy_coronagraph_checks(coronagraph, words):
    """Regression: the older layout applies the same cross-key checks to a lazuli [coronagraph] (ValueError)."""
    with pytest.raises(ValueError) as err:
        C.from_dict({"coronagraph": coronagraph})
    assert words in str(err.value)


def test_legacy_lazuli_file_maps_to_a_lazuli_coronagraph(in_tmp):
    """params.load of an older-layout lazuli file gives a lazuli [[coronagraph]] with its keys (default name lazuli_c0)."""
    legacy_toml("old.toml", coronagraph='type = "lazuli"\nwavelengths = ["632 nm"]\npupil_samples = 32\n')
    c = P.load("old.toml").coronagraphs[0]
    assert (c.type, c.name, c.values["pupil_samples"]) == ("lazuli", "lazuli_c0", 32)
    assert np.allclose(c.values["wavelengths"], [632e-9], rtol=1e-12, atol=0)


def test_example_parameter_file_loads_in_the_scene_layout(in_tmp):
    """Regression: config.write_example's older-layout file is read by params.load and migrated (as before stage 7)."""
    from toycoronagraph.params_docs import migrate
    C.write_example("example.toml")
    assert [c.name for c in P.load("example.toml").coronagraphs] == ["vortex_c2"]
    assert "[[coronagraph]]" in migrate("example.toml") and "\n# wavelengths" in C.example_text()


@pytest.mark.parametrize("key, value", [("wavelengths", ["600 nm"]), ("opd", "map.fits"), ("pupil_samples", 64)])
def test_legacy_lazuli_keys_are_refused_under_a_vortex(key, value):
    """Regression: the older layout refuses a lazuli-only key under a vortex (as the scene layout does)."""
    with pytest.raises(ValueError, match='%s: only used with type = "lazuli"' % key):
        C.from_dict({"coronagraph": {"type": "vortex", key: value}})


def test_target_with_legacy_lazuli_file_builds_lazuli_optics(in_tmp, fresh_coronagraphs):
    """Regression: Target with [coronagraph] type = "lazuli" builds lazuli optics (it used to build a vortex)."""
    import toycoronagraph.main as toy
    write_image("ring.fits", ring_image(15))
    legacy_toml("old.toml", coronagraph='type = "lazuli"\ncharge = 2\npupil_samples = 32\n')
    t = toy.Target("old.toml")
    assert (t.par.coronagraph_type, t.par.charge, t.par.target_pixel) == ("lazuli", 2, 15)
    optics = t._optics(t.par.charge)
    assert (optics.kind, optics.param("charge"), optics.param("pupil_samples")) == ("lazuli", 2, 32)
    assert t._charge(0) == 0 and t._charge(3) == 2  # 0 = the telescope alone; an odd charge falls back


@pytest.mark.parametrize("layout, words", [
    ("old", "[target] target_pixel"),
    ("new", "[grid] size = 19 or 21"),
])
def test_target_refuses_an_even_image(in_tmp, fresh_coronagraphs, layout, words):
    """Target with a lazuli coronagraph and an even image raises when it is built, with the layout's fix."""
    import toycoronagraph.main as toy
    write_image("even.fits", ring_image(20))
    if layout == "old":
        legacy_toml("params.toml", image="even")
    else:
        with open("params.toml", "w") as f:
            f.write('[observation]\ndiameter = %g\n[[source]]\nkind = "fits_image"\nfile = "even"\n[[coronagraph]]\n'
                    'type = "lazuli"\npupil_samples = 32\n' % LAZULI_D)
    with pytest.raises(P.ParamError) as err:
        toy.Target("params.toml")
    assert "odd image size (got 20)" in str(err.value) and words in str(err.value)


def test_target_follows_the_even_size_fix_of_the_message(in_tmp, fresh_coronagraphs):
    """Scene layout: [grid] size = 19 on a 20-pixel image (the error's advice) resamples it and Target builds lazuli optics."""
    import toycoronagraph.main as toy
    write_image("even.fits", ring_image(20))
    with open("params.toml", "w") as f:
        f.write('[observation]\ndiameter = %g\n[grid]\nsize = 19\n[[source]]\nkind = "fits_image"\nfile = "even"\n'
                '[[coronagraph]]\ntype = "lazuli"\npupil_samples = 32\n[compute]\nnum_cores = 1\n' % LAZULI_D)
    t = toy.Target("params.toml")
    optics = t._optics(t.par.charge)
    assert (t.par.target_pixel, t.pre_img.shape, optics.kind, optics.img_pixel) == (19, (19, 19), "lazuli", 19)


def test_target_legacy_lazuli_final_matches_observe(in_tmp, fresh_coronagraphs):
    """A tiny Target run (older layout, type = "lazuli") gives observe()'s image and saves lazuli_charge2 plots."""
    import toycoronagraph.main as toy
    write_image("ring.fits", ring_image(15))
    legacy_toml("old.toml", coronagraph='type = "lazuli"\ncharge = 2\npupil_samples = 32\n')
    final = toy.Target("old.toml").plot_final(add_planet=False, method="circular", plot_dpi=30)
    assert os.path.exists("ring_lazuli_charge2_final.png") and not os.path.exists("ring_charge2_final.png")
    scene = Scene(diameter=LAZULI_D, compute={"psf_cache_dir": "cache", "num_cores": 1})
    scene.add(FitsImage(file="ring.fits"))
    res = observe(scene, [{"type": "lazuli", "charge": 2, "pupil_samples": 32}], method="circular", verbose=False)
    reference = res.final("lazuli_c2")
    assert final.shape == (15, 15) and np.abs(final-reference).max() <= 1e-12*reference.max()


@pytest.mark.parametrize("charge, sentence", [
    (2, "This is the Lazuli telescope with a charge-2 vortex coronagraph."),
    (0, "This is the Lazuli telescope without a coronagraph."),
])
def test_target_contrast_names_the_lazuli_telescope(in_tmp, fresh_coronagraphs, capsys, charge, sentence):
    """Target.contrast says which Lazuli configuration it measures."""
    import toycoronagraph.main as toy
    write_image("ring.fits", ring_image(15))
    legacy_toml("old.toml", coronagraph='type = "lazuli"\ncharge = %d\npupil_samples = 32\n' % charge)
    t = toy.Target("old.toml")
    t.add_planet(pos=[4*PIXSCALE, 35], brightness=1e-4, mode="polar")
    capsys.readouterr()
    t.contrast(plot_contrast=False, method="circular")
    assert sentence in capsys.readouterr().out


# ---------------------------------------------------------------- n) without prysm

def test_package_works_without_prysm(in_tmp):
    """With prysm and the config packages missing, toycoronagraph and the vortex work; a lazuli table needs them."""
    code = textwrap.dedent('''
        import sys
        for name in ("prysm", "config_um", "config_um_wcc"):
            sys.modules[name] = None  # "import prysm" now raises ImportError
        import toycoronagraph, toycoronagraph.observe as O, toycoronagraph.main, toycoronagraph.config
        from toycoronagraph import params as P
        from toycoronagraph.optics import Optics
        assert "toycoronagraph.lazuli" not in sys.modules
        assert P.make({}, [], [{"type": "lazuli", "charge": 2}]).coronagraphs[0].name == "lazuli_c2"
        assert Optics.vortex(2, 15, 2.0).backend().source(0, 0).frame().shape == (15, 15)
        try:
            O.optics_for("lazuli", {q.name: q.default for q in P.CORONAGRAPHS["lazuli"][1]}, 21, 2.6, 632e-9, 3.0)
        except ImportError:
            print("lazuli needs prysm")
    ''')
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([ROOT] + [p for p in [os.environ.get("PYTHONPATH")] if p]),
               PYTHONDONTWRITEBYTECODE="1")  # write nothing outside the tmp folder
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, timeout=300)
    assert run.returncode == 0, run.stderr
    assert "lazuli needs prysm" in run.stdout
