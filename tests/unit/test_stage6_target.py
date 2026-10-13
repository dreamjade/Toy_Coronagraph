"""
Stage 6 unit tests: toycoronagraph/main.py Target end to end on tiny images - reading the target through
the parameter file (the .toml and the legacy .py layout of toycoronagraph/config.py), planets (add / delete /
move / list, integer positions), orbit plots, the method choice (symmetry, parameters, argument), radial PSF
caches, final images, their figures (+y up) and file names, planet orientation in both methods
(= scene.Grid.pixel), contrast, and the lossless VP9 planet video.

Images are N = 32 (even) and 33 (odd) pixels at 4 px per lambda/D, so every PSF is computed in about a
second; everything is written in each test's temporary folder (in_tmp). The comparisons against
brute-force exact images stay in tests/stage6_target.py.
"""
import io
import os
import shutil
import subprocess
import numpy as np
import pytest
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Circle
from astropy.io import fits
import toycoronagraph.main as toy
from toycoronagraph.config import RESOLVED_MARK
from toycoronagraph.observe import Coronagraph
from toycoronagraph.planet import planet_position
from toycoronagraph.psf import extended_frame_size, radial_source_positions
from toycoronagraph.scene import Grid

WAVELEN, DIAMETER = 0.8e-6, 2.4
PIXSCALE = WAVELEN/DIAMETER*206264.806247/4  # arcsec per pixel: 4 px per lambda/D
SIZES = [32, 33]
DPI = 20  # small figures: the tests check files and arrays, not looks

pytestmark = pytest.mark.filterwarnings("ignore:FigureCanvasAgg is non-interactive")
sizes = pytest.mark.parametrize("n", SIZES, ids=["even", "odd"])


@pytest.fixture(scope="module", autouse=True)
def module_coronagraphs():
    """Share Coronagraph objects (their in-memory spatial models and final images) among this module's tests
    only, so each spatial model is built once; files still go to each test's own folder."""
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(Coronagraph, "_instances", {})
        yield


@pytest.fixture(autouse=True)
def close_figures():
    """Close the figures a test made."""
    yield
    plt.close("all")


def images(n):
    """The test targets of size n: a symmetric 'ring' and a non-symmetric 'arc' (Jy/arcsec^2)."""
    c = (n-1)/2
    rows, cols = np.indices((n, n))
    rr = np.hypot(rows-c, cols-c)
    ring = np.exp(-0.5*((rr-0.25*n)/(0.05*n))**2)+0.3*np.exp(-0.5*((rr-0.1*n)/(0.03*n))**2)
    arc = ring*(cols > c)+0.5*np.exp(-0.5*((rows-c-0.2*n)**2+(cols-c+0.15*n)**2)/(0.04*n)**2)
    return {"ring": ring, "arc": arc}


def write_images(folder, n):
    """Write ring.fits and arc.fits (n x n, Jy/arcsec^2, with BUNIT, PIXSCALE, WAVELEN) in folder."""
    for key, img in images(n).items():
        hdr = fits.Header([("BUNIT", "Jy/arcsec^2"), ("PIXSCALE", PIXSCALE), ("WAVELEN", WAVELEN)])
        fits.PrimaryHDU(img, hdr).writeto(str(folder/(key+".fits")), overwrite=True)


def new_target(folder, n, name="ring", target="", sections="", para_file="toycoronagraph_para.toml"):
    """Write ring.fits and arc.fits (see write_images) and a parameter file reading <name>.fits with extra
    [target] lines and extra sections (num_cores = 2); return Target() read from it."""
    write_images(folder, n)
    text = '[target]\nfits_filename = "%s"\n%s\n%s\n' % (name, target, sections)
    if "[compute]" not in sections:
        text += "[compute]\nnum_cores = 2\n"
    (folder/para_file).write_text(text)
    return toy.Target(None if para_file == "toycoronagraph_para.toml" else para_file)


def peak(img):
    """(row, col) of the brightest pixel."""
    return np.array(np.unravel_index(np.argmax(img), img.shape))


def planet_r(n):
    """A planet radius (arcsec) well outside the IWA and inside the image: 0.3 N pixels."""
    return 0.3*n*PIXSCALE


# ---------------------------------------------------------------- reading the target

@sizes
def test_target_reads_the_image_through_the_parameter_file(n, in_tmp):
    """Target reads <fits_filename>.fits in Jy/arcsec^2, resolves the parameters, keeps the file unchanged."""
    t = new_target(in_tmp, n)
    text = (in_tmp/"toycoronagraph_para.toml").read_text()
    assert t.pre_img.dtype == np.float64 and np.allclose(t.pre_img, images(n)["ring"])
    assert t.par.target_pixel == n and np.isclose(t.par.psf_scale, PIXSCALE) and t.par.wavelength == WAVELEN
    assert t.file_name == "ring.fits" and t.output_prefix == "ring_" and t.scene is None
    assert t.planets == [] and t.orbits == [] and t.planets_brightness == []
    assert (in_tmp/"toycoronagraph_para.toml").read_text() == text
    assert RESOLVED_MARK in (in_tmp/"toycoronagraph_para_used.toml").read_text()


def test_target_reads_a_named_parameter_file(in_tmp):
    """Target(para_file) reads that file and saves the values used next to it as <name>_used.toml."""
    t = new_target(in_tmp, 32, name="arc", para_file="other.toml")
    assert t.par.source == "other.toml" and np.allclose(t.pre_img, images(32)["arc"])
    assert (in_tmp/"other_used.toml").exists() and not (in_tmp/"toycoronagraph_para_used.toml").exists()


def test_target_reads_a_legacy_py_parameter_file(in_tmp, capsys):
    """Without a .toml, Target reads an old toycoronagraph_para.py (with a notice), keeps its sampling and the file."""
    write_images(in_tmp, 32)
    lambda_d = WAVELEN/DIAMETER*206264.806247
    legacy = ('fits_filename = "ring"\ntarget_pixel = 32\nwavelength = %r\nD = %r\nvisual_range = %r\nnum_cores = 2\n'
              % (WAVELEN, DIAMETER, 32*PIXSCALE/lambda_d))
    (in_tmp/"toycoronagraph_para.py").write_text(legacy)
    t = toy.Target()
    assert "old-style toycoronagraph_para.py" in capsys.readouterr().out
    assert t.par.legacy is True and np.isclose(t.par.psf_scale, PIXSCALE) and t.par.target_pixel == 32
    assert np.allclose(t.pre_img, images(32)["ring"]) and t.output_prefix == "ring_"
    assert (in_tmp/"toycoronagraph_para.py").read_text() == legacy
    assert RESOLVED_MARK in (in_tmp/"toycoronagraph_para_used.toml").read_text()
    assert not (in_tmp/"toycoronagraph_para.toml").exists()


def test_target_refuses_an_image_without_a_unit(in_tmp):
    """An image without BUNIT (and input_unit "auto") cannot be put in Jy/arcsec^2: ValueError naming input_unit."""
    fits.PrimaryHDU(np.ones((16, 16))).writeto(str(in_tmp/"bare.fits"))
    (in_tmp/"toycoronagraph_para.toml").write_text('[target]\nfits_filename = "bare"\n')
    with pytest.raises(ValueError, match="input_unit"):
        toy.Target()


def test_example_target_has_no_output_prefix(in_tmp):
    """The bundled example (fits_filename "") is read at 512 x 512 and its outputs have no name in front."""
    (in_tmp/"toycoronagraph_para.toml").write_text('[target]\nfits_filename = ""\n')
    t = toy.Target()
    assert t.pre_img.shape == (512, 512) and t.output_prefix == ""


@sizes
def test_pixel_coordinates_are_centred_on_the_star(n, in_tmp):
    """xpix / ypix are pixel-centre offsets from the star at (N-1)/2 in arcsec; extent spans the pixel edges."""
    t = new_target(in_tmp, n)
    expected = (np.arange(n)-(n-1)/2)*PIXSCALE
    assert np.allclose(t.xpix, expected) and np.allclose(t.ypix, t.xpix)
    assert np.allclose(t.xpix, -t.xpix[::-1])
    assert np.allclose(t.extent, [-n/2*PIXSCALE, n/2*PIXSCALE, -n/2*PIXSCALE, n/2*PIXSCALE])


# ---------------------------------------------------------------- planets

def test_add_planet_polar_cartesian_and_moving(in_tmp):
    """Planets are stored as (r arcsec, angle rad), brightness per pixel area (Jy/psf_scale^2), orbit or None."""
    t = new_target(in_tmp, 32)
    t.add_planet([0.2, 30.0], 1e-4, mode="polar")
    t.add_planet([0.1, -0.1], 2e-4, mode="cartesian")
    t.add_planet([0.2, 0.1, 40.0, 30.0, 0.25], 3e-4, mode="moving")
    assert np.allclose(t.planets[0], [0.2, np.radians(30)])
    assert np.allclose(t.planets[1], [np.hypot(0.1, 0.1), -np.pi/4])
    assert np.allclose(t.planets[2], planet_position(0.2, 0.1, 40.0, 30.0, 0.25, mode="polar", res=1000))
    assert np.allclose(t.planets_brightness, np.array([1e-4, 2e-4, 3e-4])/PIXSCALE**2)
    assert t.orbits[0] is None and t.orbits[1] is None and np.allclose(t.orbits[2], [0.2, 0.1, 40, 30, 0.25])


@pytest.mark.parametrize("mode, pos, expected", [
    ("polar", [1, 35], [1.0, np.radians(35)]),
    ("cartesian", [1, 1], [np.sqrt(2), np.pi/4]),
    ("moving", [1, 0, 40, 30, 0], None),
    ("polar", [np.float32(1.0), np.int32(35)], [1.0, np.radians(35)]),
])
def test_add_planet_accepts_integer_positions(mode, pos, expected, in_tmp):
    """Regression: integer and numpy positions of any precision are accepted and not truncated (angle
    conversion, orbit time stored as floats)."""
    t = new_target(in_tmp, 32)
    t.add_planet(pos, 1e-4, mode=mode)
    assert len(t.planets) == 1
    if expected is not None:
        assert np.allclose(t.planets[0], expected)
    else:
        assert np.allclose(t.planets[0], planet_position(1.0, 0.0, 40.0, 30.0, 0.0, mode="polar", res=1000))
        t.planet_move(0.25, message=False)
        assert t.orbits[0][4] == 0.25


@pytest.mark.parametrize("pos, brightness, mode, words", [
    ([0.2, 30.0], 0.0, "polar", "Brightness is invalid"),
    ([0.2, 30.0], -1e-4, "polar", "Brightness is invalid"),
    ([0.2, 30.0], "1e-4", "polar", "Brightness is invalid"),
    ([0.0, 30.0], 1e-4, "polar", "Planet position is not allowed"),
    ([0.2, 30.0, 1.0], 1e-4, "polar", "Planet position is not allowed"),
    ([0.1], 1e-4, "cartesian", "planet position is not allowed"),
    ([0.2, 1.0, 40.0, 30.0, 0.0], 1e-4, "moving", "Planet position is not allowed"),
    ([0.2, 0.1, 40.0, 30.0], 1e-4, "moving", "five parameters"),
    ([0.2, 30.0], 1e-4, "spiral", "no spiral mode"),
])
def test_add_planet_rejects_invalid_input_with_a_message(pos, brightness, mode, words, in_tmp, capsys):
    """Invalid positions, brightnesses or modes print a message and add nothing."""
    t = new_target(in_tmp, 32)
    capsys.readouterr()
    t.add_planet(pos, brightness, mode=mode)
    assert words in capsys.readouterr().out
    assert t.planets == [] and t.orbits == [] and t.planets_brightness == []


@pytest.mark.parametrize("brightness", [np.float32(1e-4), np.int64(1)], ids=["float32", "int64"])
def test_add_planet_accepts_numpy_scalar_brightness(brightness, in_tmp):
    """A positive numpy scalar brightness (e.g. taken from an array) is accepted like a Python number."""
    t = new_target(in_tmp, 32)
    t.add_planet([0.2, 30.0], brightness, mode="polar")
    assert len(t.planets) == 1 and np.isclose(t.planets_brightness[0], float(brightness)/PIXSCALE**2)


@pytest.mark.parametrize("pos, mode", [
    ([0.2, None], "polar"),
    ([None, 0.1], "cartesian"),
    ([0.2, "north"], "polar"),
    ([0.2, 0.1, 40.0, 30.0, None], "moving"),
    ([0.2, np.nan], "polar"),
    (None, "polar"),
], ids=["polar None", "cartesian None", "polar text", "moving None", "polar NaN", "no position"])
def test_add_planet_rejects_non_numeric_positions_with_a_message(pos, mode, in_tmp, capsys):
    """A position entry that is not a number is reported ('not allowed') and nothing is added."""
    t = new_target(in_tmp, 32)
    capsys.readouterr()
    t.add_planet(pos, 1e-4, mode=mode)
    assert "position is not allowed" in capsys.readouterr().out
    assert t.planets == [] and t.orbits == [] and t.planets_brightness == []


def _three_planets(t):
    """Add a static polar, a static cartesian and a moving planet to t."""
    t.add_planet([0.2, 30.0], 1e-4, mode="polar")
    t.add_planet([0.1, -0.1], 2e-4, mode="cartesian")
    t.add_planet([0.2, 0.1, 40.0, 30.0, 0.25], 3e-4, mode="moving")


def test_delete_planet_removes_the_given_order_everywhere(in_tmp, capsys):
    """delete_planet(2) removes the second planet's position, brightness and orbit, then lists the rest."""
    t = new_target(in_tmp, 32)
    _three_planets(t)
    first, third = t.planets[0], t.planets[2]
    capsys.readouterr()
    t.delete_planet(2)
    out = capsys.readouterr().out
    assert "Successfully remove planet #2" in out and "Moving Planet 2" in out
    assert len(t.planets) == len(t.orbits) == len(t.planets_brightness) == 2
    assert t.planets[0] is first and t.planets[1] is third
    assert np.allclose(t.planets_brightness, np.array([1e-4, 3e-4])/PIXSCALE**2) and t.orbits[0] is None


@pytest.mark.parametrize("order", [0, 4, -1, 1.0, "1"])
def test_delete_planet_with_an_invalid_order_keeps_every_planet(order, in_tmp, capsys):
    """An order that is not 1..len(planets) (or not an int) prints a message and deletes nothing."""
    t = new_target(in_tmp, 32)
    _three_planets(t)
    capsys.readouterr()
    t.delete_planet(order)
    assert "Input order number is not existed" in capsys.readouterr().out and len(t.planets) == 3


def test_list_planets_reports_kind_position_and_brightness_in_jy(in_tmp, capsys):
    """list_planets prints each planet as Static/Moving with (x, y) in arcsec and its brightness in Jy."""
    t = new_target(in_tmp, 32)
    capsys.readouterr()
    t.list_planets()
    assert "There is no planet" in capsys.readouterr().out
    t.add_planet([0.3, 0.0], 1e-4, mode="cartesian")
    t.list_planets()
    out = capsys.readouterr().out
    assert out.startswith("Static Planet 1: (0.3, 0.0) arcsec") and "brightness: 1.00e-04 Jy" in out


def test_planet_move_specific_and_cumulate(in_tmp):
    """'specific' sets the orbital time, 'cumulate' adds to it; the position follows the orbit."""
    t = new_target(in_tmp, 32)
    t.add_planet([0.2, 0.1, 40.0, 30.0, 0.0], 1e-4, mode="moving")
    t.planet_move(0.3, mode="specific", message=False)
    assert t.orbits[0][4] == 0.3
    t.planet_move(0.2, mode="cumulate", message=False)
    assert np.isclose(t.orbits[0][4], 0.5)
    assert np.allclose(t.planets[0], planet_position(0.2, 0.1, 40.0, 30.0, 0.5, mode="polar", res=1000))


def test_planet_move_with_an_unknown_mode_keeps_the_time(in_tmp):
    """A mode other than 'cumulate' / 'specific' leaves the orbital time (and so the position) unchanged."""
    t = new_target(in_tmp, 32)
    t.add_planet([0.2, 0.1, 40.0, 30.0, 0.25], 1e-4, mode="moving")
    before = np.array(t.planets[0])
    t.planet_move(0.5, mode="jump", message=False)
    assert t.orbits[0][4] == 0.25 and np.allclose(t.planets[0], before)


def test_planet_move_refuses_static_planets_and_bad_orders(in_tmp, capsys):
    """A static planet cannot move; an invalid order is reported; nothing changes."""
    t = new_target(in_tmp, 32)
    t.add_planet([0.2, 30.0], 1e-4, mode="polar")
    before = t.planets[0].copy()
    capsys.readouterr()
    t.planet_move(0.5)
    assert "A static planet is unable to move" in capsys.readouterr().out
    t.planet_move(0.5, order=2)
    assert "Input order number is not existed" in capsys.readouterr().out
    assert np.array_equal(t.planets[0], before)


def test_planet_move_with_plot_pos_saves_the_orbit_figure(in_tmp, capsys):
    """plot_pos=True saves oribit_planet<order>.png (sic) instead of printing the move message."""
    t = new_target(in_tmp, 32)
    t.add_planet([0.2, 30.0], 1e-4, mode="polar")
    t.add_planet([0.2, 0.1, 40.0, 30.0, 0.0], 1e-4, mode="moving")
    capsys.readouterr()
    t.planet_move(0.25, order=2, plot_pos=True, plot_dpi=DPI, res=50)
    assert (in_tmp/"oribit_planet2.png").read_bytes().startswith(b"\x89PNG")
    assert "Planet has moved" not in capsys.readouterr().out and t.orbits[1][4] == 0.25


def test_plot_orbit_of_a_static_planet_prints_a_message(in_tmp, capsys):
    """Regression: plot_orbit on a static planet says it has no orbit (no crash, no figure)."""
    t = new_target(in_tmp, 32)
    t.add_planet([0.2, 30.0], 1e-4, mode="polar")
    capsys.readouterr()
    t.plot_orbit(1, plot_dpi=DPI)
    assert "A static planet has no orbit" in capsys.readouterr().out
    assert not (in_tmp/"oribit_planet1.png").exists()
    t.plot_orbit(2, plot_dpi=DPI)
    assert "Input order number is not existed" in capsys.readouterr().out


def test_plot_orbit_of_a_moving_planet_saves_the_figure(in_tmp):
    """plot_orbit on a moving planet saves oribit_planet<order>.png (sic)."""
    t = new_target(in_tmp, 32)
    t.add_planet([0.2, 30.0], 1e-4, mode="polar")
    t.add_planet([0.2, 0.1, 40.0, 30.0, 0.25], 3e-4, mode="moving")
    t.plot_orbit(2, plot_dpi=DPI, res=50)
    assert (in_tmp/"oribit_planet2.png").read_bytes().startswith(b"\x89PNG")


# ---------------------------------------------------------------- method choice and PSF caches

@sizes
def test_auto_method_follows_the_target_symmetry(n, in_tmp, capsys):
    """'auto': circular for the symmetric ring, spatial (rotate) for the arc; the choice is printed."""
    t = new_target(in_tmp, n, "ring")
    capsys.readouterr()
    assert t.check_symmetry() == ("circular", None)
    assert "Method: circular (auto; target symmetric" in capsys.readouterr().out
    t = new_target(in_tmp, n, "arc")
    capsys.readouterr()
    assert t.choose_method() == ("spatial", "rotate")
    assert "Method: spatial, rotate mode (auto; target not symmetric" in capsys.readouterr().out


def test_user_symmetry_flags_override_the_detection(in_tmp):
    """[target] symmetric = true makes the arc circular; [coronagraph] symmetric = false gives spatial direct."""
    t = new_target(in_tmp, 32, "arc", target="symmetric = true")
    assert t.choose_method(verbose=False) == ("circular", None)
    t = new_target(in_tmp, 32, "ring", sections="[coronagraph]\nsymmetric = false")
    assert t.choose_method(verbose=False) == ("spatial", "direct")


def test_method_argument_and_parameters_choose_the_method(in_tmp):
    """[method] method and the method argument set the method; an unknown one raises ValueError."""
    t = new_target(in_tmp, 32, "ring", sections='[method]\nmethod = "spatial"')
    assert t.choose_method(verbose=False) == ("spatial", "rotate")
    assert t.choose_method(method="circular", verbose=False) == ("circular", None)
    assert t.choose_method(spatial_mode="direct", verbose=False) == ("spatial", "direct")
    with pytest.raises(ValueError, match="method must be"):
        t.choose_method(method="fast", verbose=False)


@sizes
def test_radial_psfs_are_cached_under_a_name_with_the_settings(n, in_tmp, capsys, monkeypatch):
    """The radial PSF file records charge and size, holds ((N+1)//2, S, S) float32, and is reused from disk."""
    t = new_target(in_tmp, n, sections="[coronagraph]\ncharge = 4")
    path = t.radial_psfs_file()
    s = extended_frame_size(n)
    assert "_c4_N%d_" % n in os.path.basename(path) and os.path.basename(path).endswith("_frame%d.npy" % s)
    assert os.path.dirname(os.path.abspath(path)) == str(in_tmp)
    frames = t.radial_psfs()
    assert frames.shape == ((n+1)//2, s, s) and frames.dtype == np.float32
    monkeypatch.setattr(Coronagraph, "_instances", {})  # a new session: only the file remains
    capsys.readouterr()
    assert new_target(in_tmp, n, sections="[coronagraph]\ncharge = 4").radial_psfs_file() == path
    assert "Calculating" not in capsys.readouterr().out


def test_invalid_charge_falls_back_to_the_parameters(in_tmp, capsys):
    """A charge that is not positive and even is replaced by [coronagraph] charge, with a message."""
    t = new_target(in_tmp, 32)
    capsys.readouterr()
    path = t.radial_psfs_file(charge=3)
    assert "use 2 from the parameters" in capsys.readouterr().out and "_c2_" in os.path.basename(path)


def test_plot_final_with_an_invalid_charge_is_named_after_the_charge_used(in_tmp, capsys):
    """plot_final(charge=3) uses [coronagraph] charge (2) and names the figure charge2."""
    t = new_target(in_tmp, 32)
    capsys.readouterr()
    t.plot_final(charge=3, method="circular", plot_dpi=DPI)
    assert "use 2 from the parameters" in capsys.readouterr().out
    assert (in_tmp/"ring_charge2_final.png").exists() and not list(in_tmp.glob("*charge3*"))


def test_method_choice_is_printed_again_only_when_it_changes(in_tmp, capsys):
    """After check_symmetry, plot_final does not repeat the same 'Method:' lines; a different method prints its own."""
    t = new_target(in_tmp, 32)
    t.check_symmetry()
    capsys.readouterr()
    t.plot_final(plot_dpi=DPI)
    assert "Method:" not in capsys.readouterr().out
    t.plot_final(method="spatial", plot_dpi=DPI)
    assert "Method: spatial, rotate mode" in capsys.readouterr().out
    t.plot_final(method="spatial", plot_dpi=DPI)
    assert "Method:" not in capsys.readouterr().out


def test_spatial_error_estimate_is_printed_once_per_image(in_tmp, capsys, monkeypatch):
    """The spatial disk image's estimated error is printed the first time only, per (charge, mode, iwa_ignore)."""
    monkeypatch.setattr(Coronagraph, "_instances", {})  # no final image made by other tests
    t = new_target(in_tmp, 32, "arc")
    capsys.readouterr()
    t.plot_final(method="spatial", add_planet=False, plot_dpi=DPI)
    assert "Estimated error of this image" in capsys.readouterr().out
    t.plot_final(method="spatial", add_planet=False, plot_dpi=DPI)
    assert "Estimated error" not in capsys.readouterr().out
    t.plot_final(method="spatial", add_planet=False, iwa_ignore=True, plot_dpi=DPI)
    assert "Estimated error of this image" in capsys.readouterr().out


# ---------------------------------------------------------------- final images

@pytest.mark.parametrize("method, add_planet, iwa_ignore, name", [
    ("circular", False, False, "ring_charge2_final.png"),
    ("circular", True, True, "ring_charge2_with_planets_iwa_ignore_final.png"),
    ("spatial", True, False, "ring_charge2_with_planets_spatial_final.png"),
])
def test_plot_final_returns_the_image_and_saves_a_named_png(method, add_planet, iwa_ignore, name, in_tmp):
    """plot_final returns the (N, N) float64 final image and saves <name>_charge<c>[...]_final.png."""
    t = new_target(in_tmp, 33)
    t.add_planet([planet_r(33), 35.0], 1e-4, mode="polar")
    img = t.plot_final(method=method, add_planet=add_planet, iwa_ignore=iwa_ignore, plot_dpi=DPI)
    assert img.shape == (33, 33) and img.dtype == np.float64 and np.all(np.isfinite(img))
    assert (in_tmp/name).read_bytes().startswith(b"\x89PNG")


def test_plot_final_file_name_saves_exactly_there(in_tmp, capsys):
    """file_name saves exactly there (folders created); the format follows the extension, PNG without one."""
    t = new_target(in_tmp, 32)
    t.plot_final(plot_dpi=DPI, file_name=os.path.join("out", "sub", "final"))
    assert "Saved the image to %s" % os.path.join("out", "sub", "final") in capsys.readouterr().out
    assert (in_tmp/"out"/"sub"/"final").read_bytes().startswith(b"\x89PNG")
    t.plot_final(plot_dpi=DPI, file_name="final.pdf")
    assert (in_tmp/"final.pdf").read_bytes().startswith(b"%PDF")
    assert not (in_tmp/"ring_charge2_final.png").exists()


def test_plot_final_shows_the_returned_image_with_plus_y_up(in_tmp):
    """The figure shows the returned array, row 0 at the bottom (+y up) over the pixel edges; flip=False turns y over."""
    t = new_target(in_tmp, 32)
    img = t.plot_final(method="circular", plot_dpi=DPI)
    ax = plt.gcf().axes[0]
    shown = ax.images[0]
    assert np.array_equal(shown.get_array(), img) and shown.origin == "lower"
    assert np.allclose(shown.get_extent(), t.extent) and not ax.yaxis_inverted()
    plt.close("all")
    t.plot_final(method="circular", plot_dpi=DPI, flip=False)
    assert plt.gcf().axes[0].yaxis_inverted()


def test_plot_origin_shows_the_target_with_plus_y_up(in_tmp):
    """plot_origin shows pre_img with row 0 at the bottom, axes limited to the image; flip=False turns y over."""
    t = new_target(in_tmp, 33)
    t.plot_origin(plot_dpi=DPI)
    ax = plt.gcf().axes[0]
    assert np.array_equal(ax.images[0].get_array(), t.pre_img) and ax.images[0].origin == "lower"
    assert np.allclose(ax.get_xlim(), t.extent[:2]) and np.allclose(ax.get_ylim(), t.extent[2:])
    plt.close("all")
    t.plot_origin(plot_dpi=DPI, flip=False)
    assert plt.gcf().axes[0].yaxis_inverted()


def test_plot_origin_and_final_log_scale(in_tmp):
    """log=False (default): a linear scale from 0 to the peak; log=True: log_range decades below the peak, the
    pixels without light in the darkest colour (not blank)."""
    t = new_target(in_tmp, 33)
    t.plot_origin(plot_dpi=DPI)
    norm = plt.gcf().axes[0].images[0].norm
    assert type(norm).__name__ == "Normalize" and (norm.vmin, norm.vmax) == (0, t.pre_img.max())
    plt.close("all")
    img = t.plot_final(method="circular", plot_dpi=DPI, log=True, log_range=3)
    shown = plt.gcf().axes[0].images[0]
    assert type(shown.norm).__name__ == "LogNorm"
    assert shown.norm.vmax == pytest.approx(img.max()) and shown.norm.vmin == pytest.approx(img.max()*1e-3)
    assert shown.cmap.get_bad()[3] == 1 and np.allclose(shown.cmap.get_bad(), shown.cmap(0.0))
    plt.close("all")


@pytest.mark.parametrize("unit_line, label", [("", "Jy/arcsec^2"), ('unit = "mJy/arcsec^2"', "mJy/arcsec^2")])
def test_colour_bar_label_follows_the_target_unit(unit_line, label, in_tmp):
    """The colour bars of plot_origin and plot_final show [target] unit ("none": Jy/arcsec^2)."""
    t = new_target(in_tmp, 32, target=unit_line)
    t.plot_origin(plot_dpi=DPI)
    assert plt.gcf().axes[1].get_ylabel() == label
    plt.close("all")
    t.plot_final(method="circular", plot_dpi=DPI)
    assert plt.gcf().axes[1].get_ylabel() == label


def test_plot_origin_draws_a_planet_outside_the_image_only_without_boundary(in_tmp, capsys):
    """A planet beyond the image edge is left out with boundary=True and drawn grey with boundary=False."""
    t = new_target(in_tmp, 32)
    t.add_planet([1.0, 0.0], 1e-4, mode="cartesian")

    def circles():
        """The planet circles on the current figure's image axes."""
        return [a for a in plt.gcf().axes[0].get_children() if isinstance(a, Circle)]

    t.plot_origin(plot_dpi=DPI)
    assert circles() == []
    plt.close("all")
    t.plot_origin(plot_dpi=DPI, boundary=False)
    (circle,) = circles()
    assert np.allclose(circle.center, (1.0, 0.0)) and np.allclose(circle.get_facecolor(), to_rgba("#7f7f7f"))


@sizes
def test_both_methods_agree_on_a_symmetric_target(n, in_tmp):
    """On the symmetric ring the spatial final image matches the circular one within 3% (relative L1)."""
    t = new_target(in_tmp, n)
    cir = t.plot_final(method="circular", plot_dpi=DPI)
    sp = t.plot_final(method="spatial", plot_dpi=DPI)
    assert np.abs(sp-cir).sum()/cir.sum() < 0.03


@sizes
@pytest.mark.parametrize("method", ["circular", "spatial"])
def test_planet_orientation_matches_grid_pixel(n, method, in_tmp):
    """Regression: a planet at (r, angle) appears at scene.Grid.pixel in both methods (+y = increasing row)."""
    t = new_target(in_tmp, n)
    r, angle = planet_r(n), 35
    t.add_planet([r, angle], 1e-4, mode="polar")
    expected = np.array(Grid(size=n, pixel_scale=PIXSCALE, wavelength=WAVELEN, diameter=DIAMETER).pixel(r, np.radians(angle)))
    assert expected[0] > (n-1)/2  # above the star
    with_p = t.plot_final(method=method, plot_dpi=DPI)
    without = t.plot_final(method=method, add_planet=False, plot_dpi=DPI)
    assert np.abs(peak(with_p-without)-expected).max() <= 1


@pytest.mark.parametrize("method", ["circular", "spatial"])
def test_contrast_returns_planet_and_dust_light(method, in_tmp, capsys):
    """contrast gives planet > 0 and 0 <= dust (IWA ignored) <= dust, in Jy; an invalid order gives None."""
    t = new_target(in_tmp, 32)
    t.add_planet([planet_r(32), 35.0], 1e-4, mode="polar")
    planet_b, dust_b, dust_iwa = t.contrast(plot_contrast=False, method=method)
    assert "This is a charge-2 vortex coronagraph." in capsys.readouterr().out
    assert planet_b > 0 and 0 <= dust_iwa <= dust_b
    assert t.contrast(order=2, plot_contrast=False, method=method) is None
    assert "Input order number is not existed" in capsys.readouterr().out


def test_contrast_plot_shows_both_dust_to_planet_ratios_without_saving(in_tmp):
    """plot_contrast draws dust/planet * dust factor, all the dust and the IWA's left out, over [0, exozodi_limit],
    and a reference line at 1; no file."""
    t = new_target(in_tmp, 32)
    t.add_planet([planet_r(32), 35.0], 1e-4, mode="polar")
    planet_b, dust_b, dust_iwa = t.contrast(method="circular", plot_dpi=DPI, exozodi_limit=5, res=11)
    complete, ignored, one = plt.gcf().axes[0].lines
    assert complete.get_label() == "all the dust" and ignored.get_label() == "dust inside the IWA left out"
    assert np.allclose(one.get_ydata(), 1)
    assert np.allclose(complete.get_xdata(), np.linspace(0, 5, 11))
    assert np.allclose(complete.get_ydata(), dust_b/planet_b*np.linspace(0, 5, 11))
    assert np.allclose(ignored.get_ydata(), dust_iwa/planet_b*np.linspace(0, 5, 11))
    assert not list(in_tmp.glob("*.png")) and not list(in_tmp.glob("*.pdf"))


def test_contrast_circular_planet_beyond_the_radial_psfs_gives_none(in_tmp, capsys):
    """Circular method: a planet beyond the radial PSFs (outside the image) has no contrast (None, message)."""
    t = new_target(in_tmp, 32)
    t.add_planet([0.6*32*PIXSCALE, 0.0], 1e-4, mode="polar")
    assert t.contrast(plot_contrast=False, method="circular") is None
    assert "Planet is outside the range of the plot" in capsys.readouterr().out


@sizes
def test_contrast_circular_planet_just_beyond_the_last_radial_psf_gives_none(n, in_tmp, capsys):
    """Circular method: a planet beyond the last radial PSF, which plot_final leaves out, has no contrast (None)."""
    t = new_target(in_tmp, n)
    r_px = (radial_source_positions(n)[-1]+(n+1)//2)/2  # between the last radial PSF and (N+1)//2
    t.add_planet([r_px*PIXSCALE, 45.0], 1e-4, mode="polar")
    with_p = t.plot_final(method="circular", plot_dpi=DPI)
    assert np.array_equal(with_p, t.plot_final(method="circular", add_planet=False, plot_dpi=DPI))
    capsys.readouterr()
    assert t.contrast(plot_contrast=False, method="circular") is None
    assert "Planet is outside the range of the plot" in capsys.readouterr().out


def test_plot_origin_saves_named_figures_and_reports_planets_outside(in_tmp, capsys):
    """plot_origin saves <name>_origin.png, or <name>_origin_with_planets.png and names planets outside."""
    t = new_target(in_tmp, 32)
    t.plot_origin(plot_dpi=DPI)
    assert (in_tmp/"ring_origin.png").exists()
    t.add_planet([1.0, 0.0], 1e-4, mode="cartesian")  # beyond the image edge (0.275")
    capsys.readouterr()
    t.plot_origin(plot_dpi=DPI)
    assert "Planet #1 is outside the range of the plot" in capsys.readouterr().out
    assert (in_tmp/"ring_origin_with_planets.png").exists()


# ---------------------------------------------------------------- planet video

def _ffprobe(path):
    """'codec,width,height,pix_fmt' of the first video stream (ffprobe)."""
    return subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0", "-show_entries",
                           "stream=codec_name,width,height,pix_fmt", "-of", "csv=p=0", str(path)],
                          capture_output=True, text=True, check=True).stdout.strip()


@pytest.mark.skipif(shutil.which("ffprobe") is None or shutil.which("ffmpeg") is None, reason="needs ffmpeg/ffprobe")
@pytest.mark.parametrize("method, iwa_ignore, name", [
    ("circular", False, "arc_planet_video.webm"),
    ("spatial", True, "arc_planet_video_iwa_ignore_spatial.webm"),
])
def test_planet_video_is_lossless_vp9_webm_and_leaves_the_planet_in_place(method, iwa_ignore, name, in_tmp):
    """The video is VP9 on RGB (gbrp) in WebM at 10 x 8 inch * plot_dpi; the planet ends where it started."""
    t = new_target(in_tmp, 32, "arc")
    t.add_planet([planet_r(32), 0.2, 40.0, 30.0, 0.1], 5e-5, mode="moving")
    before, t_before = np.array(t.planets[0]), t.orbits[0][4]
    t.planet_video(order=1, length=1, fps=2, plot_dpi=10, method=method, iwa_ignore=iwa_ignore)
    assert _ffprobe(in_tmp/name) == "vp9,100,80,gbrp"
    assert np.allclose(t.planets[0], before) and t.orbits[0][4] == t_before


@pytest.mark.skipif(shutil.which("ffmpeg") is None, reason="needs ffmpeg")
def test_planet_video_frames_decode_to_exactly_the_rendered_pixels(in_tmp, monkeypatch):
    """Lossless: the video has length*fps frames and each decodes to the RGB pixels matplotlib rendered for it."""
    rendered = []
    grab = animation.FFMpegWriter.grab_frame

    def grab_and_keep(writer, **kwargs):
        """FFMpegWriter.grab_frame, also keeping the RGBA bytes of the frame sent to ffmpeg."""
        grab(writer, **kwargs)
        buf = io.BytesIO()
        writer.fig.savefig(buf, format=writer.frame_format, dpi=writer.dpi, **kwargs)
        rendered.append(np.frombuffer(buf.getvalue(), np.uint8).reshape(80, 100, 4))

    monkeypatch.setattr(animation.FFMpegWriter, "grab_frame", grab_and_keep)
    t = new_target(in_tmp, 32, "arc")
    t.add_planet([planet_r(32), 0.2, 40.0, 30.0, 0.1], 5e-5, mode="moving")
    t.planet_video(length=1, fps=3, plot_dpi=10, method="circular")
    raw = subprocess.run(["ffmpeg", "-v", "error", "-i", "arc_planet_video.webm", "-f", "rawvideo", "-pix_fmt", "rgb24",
                          "-"], capture_output=True, check=True).stdout
    decoded = np.frombuffer(raw, np.uint8).reshape(-1, 80, 100, 3)
    assert len(rendered) == len(decoded) == 3
    assert all(f[..., 3].min() == 255 for f in rendered)  # opaque frames: RGB carries everything
    assert np.array_equal(decoded, np.array(rendered)[..., :3])
    assert not np.array_equal(decoded[0], decoded[1])  # the planet moves


def test_planet_video_of_a_static_planet_is_refused(in_tmp, capsys):
    """A static planet has no orbit to follow: a message, and no video file."""
    t = new_target(in_tmp, 32)
    t.add_planet([planet_r(32), 35.0], 1e-4, mode="polar")
    capsys.readouterr()
    t.planet_video(length=1, fps=2, plot_dpi=10)
    assert "A static planet is unable to move" in capsys.readouterr().out
    t.planet_video(order=2, length=1, fps=2, plot_dpi=10)
    assert "Input order number is not existed" in capsys.readouterr().out
    assert not list(in_tmp.glob("*.webm"))
