"""
Stage 5 unit tests: the older-layout parameter file (toycoronagraph/config.py: from_dict defaults and
checks, load_config incl. the legacy .py file, resolve against the FITS header, used_path / write_used /
write_resolved), the input methods (toycoronagraph/inputs.py: plain images, cubes, HDUs, MCFOST cubes,
unit conversion) and the method choice (toycoronagraph/selection.py: image_asymmetry, image_symmetry,
psf_symmetry, choose_method).

Images are N = 32 (even) and 33 (odd) pixels where the centre convention matters; files are written
only in each test's temporary folder (in_tmp). The longer checks stay in tests/stage5_selection.py.
"""
import os
import tomllib
from types import SimpleNamespace
import numpy as np
import pytest
from astropy.io import fits
from toycoronagraph import DATADIR
from toycoronagraph import config as C
from toycoronagraph.config import (SCHEMA, RESOLVED_MARK, from_dict, load_config, resolve, used_path, write_used,
                                   write_resolved, write_example, example_text)
from toycoronagraph.inputs import read_plane, read_target, target_path, is_mcfost
from toycoronagraph.optics import Optics
from toycoronagraph.selection import image_asymmetry, image_symmetry, psf_symmetry, choose_method

RAD = 206264.806247
SIZES = [32, 33]
sizes = pytest.mark.parametrize("n", SIZES, ids=["even", "odd"])

# The old-style parameter file, as toycoronagraph.para used to write it (small image)
LEGACY_PY = '''
fits_filename = None
target_pixel = 32
num_cores = 2
rot_number = target_pixel
wavelength = 1.6e-6
D = 6.5
rad_to_arcsecond = 206264.806247
visual_range = 24
psf_scale = wavelength/D*rad_to_arcsecond*visual_range/target_pixel
light_speed = 299792458
jy=10**26
F_transfer = jy*wavelength/light_speed/psf_scale**2
coronagraph_type='vortex'
lyot_mask_size = 0.9
psf_range=visual_range/2
'''
MCFOST_PLANES = ["I", "Q", "U", "V", "direct star", "scattered star", "direct thermal", "scattered thermal"]


def write_fits(path, data, **header):
    """Write data to a FITS file with the given header keywords."""
    hdr = fits.Header()
    for key, value in header.items():
        hdr[key] = value
    fits.PrimaryHDU(data, hdr).writeto(path, overwrite=True)


def mcfost_cube(path, n):
    """Write a random MCFOST-like cube (8 planes, 2 azimuths, 3 inclinations, n x n) with its header; returns it."""
    cube = np.random.default_rng(0).random((8, 2, 3, n, n)).astype(np.float32)
    header = {"BUNIT": "W.m-2.pixel-1", "CDELT1": -0.01/3600, "CDELT2": 0.01/3600, "CUNIT2": "deg", "WAVE": 0.8}
    header.update({"FLUX_%d" % (i+1): name for i, name in enumerate(MCFOST_PLANES)})
    write_fits(path, cube, **header)
    return cube


def cfg_with(**target):
    """A configuration with the given [target] keys, everything else default."""
    return from_dict({"target": target} if target else {})


def targets(n):
    """Symmetric and non-symmetric test targets of size n (dicts name -> (n, n) image)."""
    c = (n-1)/2
    rows, cols = np.indices((n, n))
    rr = np.hypot(rows-c, cols-c)
    ring = lambda r: np.exp(-0.5*((r-0.3*n)/(0.06*n))**2)
    symmetric = {"smooth ring": ring(rr),
                 "sharp annulus": ((rr > 0.2*n) & (rr < 0.35*n)).astype(float),
                 "steep profile": np.exp(-rr/(0.05*n))}
    asymmetric = {"inclined": ring(np.hypot(cols-c, (rows-c)/0.95)),
                  "off-centre by 1 px": ring(np.hypot(rows-c-1, cols-c)),
                  "8-fold pattern": ring(rr)*(1+0.3*np.cos(8*np.arctan2(rows-c, cols-c))),
                  "half ring": ring(rr)*(cols > c),
                  "light in the corners": ring(rr)+3*np.exp(-0.5*(np.hypot(rows, cols)/(0.05*n))**2)}
    return symmetric, asymmetric


# ---------------------------------------------------------------- config: from_dict

def test_from_dict_empty_gives_schema_defaults():
    """An empty file gives every SCHEMA default, with the documented renames and derived keys."""
    cfg = from_dict({})
    rename = {("target", "symmetric"): "image_symmetric", ("coronagraph", "symmetric"): "psf_symmetric",
              ("coronagraph", "type"): "coronagraph_type"}
    for section, keys in SCHEMA.items():
        for key, (default, _, _) in keys.items():
            if (section, key) != ("target", "fits_filename"):
                assert getattr(cfg, rename.get((section, key), key)) == default, (section, key)
    assert cfg.fits_filename is None  # "" means the example target
    assert cfg.D == cfg.diameter == 2.4
    assert cfg.legacy is False and cfg.source == "parameters"
    assert cfg.target_pixel == cfg.wavelength == cfg.visual_range == cfg.input_unit == "auto"


def test_from_dict_keeps_given_values_and_converts_numbers():
    """Given values are kept; integers where a float is expected become floats."""
    cfg = from_dict({"method": {"method": "spatial", "error_target": 1}, "target": {"symmetric": True, "target_pixel": 64},
                     "telescope": {"wavelength": 2, "diameter": 6}, "coronagraph": {"charge": 6, "symmetric": False}})
    assert cfg.method == "spatial" and cfg.spatial_mode == "auto"
    assert cfg.error_target == 1.0 and isinstance(cfg.error_target, float)
    assert cfg.wavelength == 2.0 and isinstance(cfg.wavelength, float)
    assert cfg.D == 6.0 and isinstance(cfg.diameter, float)
    assert cfg.target_pixel == 64 and isinstance(cfg.target_pixel, int)
    assert cfg.image_symmetric is True and cfg.psf_symmetric is False and cfg.charge == 6


def test_from_dict_ignores_a_resolved_section():
    """A [resolved] section (values of an earlier run) is not read back."""
    cfg = from_dict({"resolved": {"psf_scale": 0.5, "anything": "x"}, "telescope": {"visual_range": "auto"}})
    assert cfg.visual_range == "auto" and not hasattr(cfg, "anything")


@pytest.mark.parametrize("data, words", [
    ({"methods": {}}, "unknown section"),
    ({"method": {"methd": "auto"}}, "unknown key"),
    ({"target": 5}, "must be a table"),
    ({"method": {"method": "fast"}}, "must be one of"),
    ({"target": {"symmetric": 1}}, "must be one of"),
    ({"target": {"target_pixel": "512"}}, 'must be "auto" or an int'),
    ({"target": {"target_pixel": True}}, 'must be "auto" or an int'),
    ({"compute": {"use_gpu": 1}}, "must be a bool"),
    ({"compute": {"num_cores": 2.0}}, r"\[compute\] num_cores"),
    ({"target": {"fits_slice": [5.0]}}, "list of integers"),
    ({"target": {"input_unit": "mJy"}}, "must be one of"),
    ({"target": {"mcfost_component": "dust"}}, "must be one of"),
    ({"telescope": {"wavelength": "1um"}}, 'must be "auto" or a float'),
    ({"telescope": {"diameter": True}}, "must be a float"),
    ({"coronagraph": {"charge": 3}}, "positive and even"),
    ({"coronagraph": {"charge": 0}}, "positive and even"),
    ({"coronagraph": {"type": "lyot"}}, "must be one of"),
], ids=lambda v: str(v) if isinstance(v, str) else None)
def test_from_dict_rejects_bad_files(data, words):
    """Unknown sections/keys and values outside SCHEMA raise ValueError saying what is wrong."""
    with pytest.raises(ValueError, match=words):
        from_dict(data)


def test_from_dict_errors_name_the_file():
    """Error messages about unknown sections start with the file name."""
    with pytest.raises(ValueError, match=r"^my\.toml: unknown section \[methods\]"):
        from_dict({"methods": {}}, source="my.toml")


def test_toml_value_formats_values_as_toml():
    """Lists become arrays, booleans true/false, strings are quoted, numbers use repr."""
    assert C._toml_value([5, 0, 0]) == "[5, 0, 0]"
    assert C._toml_value((True, False)) == "[true, false]"
    assert C._toml_value("auto") == '"auto"'
    assert C._toml_value(1.6e-6) == "1.6e-06" and C._toml_value(32) == "32"


# ---------------------------------------------------------------- config: files

def test_example_file_lists_every_key_and_loads_back_to_defaults(in_tmp):
    """write_example writes every SCHEMA section and key with the usage note; it loads back to the defaults."""
    path = str(in_tmp/"toycoronagraph_para.toml")
    write_example(path)
    cfg, defaults = load_config(path), from_dict({})
    for key in vars(defaults):
        if key != "source":
            assert getattr(cfg, key) == getattr(defaults, key), key
    assert cfg.source == path
    text = example_text()
    assert "Using your own target" in text and "MCFOST" in text
    for section, keys in SCHEMA.items():
        assert "[%s]" % section in text and all(key in text for key in keys)


def test_load_config_without_any_file_creates_a_minimal_one(in_tmp, capsys):
    """With neither file, a minimal toycoronagraph_para.toml is created and read (example target, "auto" size)."""
    cfg = load_config()
    assert (in_tmp/"toycoronagraph_para.toml").exists()
    assert "created a minimal one" in capsys.readouterr().out
    assert cfg.target_pixel == "auto" and cfg.fits_filename is None and cfg.charge == 2


def test_load_config_reads_a_legacy_py_file_with_a_notice(in_tmp, capsys):
    """An old toycoronagraph_para.py is read when there is no .toml; its derived values are kept by resolve."""
    (in_tmp/"toycoronagraph_para.py").write_text(LEGACY_PY)
    cfg = load_config()
    assert "old-style toycoronagraph_para.py" in capsys.readouterr().out
    assert cfg.legacy is True and cfg.source == "toycoronagraph_para.py"
    assert cfg.target_pixel == 32 and cfg.num_cores == 2 and cfg.lyot_mask_size == 0.9 and cfg.D == 6.5
    resolve(cfg, (32, 32), {"BUNIT": "W/m^2/pixel", "PIXSCALE": 0.5, "WAVELEN": 1e-6})
    assert np.isclose(cfg.psf_scale, 1.6e-6/6.5*RAD*24/32) and cfg.wavelength == 1.6e-6  # the header is not used
    assert cfg.psf_range == 12 and cfg.rot_number == 32 and cfg.method == "auto"
    assert not any("WARNING" in note for note in cfg.notes)


def test_legacy_py_file_missing_values_take_the_old_defaults(in_tmp):
    """A .py file with only target_pixel: 1 um, D 2.4 m, 32 lambda/D; psf_scale, F_transfer, psf_range derived."""
    (in_tmp/"old.py").write_text('target_pixel = 64\nfits_filename = "disk"\nunit = "mJy/arcsec^2"\n')
    cfg = load_config("old.py")
    lambda_d = 1e-6/2.4*RAD
    assert cfg.legacy is True and cfg.source == "old.py"
    assert cfg.wavelength == 1e-6 and cfg.D == 2.4 and cfg.visual_range == 32 and cfg.lyot_mask_size == 0.95
    assert np.isclose(cfg.psf_scale, lambda_d*32/64) and cfg.psf_range == 16
    assert np.isclose(cfg.F_transfer, 1e26*1e-6/299792458/cfg.psf_scale**2)
    assert cfg.num_cores == -1 and cfg.rot_number == 0 and cfg.charge == 2 and cfg.method == "auto"
    assert cfg.fits_filename == "disk" and cfg.unit == "mJy/arcsec^2"


def test_load_config_legacy_notice_only_when_verbose(in_tmp, capsys):
    """verbose=False reads the old .py file silently; an explicit .py path is read without the notice."""
    (in_tmp/"toycoronagraph_para.py").write_text(LEGACY_PY)
    assert load_config(verbose=False).legacy is True
    assert load_config("toycoronagraph_para.py").legacy is True
    assert capsys.readouterr().out == ""


def test_load_config_prefers_the_toml_file_over_the_legacy_one(in_tmp):
    """When both exist, toycoronagraph_para.toml is read and the .py file ignored."""
    (in_tmp/"toycoronagraph_para.py").write_text(LEGACY_PY)
    (in_tmp/"toycoronagraph_para.toml").write_text('[target]\ntarget_pixel = 128\n')
    cfg = load_config()
    assert cfg.target_pixel == 128 and cfg.legacy is False


@pytest.mark.parametrize("source, expected", [
    ("p.toml", "p_used.toml"),
    (os.path.join("sub", "p.toml"), os.path.join("sub", "p_used.toml")),
    ("toycoronagraph_para.py", "toycoronagraph_para_used.toml"),
    ("parameters", "toycoronagraph_para_used.toml"),
])
def test_used_path_is_next_to_the_parameter_file(source, expected):
    """The values used go to <parameter file name>_used.toml (a default name when there is no file)."""
    assert used_path(SimpleNamespace(source=source)) == expected


def _resolved_file_cfg(folder, n=16):
    """Write folder/d.fits (Jy/pixel, PIXSCALE 0.01, WAVELEN 1.6e-6) and folder/p.toml; load and read them."""
    write_fits(str(folder/"d.fits"), np.ones((n, n)), BUNIT="Jy/pixel", PIXSCALE=0.01, WAVELEN=1.6e-6)
    user = '# my settings\n[target]\nfits_filename = "d"  # the disk\n[telescope]\ndiameter = 6.5\n'
    (folder/"p.toml").write_text(user)
    cfg = load_config("p.toml")
    read_target(cfg, verbose=False)
    return cfg, user


def test_write_used_saves_the_values_used_and_leaves_the_file_alone(in_tmp):
    """write_used writes p_used.toml (valid TOML, [resolved] = the values used); p.toml is unchanged."""
    cfg, user = _resolved_file_cfg(in_tmp)
    path = write_used(cfg)
    assert path == "p_used.toml" and (in_tmp/"p.toml").read_text() == user
    text = (in_tmp/path).read_text()
    assert text.startswith("# Values used by toycoronagraph for p.toml") and text.count(RESOLVED_MARK) == 1
    used = tomllib.loads(text)["resolved"]
    assert used["target_pixel"] == 16 and used["input_unit"] == "Jy/pixel" and used["fits_slice"] == []
    assert used["psf_scale"] == cfg.psf_scale == 0.01 and used["wavelength"] == 1.6e-6
    assert np.isclose(used["visual_range"], cfg.visual_range) and used["rot_number"] == 16
    assert "# header PIXSCALE" in text and used["input"] == "d.fits (image)"


def test_write_used_twice_keeps_one_resolved_section(in_tmp):
    """Each run rewrites the _used file: one [resolved] section, with the latest values."""
    cfg, _ = _resolved_file_cfg(in_tmp)
    write_used(cfg)
    cfg.psf_scale = 0.02
    text = (in_tmp/write_used(cfg)).read_text()
    assert text.count("[resolved]") == 1 and tomllib.loads(text)["resolved"]["psf_scale"] == 0.02


def test_write_back_appends_the_section_once_and_it_is_not_read_back(in_tmp):
    """With write_back the [resolved] section replaces the previous one at the end of p.toml; "auto" stays "auto"."""
    cfg, user = _resolved_file_cfg(in_tmp)
    cfg.write_back = True
    write_used(cfg)
    write_used(cfg)
    text = (in_tmp/"p.toml").read_text()
    assert text.startswith(user.rstrip("\n")) and text.count(RESOLVED_MARK) == 1 and text.count("[resolved]") == 1
    again = load_config("p.toml")
    assert again.visual_range == "auto" and again.target_pixel == "auto" and again.diameter == 6.5


def test_write_used_for_a_legacy_py_file_writes_the_default_used_file(in_tmp):
    """A legacy .py configuration saves its values to toycoronagraph_para_used.toml; the .py file is unchanged."""
    (in_tmp/"toycoronagraph_para.py").write_text(LEGACY_PY)
    cfg = load_config(verbose=False)
    resolve(cfg, (32, 32), {"BUNIT": "Jy/pixel"})
    assert write_used(cfg) == "toycoronagraph_para_used.toml"
    text = (in_tmp/"toycoronagraph_para_used.toml").read_text()
    used = tomllib.loads(text)["resolved"]
    assert used["target_pixel"] == 32 and used["wavelength"] == 1.6e-6 and np.isclose(used["psf_scale"], cfg.psf_scale)
    assert "# old .py file" in text and (in_tmp/"toycoronagraph_para.py").read_text() == LEGACY_PY


def test_write_resolved_skips_files_that_are_not_toml(in_tmp):
    """write_resolved leaves a legacy .py parameter file untouched."""
    (in_tmp/"toycoronagraph_para.py").write_text(LEGACY_PY)
    cfg = C._from_legacy("toycoronagraph_para.py")
    resolve(cfg, (32, 32), {"BUNIT": "Jy/pixel"})
    write_resolved(cfg)
    assert (in_tmp/"toycoronagraph_para.py").read_text() == LEGACY_PY


# ---------------------------------------------------------------- config: resolve

@sizes
def test_resolve_fills_auto_values_from_the_header(n):
    """'auto' takes the size, unit, pixel scale and wavelength from the header; derived values follow."""
    cfg = cfg_with()
    notes = resolve(cfg, (n, n), {"BUNIT": "Jy/arcsec^2", "PIXSCALE": 0.01, "WAVELEN": 1.6e-6})
    lambda_d = 1.6e-6/2.4*RAD
    assert cfg.target_pixel == n and cfg.input_unit == "Jy/arcsec^2" and cfg.wavelength == 1.6e-6
    assert np.isclose(cfg.psf_scale, 0.01) and np.isclose(cfg.visual_range, n*0.01/lambda_d)
    assert np.isclose(cfg.psf_range, cfg.visual_range/2) and cfg.rot_number == n
    assert np.isclose(cfg.F_transfer, 1e26*1.6e-6/299792458/0.01**2)
    assert notes is cfg.notes and notes[0].startswith("Input:") and not any("WARNING" in x for x in notes)
    assert cfg.sources["psf_scale"] == "header PIXSCALE" and cfg.sources["target_pixel"] == "image size"


def test_resolve_reads_mcfost_style_headers():
    """CDELT in degrees and WAVE in micrometers (MCFOST) give the pixel scale in arcsec and the wavelength in m."""
    cfg = cfg_with()
    resolve(cfg, (32, 32), {"BUNIT": "W.m-2.pixel-1", "CDELT1": -0.02/3600, "CDELT2": 0.01/3600, "CUNIT2": "deg",
                            "WAVE": 1.6})
    assert cfg.input_unit == "W/m^2/pixel" and np.isclose(cfg.psf_scale, 0.01) and cfg.wavelength == 1.6e-6


def test_resolve_without_header_uses_the_old_defaults():
    """No wavelength or pixel scale in the header: 1 um and 32 lambda/D across the image."""
    cfg = cfg_with(input_unit="Jy/pixel")
    resolve(cfg, (32, 32), {})
    assert cfg.wavelength == 1e-6 and cfg.visual_range == 32 and np.isclose(cfg.psf_scale, 1e-6/2.4*RAD)


def test_resolve_explicit_values_override_the_header_with_warnings():
    """Explicit wavelength, visual_range and input_unit win over the header; each mismatch is a WARNING."""
    cfg = from_dict({"telescope": {"wavelength": 1.0e-6, "visual_range": 20.0}, "target": {"input_unit": "Jy/pixel"}})
    resolve(cfg, (32, 32), {"BUNIT": "Jy/arcsec^2", "PIXSCALE": 0.01, "WAVELEN": 1.6e-6})
    assert cfg.wavelength == 1e-6 and cfg.visual_range == 20 and cfg.input_unit == "Jy/pixel"
    assert np.isclose(cfg.psf_scale, 1e-6/2.4*RAD*20/32)
    warnings = [x for x in cfg.notes if "WARNING" in x]
    assert len(warnings) == 3
    assert any("wavelength" in w for w in warnings) and any("visual_range" in w for w in warnings)
    assert any("input_unit" in w for w in warnings)


@pytest.mark.parametrize("cunit, scale, source", [
    (None, 0.01*3600, "header CDELT1 (deg)"),   # no CUNIT: degrees
    ("arcsec", 0.01, "header CDELT1 (arcsec)"),
    ("ARCSEC", 0.01, "header CDELT1 (arcsec)"),
    ("mas", 1e-5, "header CDELT1 (mas)"),
    ("rad", 0.01*RAD, "header CDELT1 (rad)"),
    ("furlong", 1e-6/2.4*RAD*32/32, "default (no pixel scale in the header)"),
])
def test_resolve_pixel_scale_from_cdelt_and_its_unit(cunit, scale, source):
    """CDELT1 with CUNIT1 in deg (default), arcsec, mas or rad gives arcsec/pixel; an unknown unit is ignored."""
    header = {"CDELT1": -0.01}
    if cunit:
        header["CUNIT1"] = cunit
    cfg = cfg_with(input_unit="Jy/pixel")
    resolve(cfg, (32, 32), header)
    assert np.isclose(cfg.psf_scale, scale) and cfg.sources["psf_scale"] == source


@pytest.mark.parametrize("factor, warned", [(1.005, False), (1.02, True)])
def test_resolve_reports_mismatches_beyond_one_percent_only(factor, warned):
    """Explicit wavelength / pixel scale within 1% of the header are not reported; beyond it they are."""
    lambda_d = 1.6e-6/2.4*RAD
    # wavelength off by `factor`; visual_range kept, so the PSF pixel scale is 0.01 * factor too
    cfg = from_dict({"telescope": {"wavelength": 1.6e-6*factor, "visual_range": 32*0.01/lambda_d},
                     "target": {"input_unit": "Jy/pixel"}})
    resolve(cfg, (32, 32), {"PIXSCALE": 0.01, "WAVELEN": 1.6e-6})
    assert np.isclose(cfg.psf_scale, 0.01*factor)
    assert sum("WARNING" in note for note in cfg.notes) == (2 if warned else 0)


def test_resolve_without_a_header_and_with_a_matching_target_pixel():
    """header None is an empty header; an explicit target_pixel equal to the image size is accepted."""
    cfg = cfg_with(input_unit="Jy/pixel", target_pixel=16)
    resolve(cfg, (16, 16), None)
    assert cfg.target_pixel == 16 and cfg.sources["target_pixel"] == "parameters (matches the image)"
    assert cfg.wavelength == 1e-6 and cfg.visual_range == 32


def test_resolve_keeps_an_explicit_rot_number():
    """rot_number 0 means target_pixel; a positive value is kept."""
    cfg = from_dict({"method": {"rot_number": 90}, "target": {"input_unit": "Jy/pixel"}})
    resolve(cfg, (32, 32), {})
    assert cfg.rot_number == 90


@pytest.mark.parametrize("shape, target, header, words", [
    ((32, 33), {}, {"BUNIT": "Jy/pixel"}, "square and 2-D"),
    ((2, 32, 32), {}, {"BUNIT": "Jy/pixel"}, "square and 2-D"),
    ((32, 32), {"target_pixel": 33}, {"BUNIT": "Jy/pixel"}, "target_pixel = 33"),
    ((32, 32), {}, {"PIXSCALE": 0.01}, "unit of the input image is unknown"),
    ((32, 32), {}, {"BUNIT": "furlongs"}, "unit of the input image is unknown"),
], ids=["not square", "3-D", "size mismatch", "no BUNIT", "unknown BUNIT"])
def test_resolve_explains_unusable_input(shape, target, header, words):
    """A non-square image, a size mismatch, or an unknown unit raise ValueError explaining the problem."""
    with pytest.raises(ValueError, match=words):
        resolve(cfg_with(**target), shape, header)


# ---------------------------------------------------------------- inputs

@pytest.mark.parametrize("name, expected", [
    (None, os.path.join(DATADIR, "example.fits")),
    ("disk", "disk.fits"),
    (os.path.join("data", "disk.fits"), os.path.join("data", "disk.fits")),
    ("RT.fits.gz", "RT.fits.gz"),
    ("image.FIT", "image.FIT"),
    ("image.fts", "image.fts"),
])
def test_target_path_adds_fits_only_without_an_extension(name, expected):
    """target_path: the example for None, '.fits' appended when the name has no FITS extension."""
    assert target_path(SimpleNamespace(fits_filename=name)) == expected


def test_is_mcfost_needs_five_axes_four_or_eight_planes_and_flux_keywords():
    """MCFOST cubes are recognized by their 5 axes, 4 or 8 planes and the FLUX_1 keyword."""
    hdr = {"FLUX_1": "I"}
    assert is_mcfost(np.zeros((8, 1, 1, 4, 4)), hdr) and is_mcfost(np.zeros((4, 1, 1, 4, 4)), hdr)
    assert not is_mcfost(np.zeros((8, 1, 1, 4, 4)), {})
    assert not is_mcfost(np.zeros((5, 1, 1, 4, 4)), hdr) and not is_mcfost(np.zeros((4, 4)), hdr)


def test_read_plane_reads_a_2d_image_as_float64(in_tmp):
    """A plain 2-D image is returned as float64, with its header, format 'image' and no index."""
    img = np.random.default_rng(1).random((5, 5)).astype(np.float32)
    write_fits("img.fits", img, BUNIT="Jy/pixel")
    data, header, fmt, what, index = read_plane("img.fits")
    assert data.dtype == np.float64 and np.array_equal(data, img.astype(np.float64))
    assert header["BUNIT"] == "Jy/pixel" and fmt == "image" and what == "image" and index == []


def test_read_plane_reshapes_a_flattened_square_image(in_tmp):
    """A 1-D array of n*n values is reshaped to (n, n)."""
    img = np.arange(36.0).reshape(6, 6)
    write_fits("flat.fits", img.ravel())
    assert np.array_equal(read_plane("flat.fits")[0], img)


def test_read_plane_cube_needs_fits_slice(in_tmp):
    """A 3-D cube is refused without fits_slice (message names fits_slice); fits_slice picks the plane."""
    cube = np.random.default_rng(2).random((3, 6, 6))
    write_fits("cube.fits", cube)
    with pytest.raises(ValueError, match="3-D.*fits_slice"):
        read_plane("cube.fits")
    data, _, _, what, index = read_plane("cube.fits", fits_slice=[2])
    assert np.array_equal(data, cube[2]) and index == [2] and "fits_slice" in what
    with pytest.raises(ValueError, match="4 indices for the 3-D data"):
        read_plane("cube.fits", fits_slice=[0, 0, 0, 0])


def test_read_plane_reads_the_requested_hdu(in_tmp):
    """fits_hdu picks the extension; a missing HDU and an empty primary HDU are refused with a reason."""
    img = np.random.default_rng(3).random((6, 6))
    fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(img)]).writeto("two.fits")
    assert np.array_equal(read_plane("two.fits", fits_hdu=1)[0], img)
    with pytest.raises(ValueError, match="has 2 HDUs; fits_hdu = 2 does not exist"):
        read_plane("two.fits", fits_hdu=2)
    with pytest.raises(ValueError, match="0-D"):
        read_plane("two.fits", fits_hdu=0)


@pytest.mark.parametrize("kw, plane", [
    ({}, (5, 0, 0)),  # default: scattered starlight, first inclination and azimuth
    ({"mcfost_component": "total", "mcfost_inclination": 2, "mcfost_azimuth": 1}, (0, 1, 2)),
    ({"mcfost_component": "Q"}, (1, 0, 0)),
    ({"input_format": "mcfost", "mcfost_component": "direct_thermal"}, (6, 0, 0)),
    ({"fits_slice": [7, 1, 1]}, (7, 1, 1)),
])
def test_read_plane_mcfost_selects_component_azimuth_inclination(kw, plane, in_tmp):
    """MCFOST cubes (axes plane, azimuth, inclination, y, x): plane by name, the others by index."""
    cube = mcfost_cube("RT.fits.gz", 8)
    data, _, fmt, _, index = read_plane("RT.fits.gz", **kw)
    assert np.array_equal(data, cube[plane].astype(np.float64))
    assert fmt == "mcfost" and index == list(plane)


def test_read_plane_mcfost_disk_is_total_minus_direct_star(in_tmp):
    """'disk' is the total light minus the direct starlight, at the chosen inclination and azimuth."""
    cube = mcfost_cube("RT.fits.gz", 8)
    data, _, _, what, index = read_plane("RT.fits.gz", mcfost_component="disk", mcfost_inclination=2, mcfost_azimuth=1)
    assert np.allclose(data, cube[0, 1, 2].astype(np.float64)-cube[4, 1, 2])
    assert "MCFOST disk" in what and index == "total - direct_star at [1, 2]"


@pytest.mark.parametrize("kw, words", [
    ({"mcfost_inclination": 3}, "mcfost_inclination = 3 is out of range"),
    ({"mcfost_azimuth": 2}, "mcfost_azimuth = 2 is out of range"),
    ({"mcfost_inclination": -1}, "out of range"),
    ({"input_format": "image"}, "5-D"),
])
def test_read_plane_mcfost_refuses_bad_selections(kw, words, in_tmp):
    """Out-of-range MCFOST indices, or reading a cube as an image, raise ValueError with the reason."""
    mcfost_cube("RT.fits.gz", 8)
    with pytest.raises(ValueError, match=words):
        read_plane("RT.fits.gz", **kw)


def test_read_plane_mcfost_disk_needs_the_direct_star_plane(in_tmp):
    """'disk' on a 4-plane (Stokes only) MCFOST cube is refused: there is no direct_star plane to subtract."""
    write_fits("RT4.fits", np.ones((4, 1, 1, 6, 6)), FLUX_1="I", BUNIT="W.m-2.pixel-1")
    assert read_plane("RT4.fits", mcfost_component="Q")[4] == [1, 0, 0]
    with pytest.raises(ValueError, match="mcfost_component = 4 is out of range"):
        read_plane("RT4.fits", mcfost_component="disk")


def test_read_plane_mcfost_format_on_a_plain_image_is_refused(in_tmp):
    """input_format = 'mcfost' on 2-D data says it is not an MCFOST cube."""
    write_fits("img.fits", np.ones((6, 6)))
    with pytest.raises(ValueError, match="not an MCFOST image cube"):
        read_plane("img.fits", input_format="mcfost")


@sizes
@pytest.mark.parametrize("unit", ["Jy/arcsec^2", "Jy/pixel", "W/m^2/pixel"])
def test_read_target_converts_to_jy_per_arcsec2(n, unit, in_tmp):
    """read_target returns Jy/arcsec^2: Jy/pixel / psf_scale^2, nu F_nu per pixel * F_transfer."""
    img = np.random.default_rng(4).random((n, n))
    write_fits("disk.fits", img, BUNIT=unit, PIXSCALE=0.02, WAVELEN=1.6e-6)
    cfg = cfg_with(fits_filename="disk")
    out = read_target(cfg, verbose=False)
    factor = {"Jy/arcsec^2": 1.0, "Jy/pixel": 1/0.02**2, "W/m^2/pixel": cfg.F_transfer}[unit]
    assert out.shape == (n, n) and out.dtype == np.float64
    assert np.allclose(out, img*factor, rtol=1e-12)
    assert cfg.input_description == "disk.fits (image)" and cfg.used_slice == [] and cfg.target_pixel == n


def test_read_target_mcfost_cube_uses_its_header(in_tmp):
    """An MCFOST cube: plane 5 by default, nu F_nu * F_transfer, scale from CDELT2 (deg), wavelength from WAVE (um)."""
    cube = mcfost_cube("RT.fits.gz", 8)
    cfg = cfg_with(fits_filename="RT.fits.gz")
    out = read_target(cfg, verbose=False)
    assert np.isclose(cfg.psf_scale, 0.01) and cfg.wavelength == 0.8e-6 and cfg.input_unit == "W/m^2/pixel"
    assert np.isclose(cfg.F_transfer, 1e26*0.8e-6/299792458/0.01**2)
    assert np.allclose(out, cube[5, 0, 0].astype(np.float64)*cfg.F_transfer, rtol=1e-12)
    assert cfg.used_slice == [5, 0, 0]
    assert cfg.input_description == "RT.fits.gz (MCFOST scattered_star, inclination 0, azimuth 0)"


def test_read_target_reads_the_header_of_the_chosen_hdu(in_tmp):
    """With fits_hdu = 1 the unit and pixel scale come from that extension's header (not the primary's)."""
    img = np.random.default_rng(5).random((8, 8))
    ext = fits.ImageHDU(img, fits.Header([("BUNIT", "Jy/pixel"), ("PIXSCALE", 0.02)]))
    fits.HDUList([fits.PrimaryHDU(header=fits.Header([("PIXSCALE", 0.5)])), ext]).writeto("two.fits")
    cfg = cfg_with(fits_filename="two", fits_hdu=1)
    out = read_target(cfg, verbose=False)
    assert cfg.input_unit == "Jy/pixel" and np.isclose(cfg.psf_scale, 0.02)
    assert np.allclose(out, img/0.02**2)


def test_read_target_prints_the_input_summary(in_tmp, capsys):
    """verbose=True prints the input summary line first."""
    write_fits("disk.fits", np.ones((8, 8)), BUNIT="Jy/pixel", PIXSCALE=0.02)
    read_target(cfg_with(fits_filename="disk"))
    assert capsys.readouterr().out.startswith("Input: disk.fits (image), 8 x 8, Jy/pixel (header BUNIT)")


def test_read_target_without_a_unit_is_refused(in_tmp):
    """An image without BUNIT and without input_unit is refused (the message says how to fix it)."""
    write_fits("bare.fits", np.ones((8, 8)))
    with pytest.raises(ValueError, match="input_unit"):
        read_target(cfg_with(fits_filename="bare"), verbose=False)
    out = read_target(cfg_with(fits_filename="bare", input_unit="Jy/arcsec^2"), verbose=False)
    assert np.array_equal(out, np.ones((8, 8)))


def test_read_target_example_is_the_512_pixel_image():
    """The bundled example (fits_filename "") is a 512 x 512 image read at 32 lambda/D across."""
    cfg = from_dict({})
    img = read_target(cfg, verbose=False)
    assert img.shape == (512, 512) and img.dtype == np.float64 and np.all(np.isfinite(img))
    assert cfg.target_pixel == 512 and np.isclose(cfg.visual_range, 32)
    assert np.isclose(cfg.psf_scale, 1e-6/2.4*RAD*32/512)


def test_read_target_example_is_the_scattered_starlight_of_the_original_cube():
    """Regression: the example image is plane [5, 0, 0] of example_mcfost.fits times F_transfer (as before)."""
    cfg = from_dict({})
    img = read_target(cfg, verbose=False)
    with fits.open(os.path.join(DATADIR, "example_mcfost.fits")) as hdul:
        old = hdul[0].data[5, 0, 0].astype(np.float64)
    assert np.allclose(img, old*1e26*1e-6/299792458/cfg.psf_scale**2, rtol=1e-12)


# ---------------------------------------------------------------- selection: target symmetry

def test_example_target_is_measured_symmetric():
    """The bundled example passes the default symmetry test (asymmetry and corner light below 0.5%)."""
    with fits.open(os.path.join(DATADIR, "example.fits")) as hdul:
        pre = hdul[0].data.astype(float)
    asym, corners = image_asymmetry(pre)
    assert asym < 0.005 and corners < 0.005
    assert image_symmetry(pre)[0]


def test_image_asymmetry_of_an_empty_image_is_zero():
    """An image without light is reported as (0, 0)."""
    assert image_asymmetry(np.zeros((16, 16))) == (0.0, 0.0)


def test_image_asymmetry_with_light_only_in_the_corners():
    """Light only beyond the inscribed circle: no asymmetry inside it, corner fraction 1."""
    img = np.zeros((16, 16))
    img[0, 0] = img[15, 15] = 1.0
    asym, corners = image_asymmetry(img)
    assert asym == 0.0 and corners == 1.0


@sizes
def test_image_asymmetry_of_symmetric_targets_is_round_off(n):
    """Exactly symmetric targets about the star at (N-1)/2 have asymmetry at round-off level."""
    symmetric, _ = targets(n)
    for name, img in symmetric.items():
        asym, corners = image_asymmetry(img)
        assert asym < 1e-12 and corners < 0.01, name


@sizes
def test_image_asymmetry_uses_the_absolute_value(n):
    """The sign of the image does not matter (|pixel| is measured)."""
    _, asymmetric = targets(n)
    img = asymmetric["half ring"]
    assert image_asymmetry(-img) == pytest.approx(image_asymmetry(img))


@sizes
def test_image_symmetry_auto_detects_symmetric_and_asymmetric_targets(n):
    """'auto' (and None) measure the target: symmetric ones pass, inclined/off-centre/arcs/corners do not."""
    symmetric, asymmetric = targets(n)
    for name, img in symmetric.items():
        sym, why = image_symmetry(img, "auto", 0.01)
        assert sym and why.startswith("target symmetric: asymmetry"), name
        assert image_symmetry(img, None)[0], name
    for name, img in asymmetric.items():
        sym, why = image_symmetry(img, "auto", 0.01)
        assert not sym and why.startswith("target not symmetric"), name


@sizes
def test_image_symmetry_tolerance_is_inclusive_for_asymmetry_and_corner_light(n):
    """Symmetric means asymmetry <= tol and corner fraction <= tol (each one alone can fail the test)."""
    _, asymmetric = targets(n)
    inclined = asymmetric["inclined"]
    asym, corners = image_asymmetry(inclined)
    assert corners < asym
    assert image_symmetry(inclined, "auto", asym)[0] and not image_symmetry(inclined, "auto", 0.999*asym)[0]
    lit = asymmetric["light in the corners"]
    asym, corners = image_asymmetry(lit)
    assert asym < 0.01*corners  # ~symmetric inside the inscribed circle; the extra light is in a corner
    assert image_symmetry(lit, "auto", corners)[0] and not image_symmetry(lit, "auto", 0.999*corners)[0]


@sizes
def test_image_symmetry_reports_round_off_asymmetry_as_zero_percent(n):
    """The reason line gives the asymmetry and corner light in percent, round-off shown as 0%."""
    symmetric, _ = targets(n)
    sym, why = image_symmetry(symmetric["smooth ring"], "auto", 0.01)
    assert sym and why.startswith("target symmetric: asymmetry 0%, light in the corners") and why.endswith("(tol 1%)")


@sizes
def test_image_symmetry_user_setting_overrides_the_measurement(n):
    """True/False from the parameters are used as given, with the reason 'set in the parameters'."""
    symmetric, asymmetric = targets(n)
    assert image_symmetry(asymmetric["half ring"], True) == (True, "target symmetric: set in the parameters")
    assert image_symmetry(symmetric["smooth ring"], False) == (False, "target not symmetric: set in the parameters")


# ---------------------------------------------------------------- selection: PSF symmetry

def test_psf_symmetry_of_the_vortex_comes_from_the_model():
    """'auto' and None use the vortex model's declaration (symmetric) without measuring."""
    o = Optics.vortex(2, 32, 4.0)
    assert psf_symmetry(o) == (True, "PSFs symmetric: vortex model")
    assert psf_symmetry(o, setting=None)[0] is True


def test_psf_symmetry_user_setting_notes_a_contradicting_model():
    """A user setting wins; when the model says otherwise the reason says so."""
    o = Optics.vortex(2, 32, 4.0)
    sym, why = psf_symmetry(o, setting=False)
    assert sym is False and "set in the parameters" in why and "vortex model says symmetric" in why
    assert psf_symmetry(o, setting=True) == (True, "PSFs symmetric: set in the parameters")


@pytest.mark.parametrize("kind, expected", [("vortex_undeclared", True), ("vortex_offset_lyot", False)])
def test_psf_symmetry_is_measured_when_the_model_does_not_declare_it(kind, expected):
    """Undeclared models are measured: the plain vortex is symmetric, a shifted Lyot stop is not."""
    o = Optics(kind, 32, 4.0, Optics.vortex(2, 32, 4.0).params, module='extra_backends')
    sym, why = psf_symmetry(o, num_cores=1)
    assert sym is expected and "measured" in why


# ---------------------------------------------------------------- selection: method choice

@pytest.mark.parametrize("image_sym, psf_sym, allowed, method, mode", [
    (True, True, True, "circular", None),
    (False, True, True, "spatial", "rotate"),
    (True, False, False, "spatial", "direct"),
    (False, False, False, "spatial", "direct"),
    (False, False, True, "spatial", "direct"),  # rotate-capable optics, but the PSFs were said not symmetric
    (False, True, False, "spatial", "direct"),
])
def test_choose_method_auto(image_sym, psf_sym, allowed, method, mode):
    """'auto': circular only if both are symmetric; rotate only for symmetric PSFs on optics allowing it."""
    m, sm, msg = choose_method("auto", "auto", image_sym, psf_sym, allowed, ("target ...", "PSFs ..."))
    assert (m, sm) == (method, mode)
    assert msg[0] == "Method: %s%s (auto; target ...; PSFs ...)" % (method, "" if mode is None else ", %s mode" % mode)
    assert not any("WARNING" in x for x in msg)


def test_choose_method_explicit_settings_are_kept_with_warnings():
    """Explicit method / spatial_mode are kept; contradicting the symmetries gives a WARNING line."""
    m, sm, msg = choose_method("circular", "auto", False, True, True)
    assert (m, sm) == ("circular", None) and msg[1].startswith("WARNING: the circular method is only exact")
    assert "set in the parameters" in msg[0]
    m, sm, msg = choose_method("spatial", "auto", True, True, True)
    assert (m, sm) == ("spatial", "rotate") and len(msg) == 1
    m, sm, msg = choose_method("spatial", "rotate", False, False, True)
    assert sm == "rotate" and msg[1] == "WARNING: spatial_mode 'rotate' assumes symmetric PSFs"
    assert choose_method("spatial", "direct", True, True, True)[:2] == ("spatial", "direct")
    assert choose_method("circular", "rotate", True, True, False)[:2] == ("circular", None)


@pytest.mark.parametrize("kwargs, words", [
    (dict(method="fast"), "method must be 'auto', 'circular' or 'spatial'"),
    (dict(method="spatial", spatial_mode="fast"), "spatial_mode must be 'auto', 'rotate' or 'direct'"),
    (dict(method="spatial", spatial_mode="rotate", rotate_allowed=False), "rotationally symmetric"),
    (dict(method="auto", spatial_mode="rotate", image_sym=False, rotate_allowed=False), "rotationally symmetric"),
])
def test_choose_method_rejects_unknown_or_impossible_settings(kwargs, words):
    """Unknown method / spatial_mode, or 'rotate' on optics that cannot use it, raise ValueError."""
    with pytest.raises(ValueError, match=words):
        choose_method(**kwargs)
