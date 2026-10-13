"""
Unit tests for stage F0 (tests/scene_stage0_params.py): the parameter system.

Modules covered:
    toycoronagraph/params.py       Param / convert (units, "auto", choices, lambda/D), reading and layout
                                   checks, Params.check / resolve (conflict rules, grid resolution,
                                   provenance), automatic names, used_toml, include, the old layout
    toycoronagraph/params_docs.py  new_config, json_schema, describe, migrate, reference_markdown
    toycoronagraph/__main__.py     the command line (describe, new, schema)

Every test builds small parameter dicts (or tiny files in its own temporary folder) and runs in
milliseconds; the stage script keeps the generated example files and the conflict table.
"""
import copy
import os

import numpy as np
import pytest
import tomllib
from astropy import units as u

import toycoronagraph.params as P
import toycoronagraph.params_docs as D
from toycoronagraph.__main__ import main as cli

# The stage script's example: a wind around eps Eridani (every value set by the user)
EPS_ERI = {
    "observation": {"wavelength": 5e-7, "diameter": 6.0},
    "grid": {"sampling": 4, "fov": 2.0},
    "star": {"distance": 3.22, "radius": 0.73},
    "source": [{"kind": "wind", "mdot": 6e-13}],
}
LAMBDA_D = 5e-7/6.0*P.RAD_TO_ARCSEC  # arcsec, of EPS_ERI
RING = {"kind": "ring", "module": "extra_sources"}          # a physical test module (tests/extra_sources.py)
BLOB = {"kind": "cme_blob", "module": "extra_sources"}      # a relative (contrast) test module


def _data(source=None, **sections):
    """
    A copy of EPS_ERI with some keys changed.

    Args:
        source (list of dict): replaces the [[source]] tables; None = keep the wind.
        **sections: {section: {key: value}} merged into EPS_ERI; a value None removes the key.

    Returns:
        dict: the parameter data, in the file layout.
    """
    data = copy.deepcopy(EPS_ERI)
    for name, keys in sections.items():
        table = data.setdefault(name, {})
        for key, value in keys.items():
            if value is None:
                table.pop(key, None)
            else:
                table[key] = value
    if source is not None:
        data["source"] = copy.deepcopy(source)
    return data


def _resolved(data, provided=()):
    """
    Check and resolve parameter data.

    Args:
        data (dict): the parameters in the file layout.
        provided (iterable of tuple): (key, value, origin) supplied as a source would.

    Returns:
        Params: after resolve().
    """
    p = P.from_data(data)
    for key, value, origin in provided:
        p.provide(key, value, origin)
    p.resolve()
    return p


def _error(func, *args, **kwargs):
    """
    The message of the ParamError raised by func(*args, **kwargs) (the test fails if none is raised).

    Returns:
        str: the error text.
    """
    with pytest.raises(P.ParamError) as err:
        func(*args, **kwargs)
    return str(err.value)


def _write(path, text):
    """Write text to a file and return its path (str)."""
    with open(path, "w") as f:
        f.write(text)
    return str(path)


# ====================================================================== Param / convert

def test_convert_number_is_float_in_canonical_unit():
    """A plain number (int or float) is returned as a float in the declared unit."""
    value = P.convert(P._param("star.distance"), 3, "[star]")
    assert value == 3.0 and type(value) is float


@pytest.mark.parametrize("key, text, expected", [
    ("star.distance", "10.5 lyr", (10.5*u.lyr).to_value(u.pc)),
    ("star.radius", "1 Rsun", 1.0),                      # "Rsun" is accepted for R_sun
    ("grid.pixel_scale", "10 mas", 0.01),
    ("observation.wavelength", "500 nm", 5e-7),
    ("observation.diameter", "650 cm", 6.5),
])
def test_convert_unit_string_to_canonical_unit(key, text, expected):
    """A string with a unit is converted to the parameter's canonical unit with astropy."""
    assert np.isclose(P.convert(P._param(key), text, "[x]"), expected, rtol=1e-12)


def test_convert_msun_alias_in_module_parameter():
    """'Msun' is read as M_sun (a wind's mass-loss rate in M_sun/yr)."""
    from toycoronagraph.sources import get_source
    mdot = {p.name: p for p in get_source("wind").params}["mdot"]
    assert np.isclose(P.convert(mdot, "6e-13 Msun/yr", "wind"), 6e-13, rtol=1e-12)


def test_convert_wrong_unit_names_where_key_and_canonical_unit():
    """A unit of the wrong dimension: '<where> <key> = <value>: cannot convert to <unit> ... (expected ...)'."""
    with pytest.raises(ValueError) as err:
        P.convert(P._param("star.distance"), "3 kg", "[star]")
    text = str(err.value)
    assert text.startswith('[star] distance = "3 kg": cannot convert to pc') and text.endswith("(expected number [pc])")


@pytest.mark.parametrize("text", ["32 lambda/D", "32 lam/d", "  32 λ/D ", "32 LAMBDA / D"])
def test_convert_lambda_d_spellings(text):
    """'N lambda/D' (lam/D, the Greek letter, any case and spacing) gives a LambdaD of N."""
    value = P.convert(P._param("grid.fov"), text, "[grid]")
    assert isinstance(value, P.LambdaD) and float(value) == 32.0 and repr(value) == '"32 lambda/D"'


def test_convert_lambda_d_refused_where_not_declared():
    """'N lambda/D' is only allowed for parameters declared with lambda_d (not for a distance)."""
    with pytest.raises(ValueError):
        P.convert(P._param("star.distance"), "3 lambda/D", "[star]")


def test_convert_accepts_its_own_lambda_d_output():
    """Converting a converted value again keeps the lambda/D size (make() re-checks checked sections)."""
    fov = P._param("grid.fov")
    once = P.convert(fov, "8 lambda/D", "[grid]")
    assert isinstance(P.convert(fov, once, "[grid]"), P.LambdaD)


@pytest.mark.parametrize("key, value, words", [
    ("star.distance", 0, "must be > 0"),              # exclusive minimum
    ("star.distance", -1.0, "must be > 0"),
    ("grid.size", 2, "must be >= 3"),                 # inclusive minimum
    ("grid.fov", "-4 lambda/D", "must be > 0"),       # bounds apply to lambda/D sizes too
])
def test_convert_bounds(key, value, words):
    """Numbers outside min / max are refused with the bound in the message."""
    with pytest.raises(ValueError, match=words):
        P.convert(P._param(key), value, "[x]")


def test_convert_maximum():
    """A number above max is refused ('must be <= max')."""
    lyot = {p.name: p for p in P.CORONAGRAPHS["vortex"][1]}["lyot_mask_size"]
    with pytest.raises(ValueError, match="must be <= 1"):
        P.convert(lyot, 1.5, "[[coronagraph]]")


@pytest.mark.parametrize("key, value, words", [
    ("star.distance", True, "not a number"),          # a bool is not a number
    ("grid.size", True, "not an integer"),
    ("grid.size", 3.0, "not an integer"),             # a float is not an integer
    ("grid.resample", 1, "not true/false"),           # 1 is not a bool
    ("run.name", 3, "not text"),
    ("star.distance", "abc", "not a number"),
])
def test_convert_type_errors(key, value, words):
    """Values of the wrong type are refused, bools never pass as numbers."""
    with pytest.raises(ValueError, match=words):
        P.convert(P._param(key), value, "[x]")


def test_convert_choices_compare_types():
    """Choices compare with their type: True is a choice, 1 is not (1 == True in Python)."""
    symmetric = {p.name: p for p in P.CORONAGRAPHS["vortex"][1]}["symmetric"]
    assert P.convert(symmetric, True, "c") is True
    with pytest.raises(ValueError, match="not allowed"):
        P.convert(symmetric, 1, "c")


def test_convert_choice_misspelled_did_you_mean():
    """A misspelled choice gets the closest allowed one."""
    with pytest.raises(ValueError, match='did you mean "separate"'):
        P.convert(P._param("method.combine"), "seperate", "[method]")


def test_convert_auto_only_where_allowed():
    """'auto' is returned as is where allowed, refused elsewhere."""
    assert P.convert(P._param("grid.size"), "auto", "[grid]") == P.AUTO
    with pytest.raises(ValueError):
        P.convert(P._param("star.distance"), "auto", "[star]")


def test_convert_lists():
    """float_list items become floats; int_list refuses bools; auto_word is accepted as is."""
    assert P.convert(P.Param("x", type="float_list"), [1, 2], "w") == [1.0, 2.0]
    with pytest.raises(ValueError, match="not a list of the right type"):
        P.convert(P.Param("x", type="int_list"), [1, True], "w")
    assert P.convert(P.Param("x", type="str_list", auto_word="all"), "all", "w") == "all"


@pytest.mark.parametrize("key, words", [
    ("grid.fov", 'number [arcsec] or "N lambda/D"'),
    ("observation.wavelength", '"auto" or number [m]'),
    ("method.combine", 'one of "separate", "together"'),
    ("grid.resample", "true/false"),
])
def test_describe_type(key, words):
    """describe_type() words the allowed values as the messages and docs show them."""
    assert P._param(key).describe_type() == words


def test_toml_values():
    """_toml writes TOML: floats keep a decimal point, bools in lower case, LambdaD as text."""
    assert P._toml(1.0) == "1.0" and P._toml(1e20) == "1e+20" and P._toml(3) == "3"
    assert P._toml(True) == "true" and P._toml(P.LambdaD(32)) == '"32 lambda/D"'
    assert tomllib.loads("a = %s" % P._toml({"b": [1.0, "x"]})) == {"a": {"b": [1.0, "x"]}}


# ====================================================================== layout and module tables

def test_unknown_section_suggests_closest():
    """An unknown section is refused with the closest known one."""
    assert "unknown section [grdi] (did you mean [grid]?)" in _error(P.from_data, {"grdi": {}})


def test_source_must_be_a_list_of_tables():
    """[source] (a single table) is refused: use [[source]]."""
    assert "use [[source]] (a list of tables)" in _error(P.from_data, {"source": {"kind": "wind"}})


def test_section_must_be_a_table():
    """A core section given as a value is refused."""
    assert "[grid] must be a table" in _error(P.from_data, {"grid": 3})


def test_source_without_kind_lists_known_kinds():
    """A [[source]] without kind is refused with the list of kinds."""
    text = _error(P.from_data, {"source": [{"mdot": 1e-13}]})
    assert "kind is missing" in text and "wind" in text and "fits_image" in text


def test_shared_key_in_source_points_to_its_home():
    """A shared key inside a [[source]] table is refused with a pointer to its section."""
    data = _data(source=[{"kind": "wind", "mdot": 6e-13, "distance": 3.22}])
    assert "distance is a shared setting: set it in [star]" in _error(P.from_data, data)


def test_misspelled_key_did_you_mean():
    """A misspelled module key gets the closest known key."""
    data = _data(source=[{"kind": "wind", "mdto": 6e-13}])
    assert "unknown key mdto (did you mean mdot?)" in _error(P.from_data, data)


def test_misspelled_kind_did_you_mean():
    """A misspelled kind gets the closest known kind."""
    data = _data(source=[{"kind": "wnd", "mdot": 6e-13}])
    assert "unknown kind 'wnd' (did you mean wind?)" in _error(P.from_data, data)


def test_unknown_kind_without_match_asks_for_module():
    """An unknown kind with no close match points to `module` for your own modules."""
    data = _data(source=[{"kind": "zzqqxx"}])
    assert 'set module = "package.module"' in _error(P.from_data, data)


def test_module_that_cannot_be_imported():
    """module = an unimportable name raises ValueError naming the module."""
    with pytest.raises(ValueError, match="could not be imported"):
        P.from_data({"source": [{"kind": "nothing_here", "module": "no_such_module_xyz"}]})


def test_all_problems_listed_at_once():
    """Several mistakes give one ParamError listing each one on its own line."""
    data = _data(source=[{"kind": "wind", "mdto": 6e-13, "distance": 3.2}], grdi={})
    text = _error(P.from_data, data)
    assert "[grdi]" in text and "mdto" in text and "shared setting" in text
    assert text.startswith("Problems in the parameters:\n")


# ====================================================================== check(): missing and unused

def test_required_module_parameter_missing():
    """A module parameter without a default must be set ('needs mdot')."""
    assert "[[source]] wind needs mdot (M_sun/yr)" in _error(P.from_data, _data(source=[{"kind": "wind"}]))


def test_required_shared_key_missing():
    """A wind without [star] distance: 'wind needs [star] distance (pc)'."""
    assert "wind needs [star] distance (pc)" in _error(P.from_data, _data(star={"distance": None}))


def test_mixed_families_need_the_star_flux():
    """Contrast and Jy components together without temperature / flux_jy: the star's flux is needed."""
    data = _data(source=[EPS_ERI["source"][0], {"kind": "fits_image", "file": "disk.fits", "input_unit": "Jy/arcsec^2"}])
    text = _error(P.from_data, data)
    assert "the star's flux is needed" in text and "wind: contrast" in text and "fits_image: physical" in text


def test_blackbody_needs_radius_and_distance():
    """With the flux from temperature, the blackbody needs [star] radius too."""
    data = {"grid": {"size": 9, "pixel_scale": 0.1}, "star": {"distance": 3.0, "temperature": 5000},
            "source": [RING, BLOB]}
    assert "the blackbody flux (temperature) needs [star] radius" in _error(P.from_data, data)


def test_flux_jy_alone_combines_families():
    """flux_jy is enough to combine contrast and Jy components (no radius needed)."""
    data = {"grid": {"size": 9, "pixel_scale": 0.1}, "star": {"distance": 3.0, "flux_jy": 2.0}, "source": [RING, BLOB]}
    p = P.from_data(data)
    assert p.star_flux_needed() and "unit conversion" in p.used_by("star.flux_jy")
    assert "unit conversion" not in p.used_by("star.distance")


def test_unused_star_key_is_a_note():
    """A [star] key no enabled component uses gives a note (not an error)."""
    p = P.from_data(_data(star={"temperature": 5080}))
    assert any(n.startswith("Note: [star] temperature is not used") for n in p.notes)


def test_display_unit_of_other_family_needs_star_flux():
    """[display] unit in contrast for Jy components needs the star's flux."""
    data = {"grid": {"size": 9, "pixel_scale": 0.1}, "display": {"unit": "contrast_pixel"}, "source": [RING]}
    assert "the star's flux is needed" in _error(P.from_data, data)


def test_disabled_component_needs_nothing():
    """A disabled wind does not need [star] distance."""
    data = _data(star={"distance": None}, source=[{"kind": "wind", "mdot": 6e-13, "enabled": False}])
    p = P.from_data(data)
    assert p.enabled() == [] and len(p.components) == 1


def test_check_twice_gives_the_same_notes():
    """Regression: check() again does not repeat its notes."""
    p = P.from_data(_data(star={"temperature": 5080}))
    first = list(p.notes)
    p.check()
    p.check()
    assert p.notes == first and p.load_notes == first and len(first) == 1


def test_resolve_and_check_again_give_the_same_notes():
    """Regression: resolve() and check() again give the same notes (none added, none lost)."""
    p = _resolved(_data(grid={"pixel_scale": LAMBDA_D/4}, star={"temperature": 5080}))
    notes = list(p.notes)
    assert any("agree" in n for n in notes) and any("temperature is not used" in n for n in notes)
    p.check()
    p.resolve()
    p.resolve()
    assert p.notes == notes


def test_resolve_again_uses_the_new_provided_values():
    """resolve() starts again from the user's settings: new provided values replace the old result."""
    p = P.from_data({"source": [{"kind": "fits_image", "file": "disk.fits"}]})
    p.provided = [("grid.size", 64, "image"), ("grid.pixel_scale", 0.01, "header A")]
    p.resolve()
    p.provided = [("grid.size", 32, "image"), ("grid.pixel_scale", 0.02, "header B")]
    p.resolve()
    assert (p.get("grid.size"), p.get("grid.pixel_scale"), p.origin("grid.pixel_scale")) == (32, 0.02, "header B")


# ====================================================================== resolve(): conflicts

def test_distance_conflict_is_an_error():
    """Provided distance off by more than 0.5%: CONFLICT error naming both values and the fix."""
    text = _error(_resolved, _data(), [("star.distance", 3.2, "cme_imaging Star")])
    assert "CONFLICT star.distance: 3.22 (set in parameters) vs 3.2 (cme_imaging Star)" in text
    assert 'on_conflict = "warn"' in text


def test_on_conflict_warn_accepts_with_a_warning():
    """[checks] on_conflict = "warn" turns the conflict into an accepted warning; the user's value wins."""
    p = _resolved(_data(checks={"on_conflict": "warn"}), [("star.distance", 3.2, "cme_imaging Star")])
    assert any(n.startswith("WARNING (accepted") for n in p.notes) and p.get("star.distance") == 3.22


def test_distance_within_tolerance_is_silent():
    """A provided distance within 0.5% gives no message."""
    p = _resolved(_data(), [("star.distance", 3.225, "cme_imaging Star")])
    assert p.notes == [] and p.get("star.distance") == 3.22


def test_wavelength_mismatch_is_a_warning():
    """A header wavelength differing from the user's is a warning (the run goes on)."""
    p = _resolved(_data(), [("observation.wavelength", 6e-7, "header WAVELEN of disk.fits")])
    assert any(n.startswith("WARNING: observation.wavelength") and "6e-07" in n for n in p.notes)


def test_provided_pixel_scale_is_not_a_conflict():
    """An image's own pixel scale differing from [grid] is no conflict (the scene resamples)."""
    p = _resolved(_data(), [("grid.pixel_scale", 0.01, "header PIXSCALE of disk.fits")])
    assert p.notes == [] and np.isclose(p.get("grid.pixel_scale"), LAMBDA_D/4)


def test_pixel_scale_and_sampling_disagree():
    """[grid] pixel_scale and sampling that disagree: 'set only one'."""
    assert "set only one" in _error(_resolved, _data(grid={"pixel_scale": 0.01}))


def test_pixel_scale_and_sampling_agree_is_a_note():
    """[grid] pixel_scale and sampling that agree: a note that one is enough."""
    p = _resolved(_data(grid={"pixel_scale": LAMBDA_D/4}))
    assert "Note: [grid] pixel_scale and sampling agree; one of them is enough" in p.notes


def test_size_pixel_scale_and_fov_disagree():
    """[grid] size x pixel_scale far from fov: 'set two of the three'."""
    data = _data(grid={"sampling": None, "pixel_scale": 0.01, "size": 101})
    assert "set two of the three" in _error(_resolved, data)


def test_size_derived_from_pixel_scale_and_fov_is_accepted_when_set():
    """Writing out the size that pixel_scale = 0.01 and fov = 2.0 give (201, made odd) is accepted."""
    derived = _resolved(_data(grid={"sampling": None, "pixel_scale": 0.01})).get("grid.size")
    assert derived == 201
    p = _resolved(_data(grid={"sampling": None, "pixel_scale": 0.01, "size": derived}))
    assert p.get("grid.size") == 201


def test_temperature_and_flux_jy_disagree():
    """[star] temperature and flux_jy whose blackbody disagrees: 'set only one'."""
    text = _error(_resolved, _data(star={"temperature": 5080, "flux_jy": 1.0}))
    assert "[star] flux_jy = 1 Jy but the blackbody of temperature = 5080 K" in text and "set only one" in text


def test_temperature_and_flux_jy_without_radius():
    """temperature and flux_jy together without a radius to compare them: 'not both'."""
    data = {"star": {"temperature": 5080, "flux_jy": 1.0, "distance": 3.0}, "grid": {"size": 9, "pixel_scale": 0.1}}
    assert "set temperature OR flux_jy, not both" in _error(_resolved, data)


def test_temperature_and_matching_flux_jy_accepted():
    """temperature and the flux_jy of its blackbody (within 0.1%) are accepted."""
    bb = P.blackbody_flux_jy(5080, 0.73, 3.22, 5e-7)
    p = _resolved(_data(star={"temperature": 5080, "flux_jy": bb*1.0005}))
    assert np.isclose(p.get("star.flux_jy"), bb*1.0005)


# ====================================================================== resolve(): values, grid, provenance

def test_resolution_order_user_then_provided_then_default():
    """The user's value wins over a provided one, which wins over the default."""
    data = {"observation": {"wavelength": 1.6e-6}, "source": [{"kind": "fits_image", "file": "disk.fits"}]}
    p = _resolved(data, [("grid.size", 33, "image size"), ("grid.pixel_scale", 0.01, "header PIXSCALE"),
                         ("observation.wavelength", 1.6e-6, "header WAVELEN")])
    assert p.origin("observation.wavelength") == "set in parameters"
    assert (p.get("grid.size"), p.origin("grid.size")) == (33, "image size")
    assert (p.get("observation.diameter"), p.origin("observation.diameter")) == (2.4, "default")


@pytest.mark.parametrize("fov", [2.0, 2.003, 1.0, 0.55])
def test_grid_from_sampling_and_fov_is_odd(fov):
    """sampling + fov: pixel scale lambda/D / sampling, size made odd (within 1.5 pixels of fov)."""
    p = _resolved(_data(grid={"fov": fov}))
    scale, size = p.get("grid.pixel_scale"), p.get("grid.size")
    assert np.isclose(scale, LAMBDA_D/4) and size % 2 == 1 and isinstance(size, int)
    assert abs(size*scale - fov) <= 1.5*scale and p.origin("grid.size") == "from fov / pixel_scale (odd)"
    assert np.isclose(p.get("grid.fov"), size*scale)


def test_grid_fov_in_lambda_d():
    """fov = "N lambda/D" is converted with the wavelength and diameter; its origin says so."""
    p = _resolved({"observation": {"wavelength": 1e-6, "diameter": 2.0}, "grid": {"sampling": 2, "fov": "8 lambda/D"}})
    lam_d = 1e-6/2.0*P.RAD_TO_ARCSEC
    assert p.get("grid.size") == 17 and np.isclose(p.get("grid.fov"), 17*lam_d/2)
    assert p.origin("grid.fov") == "set in parameters (8 lambda/D)"


def test_image_supplies_the_auto_grid():
    """With [grid] "auto", the image's size and header pixel scale are used, with their origins."""
    p = _resolved({"source": [{"kind": "fits_image", "file": "disk.fits"}]},
                  [("grid.size", 255, "image size of disk.fits"), ("grid.pixel_scale", 0.01, "header PIXSCALE of disk.fits"),
                   ("observation.wavelength", 1.6e-6, "header WAVELEN of disk.fits")])
    assert (p.get("grid.size"), p.get("grid.pixel_scale"), p.get("observation.wavelength")) == (255, 0.01, 1.6e-6)
    assert p.origin("grid.pixel_scale") == "header PIXSCALE of disk.fits"


def test_size_only_keeps_the_image_field():
    """[grid] size alone (resolution study): the image's field is kept, the pixel scale follows."""
    p = _resolved({"grid": {"size": 64}, "source": [{"kind": "fits_image", "file": "disk.fits"}]},
                  [("grid.size", 128, "image"), ("grid.pixel_scale", 0.01, "header")])
    assert np.isclose(p.get("grid.pixel_scale"), 0.02) and np.isclose(p.get("grid.fov"), 1.28)
    assert p.origin("grid.pixel_scale") == "the image's field / size"


@pytest.mark.parametrize("n_image, expected", [(128, 98), (127, 99)])
def test_pixel_scale_only_keeps_the_image_field_and_parity(n_image, expected):
    """[grid] pixel_scale alone: the size covers the image's field, with the image's parity."""
    p = _resolved({"grid": {"pixel_scale": 0.013}, "source": [{"kind": "fits_image", "file": "disk.fits"}]},
                  [("grid.size", n_image, "image"), ("grid.pixel_scale", 0.01, "header")])
    assert p.get("grid.size") == expected and p.origin("grid.size") == "the field of view of the image"


def test_nothing_gives_wavelength_or_pixel_scale():
    """No wavelength and no pixel scale anywhere: 1e-6 m and 32 lambda/D assumed, each with a note."""
    p = _resolved({"source": [{"kind": "fits_image", "file": "disk.fits"}]}, [("grid.size", 64, "image")])
    assert p.get("observation.wavelength") == 1e-6
    assert np.isclose(p.get("grid.fov"), 32*1e-6/2.4*P.RAD_TO_ARCSEC)
    assert any("1e-06 m assumed" in n for n in p.notes) and any("32 lambda/D" in n for n in p.notes)


def test_not_enough_to_define_the_grid():
    """No grid settings and no image: an error saying what to set."""
    assert "not enough to define the image" in _error(_resolved, {})


def test_origin_texts():
    """origin(): 'set in parameters' (Python), 'default', or the provided origin."""
    p = _resolved(_data())
    assert p.origin("star.distance") == "set in parameters"
    assert p.origin("method.combine") == "default" and p.origin("grid.pixel_scale") == "from sampling = 4"


def test_table_shows_values_origins_and_users():
    """table(): one line per resolved key with its unit, origin and users."""
    table = _resolved(_data()).table()
    lines = table.splitlines()
    assert lines[0].split() == ["quantity", "value", "from", "used", "by"]
    distance = next(line for line in lines if line.startswith("star.distance"))
    assert "3.22 pc" in distance and "set in parameters" in distance and distance.rstrip().endswith("wind")
    assert "from sampling = 4" in table and "the PSFs" in table and "the image grid" in table


# ====================================================================== automatic names, coronagraphs

def test_automatic_names_numbered_per_kind():
    """Regression: automatic names count per kind (ring, cme_blob, ring_2), not per position overall."""
    data = {"grid": {"size": 9, "pixel_scale": 0.1}, "star": {"distance": 3.0, "flux_jy": 2.0},
            "source": [RING, BLOB, RING, dict(RING, name="outer"), RING]}
    assert [c.name for c in P.from_data(data).components] == ["ring", "cme_blob", "ring_2", "outer", "ring_4"]


def test_duplicate_component_names_refused():
    """Two components with the same name are refused."""
    data = {"grid": {"size": 9, "pixel_scale": 0.1}, "source": [dict(RING, name="a"), dict(RING, name="a")]}
    assert "two components are named 'a'" in _error(P.from_data, data)


@pytest.mark.parametrize("name", ["total", "input"])
def test_reserved_component_names_refused(name):
    """"total" and "input" name the sum and the input row in plots and saved files: refused as component names."""
    data = {"grid": {"size": 9, "pixel_scale": 0.1}, "source": [dict(RING, name=name)]}
    assert "cannot be named %r" % name in _error(P.from_data, data)


def test_default_coronagraph_and_auto_names():
    """No [[coronagraph]]: one vortex of charge 2; automatic names are vortex_c<charge>."""
    assert [c.name for c in P.from_data(_data()).coronagraphs] == ["vortex_c2"]
    p = P.from_data(dict(_data(), coronagraph=[{"charge": 6}, {"type": "vortex", "charge": 2, "name": "mine"}]))
    assert [c.name for c in p.coronagraphs] == ["vortex_c6", "mine"] and p.coronagraphs[0].values == {"charge": 6}


@pytest.mark.parametrize("tables, words", [
    ([{"charge": 3}], "a vortex charge must be even"),
    ([{"type": "lyot"}], "(lyot): unknown type; known: vortex"),
    ([{"charge": 2}, {"charge": 2}], "two [[coronagraph]] tables have the same name"),
    ([{"charge": 2, "lyot": 0.9}], "unknown key lyot"),
])
def test_coronagraph_table_errors(tables, words):
    """Odd charge, unknown type, duplicate names and unknown keys in [[coronagraph]] are refused."""
    assert words in _error(P.from_data, dict(_data(), coronagraph=tables))


# ====================================================================== used_toml

def test_used_toml_reproduces_the_resolved_values():
    """used_toml() is valid TOML and, read back, resolves to the same values."""
    p = _resolved(_data())
    back = _resolved(tomllib.loads(p.used_toml()))
    assert {k: v for k, (v, _) in back.resolved.items()} == {k: v for k, (v, _) in p.resolved.items()}
    assert [c.name for c in back.components] == ["wind"] and back.coronagraphs[0].name == "vortex_c2"


def test_used_toml_keeps_grid_resample():
    """Regression: settings of shared sections that are not resolved values ([grid] resample) are kept."""
    used = _resolved(_data(grid={"resample": False})).used_toml()
    assert tomllib.loads(used)["grid"]["resample"] is False


def test_used_toml_writes_the_module_of_own_kinds():
    """A component of your own module gets its module = line (built-in kinds do not)."""
    data = {"grid": {"size": 9, "pixel_scale": 0.1}, "source": [RING, {"kind": "planet", "pos": [0.1, 0.0], "brightness": 1e-5}]}
    sources = tomllib.loads(_resolved(data).used_toml())["source"]
    assert sources[0]["module"] == "extra_sources" and "module" not in sources[1]


def test_used_toml_records_provided_origins_and_units():
    """Each resolved value is written with its unit and origin as a comment."""
    used = _resolved({"source": [{"kind": "fits_image", "file": "disk.fits"}]},
                     [("grid.size", 33, "image size of disk.fits"), ("grid.pixel_scale", 0.01, "header PIXSCALE of disk.fits")]
                     ).used_toml()
    line = next(line for line in used.splitlines() if line.startswith("pixel_scale"))
    assert "# arcsec; header PIXSCALE of disk.fits" in line
    assert next(line for line in used.splitlines() if line.startswith("size")).endswith("# image size of disk.fits")


def test_used_toml_lists_every_setting_of_components_and_coronagraphs():
    """Components carry all their settings (defaults too), coronagraphs all their parameters."""
    data = tomllib.loads(_resolved(_data()).used_toml())
    wind = data["source"][0]
    assert wind["v"] == 400.0 and wind["enabled"] is True and wind["scale"] == 1.0 and wind["method"] == "auto"
    assert data["coronagraph"][0] == {"type": "vortex", "name": "vortex_c2", "charge": 2, "lyot_mask_size": 0.95,
                                      "symmetric": "auto"}


# ====================================================================== files: load, include

def test_load_records_the_file_and_leaves_it_unchanged(in_tmp):
    """load(): values set in the file say so; the file is byte-identical after resolve and write_used."""
    path = _write(in_tmp/"run.toml", '[star]\ndistance = "10.5 lyr"\nradius = 0.73\n[observation]\nwavelength = 5e-7\n'
                                    'diameter = 6.0\n[grid]\nsampling = 4\nfov = 2.0\n[[source]]\nkind = "wind"\nmdot = 6e-13\n')
    before = open(path, "rb").read()
    p = P.load(path)
    p.resolve()
    p.write_used(str(in_tmp/"run_used.toml"))
    assert open(path, "rb").read() == before and p.origin("star.distance") == "set in run.toml"
    assert np.isclose(p.get("star.distance"), (10.5*u.lyr).to_value(u.pc))
    assert tomllib.loads(open(in_tmp/"run_used.toml").read())["star"]["radius"] == 0.73


def test_include_reads_shared_files_relative_to_the_includer(in_tmp):
    """include = [...] reads the files first (relative to the including file); origins name them."""
    os.makedirs(in_tmp/"shared")
    _write(in_tmp/"shared"/"star.toml", "[star]\ndistance = 3.22\nradius = 0.73\n")
    _write(in_tmp/"shared"/"wind.toml", '[[source]]\nkind = "wind"\nmdot = 6e-13\n')
    path = _write(in_tmp/"run.toml", 'include = ["shared/star.toml", "shared/wind.toml"]\n[grid]\nsampling = 4\nfov = 2.0\n'
                                    '[[source]]\nkind = "wind"\nmdot = 1e-13\nname = "slow"\n')
    p = P.load(path)
    assert p.origin("star.distance") == "set in star.toml" and p.origin("grid.sampling") == "set in run.toml"
    assert [c.name for c in p.components] == ["wind", "slow"]  # [[source]] lists are concatenated


def test_include_key_set_in_two_files(in_tmp):
    """A [section] key set in an included file and the including one is an error naming both."""
    _write(in_tmp/"star.toml", "[star]\ndistance = 3.22\n")
    path = _write(in_tmp/"run.toml", 'include = "star.toml"\n[star]\ndistance = 3.3\n')
    text = _error(P.load, path)
    assert "[star] distance is set in both" in text and "star.toml" in text and "run.toml" in text


def test_include_loop(in_tmp):
    """Files including each other: 'include loop'."""
    _write(in_tmp/"a.toml", 'include = ["b.toml"]\n')
    _write(in_tmp/"b.toml", 'include = ["a.toml"]\n')
    assert "include loop" in _error(P.load, str(in_tmp/"a.toml"))


def test_include_missing_file(in_tmp):
    """An included file that does not exist: 'file not found'."""
    path = _write(in_tmp/"run.toml", 'include = ["nope.toml"]\n')
    assert "include = 'nope.toml': file not found" in _error(P.load, path)


def test_invalid_toml_is_a_param_error(in_tmp):
    """A file that is not valid TOML raises ParamError naming the file."""
    path = _write(in_tmp/"bad.toml", "[grid\nsize = 3\n")
    assert "bad.toml: not valid TOML" in _error(P.load, path)


# ====================================================================== the old layout

@pytest.mark.parametrize("data, legacy", [
    ({"target": {}}, True), ({"telescope": {}}, True), ({"resolved": {}}, True),
    ({"coronagraph": {"charge": 2}}, True), ({"coronagraph": [{"charge": 2}]}, False), ({"grid": {}}, False),
])
def test_is_legacy(data, legacy):
    """The old layout: [target], [telescope], [resolved], or a single [coronagraph] table."""
    assert P.is_legacy(data) is legacy


def test_legacy_dict_is_mapped_with_a_note():
    """An old-layout dict is mapped to the new layout (values kept) with a note about migrate."""
    old = {"target": {"fits_filename": "", "target_pixel": 64}, "telescope": {"diameter": 6.5, "visual_range": 24},
           "coronagraph": {"charge": 4}, "method": {"error_target": 0.01}}
    p = P.from_data(old)
    assert [c.source.kind for c in p.components] == ["fits_image"] and p.components[0].source.settings["file"] == "example"
    assert p.get("grid.size") == 64 and isinstance(p.get("grid.fov"), P.LambdaD) and float(p.get("grid.fov")) == 24
    assert p.get("observation.diameter") == 6.5 and p.coronagraphs[0].name == "vortex_c4"
    assert p.get("method.error_target") == 0.01 and p.notes[0].startswith("Old parameter layout")


def _mcfost_like(path):
    """Write a tiny file that looks like an MCFOST cube (5 axes, FLUX_1 in the header)."""
    from astropy.io import fits
    header = fits.Header()
    header["FLUX_1"] = "I"
    fits.PrimaryHDU(np.zeros((4, 1, 1, 5, 5)), header).writeto(path)


def test_legacy_file_detects_an_mcfost_cube_next_to_it(in_tmp):
    """An old-layout file naming an MCFOST cube (relative to the file) becomes an mcfost source."""
    _mcfost_like(str(in_tmp/"cube.fits"))
    path = _write(in_tmp/"old.toml", '[target]\nfits_filename = "cube"\n')
    assert [c.source.kind for c in P.load(path).components] == ["mcfost"]


def test_legacy_dict_in_python_finds_relative_files_in_the_working_directory(in_tmp):
    """An old-layout dict given in Python (path None) looks for a relative FITS file in the working directory."""
    _mcfost_like(str(in_tmp/"cube.fits"))
    assert [c.source.kind for c in P.from_data({"target": {"fits_filename": "cube"}}).components] == ["mcfost"]


# ====================================================================== params_docs: new_config

def test_new_config_holds_only_the_chosen_modules():
    """A wind-only file: [star] distance and radius, no temperature, no [method], no image keys."""
    text = D.new_config(["wind"], ["vortex:6"])
    for absent in ["temperature", "flux_jy", "[method]", "fits", "mcfost", "cme"]:
        assert absent not in text, absent
    assert text.count("[[source]]") == 1 and "distance = 3.22" in text and "radius = 0.73" in text
    assert "mdot = 6e-13" in text and "(required)" in text


def test_new_config_image_only_has_auto_grid_and_no_star():
    """An image-only file: grid "auto" from the image, no [star]."""
    text = D.new_config(["fits_image"], ["vortex:2"])
    assert "[star]" not in text and 'pixel_scale = "auto"' in text and 'size = "auto"' in text


def test_new_config_mixed_families_add_the_star_flux_and_method():
    """Contrast + Jy modules (wind + planet): temperature, radius and distance for the unit conversion, and [method]."""
    text = D.new_config(["wind", "planet"], ["vortex:2"])
    assert "temperature = 5080.0" in text and "unit conversion" in text and 'combine = "separate"' in text


@pytest.mark.parametrize("sources, coronagraphs, full", [
    (["fits_image"], ["vortex:2"], False),
    (["wind", "cme"], ["vortex:2", "vortex:6"], False),
    (["fits_image", "wind", "planet"], [":4"], True),
    ([{"kind": "wind", "mdot": 1e-12}], [{"type": "vortex", "charge": 8}], False),
])
def test_new_config_loads_back(sources, coronagraphs, full):
    """Every generated file is valid TOML that from_data accepts, with the chosen modules and coronagraphs."""
    p = P.from_data(tomllib.loads(D.new_config(sources, coronagraphs, full=full)))
    kinds = [s if isinstance(s, str) else s["kind"] for s in sources]
    assert [c.source.kind for c in p.components] == kinds and len(p.coronagraphs) == len(coronagraphs)


def test_new_config_values_given_replace_defaults_and_names():
    """Values given in a source dict are written instead of the defaults; coronagraph charges set the names."""
    text = D.new_config([{"kind": "wind", "mdot": 1e-12, "name": "slow_wind"}], ["vortex:2", "vortex:6"])
    p = P.from_data(tomllib.loads(text))
    assert p.components[0].source.settings["mdot"] == 1e-12 and p.components[0].name == "slow_wind"
    assert [c.name for c in p.coronagraphs] == ["vortex_c2", "vortex_c6"]


def test_new_config_writes_the_module_line_for_own_modules():
    """Regression: a module of your own gets its module = "..." line, so the file loads back."""
    text = D.new_config([{"kind": "ring", "name": "inner", "module": "extra_sources"}], ["vortex:2"])
    assert 'module = "extra_sources"' in text and 'name = "inner"' in text
    assert [c.name for c in P.from_data(tomllib.loads(text)).components] == ["inner"]


def test_new_config_unknown_kind():
    """An unknown source kind raises ParamError listing the known kinds."""
    assert "unknown source kind 'fits_imag'" in _error(D.new_config, ["fits_imag"], ["vortex:2"])


@pytest.mark.parametrize("coronagraphs", [["lyot:2"], [{"type": "lyot"}]])
def test_new_config_unknown_coronagraph_type(coronagraphs):
    """Regression: an unknown coronagraph type raises ParamError (not KeyError)."""
    assert "unknown coronagraph type 'lyot'; known: vortex" in _error(D.new_config, ["wind"], coronagraphs)


def test_new_config_full_lists_every_section_with_optional_keys_commented():
    """full=True: every core section, optional keys commented out; the file still loads."""
    text = D.new_config(["fits_image", "wind", "planet"], ["vortex:2"], full=True)
    for section in P.SECTIONS:
        assert "[%s]" % section in text, section
    assert "# hdu = 0" in text and "# use_gpu = false" in text
    P.from_data(tomllib.loads(text))


def test_new_config_path_writes_the_file_and_its_schema(in_tmp):
    """path=...: the file and the JSON schema next to it; the header points to the schema."""
    text = D.new_config(["wind"], ["vortex:2"], path=str(in_tmp/"eps.toml"))
    assert open(in_tmp/"eps.toml").read() == text and text.startswith("#:schema ./%s" % P.SCHEMA_FILE)
    assert (in_tmp/P.SCHEMA_FILE).exists()
    assert sorted(os.listdir(in_tmp)) == sorted(["eps.toml", P.SCHEMA_FILE])


def test_new_config_without_header():
    """header=False starts with the first section."""
    assert D.new_config(["wind"], ["vortex:2"], header=False).startswith("[observation]")


# ====================================================================== params_docs: json_schema

@pytest.fixture(scope="module")
def validator():
    """A Draft-07 validator of json_schema() (the test is skipped without jsonschema)."""
    jsonschema = pytest.importorskip("jsonschema")
    schema = D.json_schema()
    jsonschema.Draft7Validator.check_schema(schema)
    return jsonschema.Draft7Validator(schema)


def _schema_errors(validator, text):
    """The paths ("a/0/b") of the schema errors of a TOML text, sorted."""
    return sorted("/".join(map(str, e.path)) for e in validator.iter_errors(tomllib.loads(text)))


@pytest.mark.parametrize("sources", [["wind", "cme"], ["fits_image"], ["fits_image", "wind", "planet"]])
@pytest.mark.parametrize("full", [False, True])
def test_json_schema_accepts_generated_files(validator, sources, full):
    """The schema accepts every generated file."""
    assert _schema_errors(validator, D.new_config(sources, ["vortex:2", "vortex:6"], full=full)) == []


def test_json_schema_rejects_bad_values(validator):
    """Unknown keys, wrong types, values out of range and shared keys in a module table are flagged."""
    bad = ('[grid]\nsize = "big"\nsizes = 3\n[[source]]\nkind = "wind"\nmdot = -1\ndistance = 3.2\n'
           '[[coronagraph]]\ntype = "vortex"\ncharge = "two"\n')
    assert {"grid", "grid/size", "source/0", "source/0/mdot", "coronagraph/0/charge"} <= set(_schema_errors(validator, bad))


@pytest.mark.parametrize("text, errors", [
    ('[[source]]\nkind = "mine"\n', ["source/0"]),                   # an unknown kind needs module
    ('[[source]]\nkind = "mine"\nmodule = "my.mod"\n', []),
    ('[star]\ndistance = "3.22 pc"\n[grid]\nfov = "32 lambda/D"\n', []),  # unit strings accepted
    ('[star]\ndistance = 0\n', ["star/distance"]),                   # exclusive minimum
    ('[[coronagraph]]\ntype = "lyot"\n', ["coronagraph/0/type"]),
])
def test_json_schema_cases(validator, text, errors):
    """Own modules need `module`; values with units pass; bounds and coronagraph types are checked."""
    assert _schema_errors(validator, text) == errors


# ====================================================================== params_docs: describe, reference, migrate

@pytest.mark.parametrize("query, words", [
    ("star.distance", "needed by: wind, cme"),
    ("star.distance", "if values disagree (> 0.5%): error"),
    ("[star] ", "[star]: the star / system"),
    ("[grid]", "first image's own grid"),
    ("grid.resample", "resolution study"),
    ("wind", "needs: star.distance, star.radius"),
    ("wind.mdot", "[[source]] (kind = \"wind\") mdot  (number [M_sun/yr])"),
    ("wind.scale", "multiply this component by a factor"),
    ("charge", "[[coronagraph]] (type = \"vortex\") charge"),
    ("vortex", "a vortex coronagraph with a Lyot stop"),
    ("mdt", "Did you mean: mdot?"),
])
def test_describe(query, words):
    """describe() explains keys, sections, modules and coronagraphs, and suggests close names."""
    assert words in D.describe(query)


def test_describe_nothing_matches():
    """An unknown name with no close match: 'Nothing called ...'."""
    assert D.describe("zzzzzz") == "Nothing called 'zzzzzz'."


def test_reference_markdown_covers_every_section_kind_and_rule():
    """reference_markdown(): every section, every built-in kind, the coronagraphs and the conflict rules."""
    text = D.reference_markdown()
    for name in P.SECTIONS:
        assert "## [%s]" % name in text
    for kind in ("fits_image", "mcfost", "scene_fits", "wind", "cme", "planet"):
        assert '### kind = "%s"' % kind in text
    assert '### type = "vortex"' in text and "| `star.distance` | 0.5% | error |" in text


OLD_TOML = ('[target]\nfits_filename = ""\ntarget_pixel = 64\n[telescope]\ndiameter = 6.5\nvisual_range = 24\n'
            '[coronagraph]\ncharge = 4\n[method]\nerror_target = 0.01\n[resolved]\npsf_scale = 0.1\n')


def test_migrate_old_toml(in_tmp):
    """migrate(): the new layout keeps the values set, drops [resolved], writes the schema; the old file is unchanged."""
    old = _write(in_tmp/"old.toml", OLD_TOML)
    text = D.migrate(old, str(in_tmp/"new.toml"))
    assert open(old).read() == OLD_TOML and (in_tmp/P.SCHEMA_FILE).exists() and "resolved" not in text
    p = P.load(str(in_tmp/"new.toml"))
    assert p.get("grid.size") == 64 and float(p.get("grid.fov")) == 24 and isinstance(p.get("grid.fov"), P.LambdaD)
    assert p.get("observation.diameter") == 6.5 and p.coronagraphs[0].name == "vortex_c4"
    assert p.get("method.error_target") == 0.01 and p.components[0].source.settings["file"] == "example"


def test_migrate_refuses_a_new_layout_file(in_tmp):
    """migrate() of a file already in the new layout raises ParamError."""
    path = _write(in_tmp/"new.toml", D.new_config(["wind"], ["vortex:2"], schema=False))
    assert "already in the new layout" in _error(D.migrate, path)


def test_migrate_old_python_file(in_tmp):
    """migrate() of an old .py parameter file gives the same settings in the new layout."""
    path = _write(in_tmp/"toycoronagraph_para.py", "target_pixel = 128\nwavelength = 1.6e-6\nD = 6.5\nvisual_range = 24\n"
                                                   "lyot_mask_size = 0.9\nnum_cores = 2\n")
    p = P.from_data(tomllib.loads(D.migrate(path)))
    assert (p.get("grid.size"), float(p.get("grid.fov")), p.get("observation.diameter"), p.get("observation.wavelength"),
            p.coronagraphs[0].values["lyot_mask_size"], p.get("compute.num_cores")) == (128, 24.0, 6.5, 1.6e-6, 0.9, 2)


# ====================================================================== the command line

def test_cli_describe_prints(capsys):
    """python -m toycoronagraph describe KEY prints describe(KEY)."""
    cli(["describe", "star.distance"])
    assert capsys.readouterr().out.strip() == D.describe("star.distance")


def test_cli_new_and_schema_write_files(in_tmp, capsys):
    """`new` writes a file that loads (and its schema); `schema -o` writes the schema."""
    cli(["new", "--sources", "wind,planet", "--coronagraphs", "vortex:2,vortex:6", "-o", "x.toml"])
    cli(["schema", "-o", "s.json"])
    p = P.load("x.toml")
    assert [c.source.kind for c in p.components] == ["wind", "planet"] and (in_tmp/"s.json").exists()
    assert (in_tmp/P.SCHEMA_FILE).exists() and "wrote x.toml" in capsys.readouterr().out


def test_cli_new_does_not_overwrite_without_force(in_tmp):
    """`new` refuses to overwrite an existing file unless --force is given."""
    _write(in_tmp/"mine.toml", "x = 1\n")
    with pytest.raises(SystemExit, match="use --force"):
        cli(["new", "-o", "mine.toml"])
    assert open(in_tmp/"mine.toml").read() == "x = 1\n"
