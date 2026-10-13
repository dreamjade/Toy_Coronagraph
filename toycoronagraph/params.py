"""
Parameters of a scene run: declared once, read from TOML, checked, resolved with provenance.

Every quantity has ONE home:
    [observation]    wavelength, diameter                       (the observation)
    [grid]           size, pixel_scale | sampling, fov, resample (the simulated image)
    [star]           distance, radius, temperature | flux_jy     (the star / system)
    [[source]]       kind + that module's own parameters only (input_* keys describe its file)
    [[coronagraph]]  type + the optics' parameters
    [method] [compute] [checks] [run] [display]
Modules declare the shared quantities they need (`requires`) and those they can supply
(`provides`, e.g. a FITS header's pixel scale); a shared key inside a [[source]] table is
refused with a pointer to its home. Values are resolved user > provided > default, every
disagreement is checked against a per-quantity rule (error / warning / resample), and the
origin of each value is kept for the printed table and the params_used.toml file.

The steps:
    load(path)        reads the file and its `include = [...]` files (_read_toml), then from_data()
    from_data(data)   checks the layout ([[source]] / [[coronagraph]] lists, known sections), builds a
                      Component per [[source]] table, then make()
    make(...)         checks the core sections and the coronagraphs, names the components, and returns
                      a Params after Params.check() (what is missing / unused)
    Params.resolve()  the final shared values (user > provided by the sources > default), with the
                      conflict rules; Params.table() / used_toml() show them with their origins
Every function collects all the problems it finds and raises one ParamError listing them, one per line.

Everything users read is generated from these declarations by toycoronagraph.params_docs
(example files, the JSON schema, the reference tables, describe, migrate). The old layout
([target], [telescope], a single [coronagraph]) is read and mapped to the new one (legacy_to_new).
"""
import os
import re
import json
import difflib
import datetime
from dataclasses import dataclass
from types import SimpleNamespace
try:
    import tomllib # Python >= 3.11
except ModuleNotFoundError: # pragma: no cover
    import tomli as tomllib

AUTO = "auto"
RAD_TO_ARCSEC = 206264.806247
SCHEMA_FILE = "toycoronagraph.schema.json"


# ====================================================================== declarations

@dataclass(frozen=True)
class Param:
    """
    The declaration of one parameter: how it is checked (convert()) and documented (params_docs).

    Attributes:
        name (str): key in the TOML file.
        default (object): value when not set; None = no default (required when it is needed).
        unit (str): canonical unit (astropy syntax), or None; numbers are in this unit, strings with a
            unit ("3.22 pc") are converted to it.
        doc (str): short description (the comment in generated files).
        type (str): "float", "int", "bool", "str", "int_list", "float_list", "str_list" or "table"
            (ignored when choices are given).
        choices (tuple): allowed values, if limited (compared with their type: 1 is not true).
        min, max (float): bounds of a number, or None; min is exclusive if exclusive_min.
        exclusive_min (bool): the value must be > min (else >= min).
        auto (bool): "auto" is also allowed.
        auto_word (str): another allowed word (e.g. "all"), or None.
        lambda_d (bool): "32 lambda/D" is also allowed (read as a LambdaD, converted to arcsec when the
            wavelength and diameter are known).
        advanced (bool): commented out in generated files (or left out unless full=True).
        example (object): value written in generated files when there is no default.
        help (str): longer description (editor hover, describe()).
    """
    name: str
    default: object = None
    unit: str = None
    doc: str = ""
    type: str = "float"
    choices: tuple = ()
    min: float = None
    max: float = None
    exclusive_min: bool = False
    auto: bool = False          # "auto" allowed
    auto_word: str = None       # another allowed word (e.g. "all")
    lambda_d: bool = False      # "32 lambda/D" allowed (converted when wavelength and diameter are known)
    advanced: bool = False      # commented out in generated files
    example: object = None      # value written in generated files when there is no default
    help: str = ""

    def describe_type(self):
        """
        The allowed values in words, as shown in error messages and the docs.

        Returns:
            str: e.g. '"auto" or number [m]', 'one of "separate", "together"', 'number [arcsec] or "N lambda/D"'.
        """
        if self.choices:
            return "one of " + ", ".join(_toml(c) for c in self.choices)
        words = {"float": "number", "int": "integer", "bool": "true/false", "str": "text", "table": "table",
                 "int_list": "list of integers", "float_list": "list of numbers", "str_list": "list of text"}
        text = words[self.type]
        if self.auto:
            text = '"auto" or ' + text
        if self.auto_word:
            text = '"%s" or %s' % (self.auto_word, text)
        if self.unit:
            text += " [%s]" % self.unit
        if self.lambda_d:
            text += ' or "N lambda/D"'
        return text


class LambdaD(float):
    """A size given in lambda/D (resolved to arcsec when the wavelength and diameter are known)."""
    def __repr__(self):
        """The TOML form, e.g. '"32 lambda/D"'."""
        return '"%g lambda/D"' % float(self)


IMAGE_UNITS = ("auto", "Jy/arcsec^2", "Jy/pixel", "contrast_pixel", "contrast_peak")

# Core sections: name -> (doc, params)
SECTIONS = {
    "observation": ("the observation", (
        Param("wavelength", AUTO, unit="m", doc='observing wavelength; "auto" = from the image header', auto=True,
              min=0, exclusive_min=True, example=1e-6,
              help="Used by the PSFs of every coronagraph and by the Jy conversions. \"auto\" takes it from the "
                   "first image source's header (WAVELEN in m, or MCFOST's WAVE in um)."),
        Param("diameter", 2.4, unit="m", doc="telescope aperture diameter", min=0, exclusive_min=True),
    )),
    "grid": ("the simulated image, star at pixel (N-1)/2. By default the first image's own grid; modules that "
             "compute their image draw on it, other images are resampled onto it", (
        Param("size", AUTO, doc='pixels on a side; "auto" = the first image source', type="int", auto=True, min=3,
              help="Odd sizes put a pixel exactly on the star (recommended for winds and CMEs)."),
        Param("pixel_scale", AUTO, unit="arcsec", doc='arcsec per pixel; "auto" = the first image source',
              auto=True, min=0, exclusive_min=True,
              help="Set this OR sampling. Image sources with another pixel scale are resampled onto the grid."),
        Param("sampling", None, doc="pixels per lambda/D (set this OR pixel_scale)", min=0, exclusive_min=True,
              example=4.0),
        Param("fov", None, unit="arcsec", doc='image width, e.g. 2.0 or "32 lambda/D"', min=0, exclusive_min=True,
              lambda_d=True, example=LambdaD(32),
              help="With pixel_scale or sampling it sets the size (made odd). If size is also set, the three must agree."),
        Param("resample", True, doc="fit images covering another field onto the grid; false: same field only",
              type="bool",
              help="Images read from files (fits_image, mcfost, scene_fits; a cme model given its own dim/extent) are put "
                   "on the grid, conserving the flux (a note says so). An image covering the grid's field of view "
                   "(pixel_scale x size) at another resolution is always rebinned. Other images are cropped, padded "
                   "or rescaled with true; with false they stop the run with an error listing each one, so no "
                   "image's field is changed. For a resolution study, set only [grid] size: the image's field is "
                   "kept and the pixel scale follows. Differences below 0.01 pixel at the edge (e.g. a rounded "
                   "header value) are ignored."),
    )),
    "star": ("the star / system", (
        Param("distance", None, unit="pc", doc="distance to the star", min=0, exclusive_min=True, example=3.22),
        Param("radius", None, unit="R_sun", doc="stellar radius", min=0, exclusive_min=True, example=0.73),
        Param("temperature", None, unit="K", doc="effective temperature (blackbody; set this OR flux_jy)", min=0,
              exclusive_min=True, example=5080.0),
        Param("flux_jy", None, unit="Jy", doc="flux density at the wavelength (set this OR temperature)", min=0,
              exclusive_min=True, example=100.0),
    )),
    "method": ("how the images go through the coronagraphs", (
        Param("combine", "separate", doc='"separate": each component its own pass (saved); "together": one pass',
              choices=("separate", "together"),
              help="\"separate\" (default) sends each component through each coronagraph and sums them, so each "
                   "contribution is saved and each component can use its best method. \"together\" sums the "
                   "components first: one pass, faster, no per-component output."),
        Param("method", AUTO, doc='"auto", "circular" or "spatial"', choices=(AUTO, "circular", "spatial")),
        Param("spatial_mode", AUTO, doc='"auto", "rotate" or "direct"', choices=(AUTO, "rotate", "direct"),
              advanced=True),
        Param("error_target", 0.02, doc="target relative error of the spatial method", min=0, exclusive_min=True),
        Param("error_metric", "total", doc='"total" (of all light) or "local" (per annulus; for high contrast)',
              choices=("total", "local"), advanced=True),
        Param("symmetry_tol", 0.01, doc="asymmetry above which a component counts as not symmetric", min=0,
              advanced=True),
        Param("rot_number", 0, doc="circular method: number of rotations; 0 = grid size", type="int", min=0,
              advanced=True),
    )),
    "compute": ("computing", (
        Param("num_cores", -1, doc="CPU cores; -1 = all", type="int"),
        Param("use_gpu", False, doc="use CuPy for the spatial convolutions", type="bool", advanced=True),
        Param("max_kernel_gb", 4.0, doc="largest spatial PSF cube to build, in GB", min=0, exclusive_min=True,
              advanced=True),
        Param("psf_cache_dir", "", doc='folder for computed PSFs; "" = the working directory', type="str",
              advanced=True),
    )),
    "checks": ("consistency checks", (
        Param("on_conflict", "error", doc='"error", or "warn" to accept disagreeing values', choices=("error", "warn")),
    )),
    "run": ("the run and its saved output", (
        Param("name", AUTO, doc='name of the results folder; "auto" = the parameter file\'s name', type="str"),
        Param("output_dir", "results", doc="where run folders are written", type="str"),
        Param("overwrite", False, doc="overwrite an existing run folder (else _2, _3, ...)", type="bool"),
        Param("write_back", False, doc="also append the values used to this file (old behaviour)", type="bool",
              advanced=True),
    )),
    "display": ("plots", (
        Param("unit", AUTO, doc='"auto", "Jy/arcsec^2", "Jy/pixel", "contrast_pixel" or "contrast_peak"',
              choices=IMAGE_UNITS,
              help="Only changes how results are shown and exported. contrast_peak = surface brightness / the "
                   "peak of the star's PSF through the same optics without the coronagraph's mask (the Lyot stop "
                   "kept; per coronagraph); input images: the ideal aperture's peak."),
        Param("log", True, doc="logarithmic scales: images' colours, profiles' axes", type="bool"),
        Param("log_range", 5.0, doc="log scales: decades shown below the top (profiles: 3 more)", min=0,
              exclusive_min=True),
        Param("axis", "arcsec", doc='plot axes in "arcsec", "mas", "lambda/D", "au" (needs [star] distance) or "pixel"',
              choices=("arcsec", "mas", "lambda/D", "au", "pixel")),
        Param("dpi", 150, doc="resolution of saved figures, dots per inch", type="int", min=1),
        Param("label", "", doc='colour bar label; "" = the unit', type="str", advanced=True),
    )),
}
# Sections whose values are resolved (Params.resolve) and saved with their origins
SHARED_SECTIONS = ("observation", "grid", "star")
# Keys allowed at the top of a file, outside any section (not used by the code: _read_toml handles include)
TOP_LEVEL = ("include",)

# Keys every [[source]] table may have
SOURCE_COMMON = (
    Param("kind", None, doc="which module", type="str"),
    Param("name", AUTO, doc='component name (FITS extension, plots); "auto" = the kind', type="str"),
    Param("enabled", True, doc="false switches this component off", type="bool"),
    Param("scale", 1.0, doc="multiply this component by a factor", min=0),
    Param("method", AUTO, doc='override for this component: "circular" or "spatial"', choices=(AUTO, "circular", "spatial"),
          advanced=True),
    Param("module", "", doc="Python module defining this kind (your own modules)", type="str", advanced=True),
)

# Coronagraph types: type -> (doc, params); future backends add theirs with register_coronagraph (and an
# Optics builder in toycoronagraph.observe.OPTICS_BUILDERS)
CORONAGRAPHS = {
    "vortex": ("a vortex coronagraph with a Lyot stop", (
        Param("charge", 2, doc="vortex charge (positive, even)", type="int", min=2),
        Param("lyot_mask_size", 0.95, doc="Lyot stop diameter relative to the aperture", min=0, max=1,
              exclusive_min=True),
        Param("symmetric", AUTO, doc="circularly symmetric PSFs? \"auto\" = from the optics model",
              choices=(AUTO, True, False), advanced=True),
    )),
    "lazuli": ("the Lazuli space telescope (UASAL), alone or with a vortex coronagraph, computed with prysm", (
        Param("charge", 0, doc="0 = the telescope alone; else the vortex charge (positive, even)", type="int", min=0,
              help="With a charge, an ideal vortex focal-plane mask of that charge and a circular Lyot stop "
                   "(lyot_mask_size) follow the telescope. Needs prysm and the Lazuli config packages "
                   "(installed from git, see the README)."),
        Param("lyot_mask_size", 0.95, doc="Lyot stop diameter relative to the aperture (used with a charge)", min=0,
              max=1, exclusive_min=True),
        Param("wavelengths", AUTO, unit="m", doc='wavelengths of the PSF, e.g. ["600 nm", "650 nm"]; "auto" = '
              "[observation] wavelength", type="float_list", auto=True, min=0, exclusive_min=True,
              help="The PSF is the weighted sum of one PSF per wavelength (weights: throughput). Numbers are in m, "
                   "text with a unit is converted. Positions and the image grid stay in units of [observation] "
                   "wavelength / diameter."),
        Param("throughput", AUTO, doc='relative weight of each wavelength; "auto" = the Lazuli mirrors, r\' filter '
              "and detector QE", type="float_list", auto=True, min=0,
              help="One weight per wavelength (normalised to sum 1). \"auto\" multiplies the M1 aluminium and M2-M4 "
                   "silver reflectances, the SDSS r' filter and the WCC detector QE of config_um / config_um_wcc, "
                   "times the band each wavelength stands for (halfway to its neighbours), so the PSF is the band "
                   "integral however the wavelengths are spaced."),
        Param("opd", "", doc='wavefront error map (FITS or .npy; nm unless its BUNIT says otherwise); "" = ideal optics',
              type="str",
              help="Lazuli's aberrations. Without a map (\"\") the optics are ideal: a perfect telescope and, with "
                   "a charge, an ideal coronagraph. A square map over the aperture diameter; its size sets "
                   "pupil_samples. It makes the PSFs not circularly symmetric, so images use the spatial method in "
                   "direct mode."),
        Param("pupil_samples", AUTO, doc='pupil samples across the aperture; "auto" = the OPD map\'s size, else '
              "config_um_wcc's beam_sampling", type="int", auto=True, min=16, advanced=True,
              help="The focal field repeats every pupil_samples lambda/D (Lazuli's, at the shortest wavelength), "
                   "so each source would reappear that far away: the run stops with the smallest sampling that "
                   "clears the PSF frames when it is too small. Above that, fewer samples are faster and less "
                   "accurate (see the README)."),
        Param("symmetric", AUTO, doc="circularly symmetric PSFs? \"auto\" = from the optics model (not with an opd)",
              choices=(AUTO, True, False), advanced=True),
    )),
}
# Keys every [[coronagraph]] table may have
CORONAGRAPH_COMMON = (
    Param("type", "vortex", doc="coronagraph model", type="str"),
    Param("name", AUTO, doc='name in outputs; "auto" = e.g. vortex_c2', type="str"),
)

# Conflict rules: quantity -> (relative tolerance, consequence, why)
CONFLICT_RULES = {
    "star.distance": (0.005, "error", "it changes the angular size and brightness of everything"),
    "star.radius": (0.005, "error", "it changes where the wind starts and the stellar flux"),
    "star.temperature": (0.005, "error", "it changes the stellar flux"),
    "observation.wavelength": (0.01, "warning", "the image's brightness is at its own wavelength, the PSFs at the scene's"),
}
# Relative tolerance of over-specified settings ([grid] pixel_scale + sampling, size x pixel_scale vs fov
# (at least 1.5 pixels there), [star] temperature + flux_jy)
GROUP_TOL = 0.001


def register_coronagraph(type_name, doc, params):
    """
    Declare the parameters of a new coronagraph type (with its optics backend).

    Args:
        type_name (str): the value of `type` in [[coronagraph]] tables.
        doc (str): one line describing it.
        params (iterable of Param): its own parameters (besides CORONAGRAPH_COMMON's type and name).
    """
    CORONAGRAPHS[type_name] = (doc, tuple(params))


def coronagraph_problems(ctype, values, refused=()):
    """
    The checks of a coronagraph table that involve several keys (each value is checked on its own by convert()).

    Args:
        ctype (str): the coronagraph type.
        values (dict): the checked values set in the table; unset keys are absent (their defaults are valid).
        refused (iterable of str): keys set in the table but refused by their own check (so absent from values,
            where they would read as their defaults): the checks involving them are left out.

    Returns:
        list of str: the problems, e.g. "charge = 3: a vortex charge must be even"; [] if none.
    """
    problems = []
    if ctype in ("vortex", "lazuli") and values.get("charge", 0) % 2:
        problems.append("charge = %d: a vortex charge must be even" % values["charge"])
    if ctype == "lazuli" and not {"wavelengths", "throughput"} & set(refused):
        wavelengths, weights = values.get("wavelengths", AUTO), values.get("throughput", AUTO)
        n = 1 if wavelengths == AUTO else len(wavelengths)
        if n == 0:
            problems.append('wavelengths = []: give at least one wavelength, or "auto"')
        elif weights != AUTO and len(weights) != n:
            problems.append('throughput has %d values for %d wavelength(s): give one per wavelength, or "auto"'
                            % (len(weights), n))
        elif weights != AUTO and not sum(weights) > 0:
            problems.append("throughput: the weights are all 0")
    return problems


def _sources():
    """
    The toycoronagraph.sources module, imported on first use (it imports this module).

    Returns:
        module: toycoronagraph.sources.
    """
    from toycoronagraph import sources
    return sources


def _source_class(kind, module=""):
    """
    The Source class of a kind.

    Args:
        kind (str): the [[source]] kind.
        module (str): Python module to import first if the kind is not registered yet; "" = none.

    Returns:
        type or None: the Source subclass, or None if the kind is unknown.

    Raises:
        ValueError: module is set but cannot be imported.
    """
    try:
        return _sources().get_source(kind, module or None)
    except KeyError:
        return None
    except ImportError as err:
        raise ValueError("module = %r for kind %r could not be imported: %s" % (module, kind, err))


# ====================================================================== values

# "32 lambda/D" (also lam/D, λ/D, any case); group 1 is the number
_LAMBDA_D = re.compile(r"^\s*([-+0-9.eE]+)\s*(lambda|lam|λ)\s*/\s*d\s*$", re.I)
# A number followed by a unit, e.g. "3.22 pc", "1e-6 m" (also the JSON schema's pattern)
_QUANTITY = re.compile(r"^\s*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?\s*[A-Za-z_(].*$")


def _toml(value):
    """
    A value written as TOML (inline).

    Args:
        value (object): a LambdaD, bool, str, dict, list / tuple, float or int (others: their repr).

    Returns:
        str: e.g. 'true', '"text"', '1.0', '[1, 2]', '{a = 1}', '"32 lambda/D"'.
    """
    if isinstance(value, LambdaD):
        return '"%g lambda/D"' % float(value)
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, dict):
        return "{%s}" % ", ".join("%s = %s" % (k, _toml(v)) for k, v in value.items())
    if isinstance(value, (list, tuple)):
        return "[%s]" % ", ".join(_toml(v) for v in value)
    if isinstance(value, float):
        text = repr(value)
        # a float must look like one: "1.0", not "1" ("1e+20", "nan" and "inf" already do)
        return text if any(c in text for c in ".en") else text+".0"
    return repr(value)


def _is_number(v):
    """
    Whether v is an int or float (a bool is not a number here).

    Args:
        v (object): any value.

    Returns:
        bool
    """
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def convert(param, value, where):
    """
    Check one value against its declaration and return it in canonical form.

    In order: "auto" (if allowed), the auto_word, the choices, then the type. A float (and each item
    of a float_list) may be a number, "N lambda/D" (if lambda_d; returned as a LambdaD, which converts
    to itself) or a string with a unit (converted to param.unit with astropy; "Msun" / "Rsun" are
    accepted). Numbers are then checked against min / max.

    Args:
        param (Param): the declaration.
        value (object): the value as read from TOML (or given in Python).
        where (str): where it is set, starting the message (e.g. "[grid]", "[[source]] #1 (wind)").

    Returns:
        object: "auto", the auto_word, the matching choice, a bool / str / dict, a list (float_list items
            as floats in param.unit), an int, a float in param.unit, or a LambdaD.

    Raises:
        ValueError: '<where> <name> = <value>: <why> (expected <param.describe_type()>)'.
    """
    p = param

    def fail(why):
        """The ValueError for this value: why it is refused and what is expected."""
        return ValueError("%s %s = %s: %s (expected %s)" % (where, p.name, _toml(value), why, p.describe_type()))

    def bounded(number, item=""):
        """
        A number checked against p.min / p.max.

        Args:
            number (int or float): the number.
            item (str): names a list item at the start of the message ("" for the value itself).

        Returns:
            int or float: number.

        Raises:
            ValueError: out of bounds (from fail()).
        """
        if p.min is not None and (number < p.min or (p.exclusive_min and number == p.min)):
            raise fail("%smust be %s %g" % (item, ">" if p.exclusive_min else ">=", p.min))
        if p.max is not None and number > p.max:
            raise fail("%smust be <= %g" % (item, p.max))
        return number

    def real(v, item=""):
        """
        A float value: a number, "N lambda/D" or a LambdaD (if p.lambda_d) or a string with a unit (if p.unit).

        Args:
            v (object): the value, or one item of a float_list.
            item (str): names a list item at the start of the message ("" for the value itself).

        Returns:
            float: in p.unit, within the bounds; a LambdaD (in lambda/D) for "N lambda/D" or a LambdaD.

        Raises:
            ValueError: not a number, not convertible to p.unit, or out of bounds (from fail()).
        """
        if isinstance(v, LambdaD) and p.lambda_d:  # converted already (make() checks checked sections again)
            number = v
        elif _is_number(v):
            number = float(v)
        elif isinstance(v, str) and p.lambda_d and _LAMBDA_D.match(v):
            number = LambdaD(float(_LAMBDA_D.match(v).group(1)))
        elif isinstance(v, str) and p.unit and _QUANTITY.match(v):
            from astropy import units as u
            try:
                number = float(u.Quantity(v.replace("Msun", "M_sun").replace("Rsun", "R_sun")).to_value(p.unit))
            except Exception as err:
                raise fail("%scannot convert to %s (%s)" % (item, p.unit, err))
        else:
            raise fail(item + "not a number")
        return bounded(number, item)

    if value == AUTO and (p.auto or AUTO in p.choices):
        return AUTO
    if p.auto_word and value == p.auto_word:
        return value
    if p.choices:
        for c in p.choices:
            if value == c and isinstance(value, bool) == isinstance(c, bool):  # 1 == True, but 1 is not a choice true
                return c
        close = difflib.get_close_matches(str(value), [str(c) for c in p.choices], 1)
        raise fail("not allowed" + (' (did you mean "%s"?)' % close[0] if close else ""))
    if p.type == "bool":
        if isinstance(value, bool):
            return value
        raise fail("not true/false")
    if p.type == "str":
        if isinstance(value, str):
            return value
        raise fail("not text")
    if p.type == "table":
        if isinstance(value, dict):
            return value
        raise fail("not a table")
    if p.type == "float_list":
        if isinstance(value, list):
            return [real(v, "item %s: " % _toml(v)) for v in value]
        raise fail("not a list")
    if p.type.endswith("_list"):
        item = {"int_list": lambda v: isinstance(v, int) and not isinstance(v, bool),
                "str_list": lambda v: isinstance(v, str)}[p.type]
        if isinstance(value, list) and all(item(v) for v in value):
            return list(value)
        raise fail("not a list of the right type")
    if p.type == "int":
        if isinstance(value, int) and not isinstance(value, bool):
            return bounded(value)
        raise fail("not an integer")
    return real(value) # float, with an optional unit


# ====================================================================== reading

class ParamError(ValueError):
    """Problems in the parameters (all of them, one per line)."""


def _read_toml(path, seen=()):
    """
    Parse a parameter file and the files it includes (`include = "a.toml"` or `include = [...]`).

    The included files are read first (in order, recursively; relative paths are relative to the
    including file), then the file's own values. [[...]] lists are concatenated; a key of a [section]
    set in two files is an error.

    Args:
        path (str): the file.
        seen (tuple of str): absolute paths of the including files (to detect include loops).

    Returns:
        tuple: (data, origins): data (dict) is the merged file content, without `include`; origins (dict)
            is {(section, key): absolute path of the file setting it}, for the keys of [section] tables.

    Raises:
        ParamError: include loop, included file not found, invalid TOML, or a key set in two files.
    """
    path = os.path.abspath(path)
    if path in seen:
        raise ParamError("include loop: %s" % " -> ".join(list(seen)+[path]))
    with open(path, "rb") as f:
        try:
            data = tomllib.load(f)
        except tomllib.TOMLDecodeError as err:
            raise ParamError("%s: not valid TOML: %s" % (path, err))
    merged, origins = {}, {}
    includes = data.pop("include", [])
    if isinstance(includes, str):
        includes = [includes]
    for inc in includes:
        inc_path = inc if os.path.isabs(inc) else os.path.join(os.path.dirname(path), inc)
        if not os.path.exists(inc_path):
            raise ParamError("%s: include = %r: file not found" % (path, inc))
        sub, sub_origins = _read_toml(inc_path, seen+(path,))
        _merge(merged, origins, sub, sub_origins)
    own_origins = {(section, key): path for section, value in data.items() if isinstance(value, dict) for key in value}
    _merge(merged, origins, data, own_origins)
    return merged, origins


def _merge(merged, origins, new, new_origins):
    """
    Merge one file's content into what has been read so far (in place).

    Lists ([[source]], [[coronagraph]]) are concatenated, tables are merged key by key, other
    top-level values are replaced.

    Args:
        merged (dict): the content read so far; updated.
        origins (dict): {(section, key): file path} of merged; updated.
        new (dict): the content to add.
        new_origins (dict): {(section, key): file path} of new.

    Raises:
        ParamError: keys set in both (all of them, one per line).
    """
    conflicts = []
    for section, value in new.items():
        if isinstance(value, list):
            merged.setdefault(section, [])
            merged[section] += value
        elif isinstance(value, dict):
            target = merged.setdefault(section, {})
            for key, v in value.items():
                if key in target:
                    conflicts.append("[%s] %s is set in both %s and %s: keep it in one file"
                                     % (section, key, origins[(section, key)], new_origins[(section, key)]))
                target[key] = v
                origins[(section, key)] = new_origins.get((section, key))
        else:
            merged[section] = value
    if conflicts:
        raise ParamError("\n".join(conflicts))


def is_legacy(data):
    """
    Whether file content is in the old layout: [target], [telescope] or [resolved] (written back by
    old runs), or a single [coronagraph] table instead of a [[coronagraph]] list.

    Args:
        data (dict): the parsed file.

    Returns:
        bool
    """
    return any(k in data for k in ("target", "telescope", "resolved")) or isinstance(data.get("coronagraph"), dict)


def _where_is(key):
    """
    The homes of a key name among the core sections.

    Args:
        key (str): a key name, e.g. "distance".

    Returns:
        list of str: the SECTIONS declaring it, e.g. ["star"]; [] if none.
    """
    return [s for s, (_, ps) in SECTIONS.items() if any(p.name == key for p in ps)]


def _check_table(table, params, where, errors):
    """
    Check and convert the values of one table (the input is not modified).

    Unknown keys are errors: a core section's key ("distance" in a [[source]] table) points to its
    home, other keys get the closest known name or the list of known keys.

    Args:
        table (dict): {key: value} as read.
        params (iterable of Param): the keys allowed in this table.
        where (str): the table's name in messages, e.g. "[grid]" or "[[source]] #1 (wind)".
        errors (list of str): problems found are appended to it.

    Returns:
        dict: {key: converted value} for the valid keys (see convert()).
    """
    known = {p.name: p for p in params}
    out = {}
    for key, value in table.items():
        if key not in known:
            homes = _where_is(key)
            if homes:
                errors.append("%s %s: %s is a shared setting: set it in [%s]" % (where, key, key, homes[0]))
            else:
                close = difflib.get_close_matches(key, list(known), 1, 0.6)
                errors.append("%s: unknown key %s%s" % (where, key, (" (did you mean %s?)" % close[0]) if close else
                                                         "; known: " + ", ".join(known)))
            continue
        try:
            out[key] = convert(known[key], value, where)
        except ValueError as err:
            errors.append(str(err))
    return out


# ====================================================================== building and checking

@dataclass
class Component:
    """
    One source module in a scene, and how it is used.

    Attributes:
        source (Source): the module instance (toycoronagraph.sources), holding its own settings.
        name (str): name in outputs; "auto" = the kind, numbered from the second of a kind (set by make()).
        scale (float): factor multiplying its image.
        enabled (bool): False leaves it out of the run.
        method (str): "auto", "circular" or "spatial" (override of [method] method for this component).
    """
    source: object
    name: str = AUTO
    scale: float = 1.0
    enabled: bool = True
    method: str = AUTO


def load(path):
    """
    Read and check a parameter file (an old-layout file is mapped to the new layout).

    Args:
        path (str): the TOML file (its `include = [...]` files are read too).

    Returns:
        Params: checked values with their origins; Params.resolve() then fills in what the
            sources provide.

    Raises:
        ParamError: every problem found, one per line.
    """
    data, origins = _read_toml(path)
    return from_data(data, path, origins)


def from_data(data, path=None, origins=None):
    """
    Check parameters given as a dict in the file layout (see load()).

    Args:
        data (dict): the file content: {section: {key: value}}, "source" / "coronagraph": lists of tables.
            An old-layout dict (is_legacy) is mapped with legacy_to_new() first (and a note added).
        path (str): the file, for messages, origins and relative paths in the old layout; None = given in Python
            (relative paths in the old layout are then in the working directory).
        origins (dict): {(section, key): file path} from _read_toml(); None = all from path.

    Returns:
        Params: see make().

    Raises:
        ParamError: every problem found, one per line (the layout and [[source]] tables first; then those
            of make()).
    """
    notes = []
    if is_legacy(data):
        data, origins = legacy_to_new(data, os.path.dirname(os.path.abspath(path)) if path else os.getcwd()), {}
        notes.append("Old parameter layout: read and mapped to the new one (python -m toycoronagraph migrate %s "
                     "writes it in the new layout)" % (path or "FILE"))
    errors = []
    for section, value in data.items():
        if section in ("source", "coronagraph"):
            if not isinstance(value, list):
                errors.append("use [[%s]] (a list of tables), one per %s" % (section, section))
        elif section not in SECTIONS:
            close = difflib.get_close_matches(section, list(SECTIONS)+["source", "coronagraph", "include"], 1)
            errors.append("unknown section [%s]%s" % (section, (" (did you mean [%s]?)" % close[0]) if close else ""))
        elif not isinstance(value, dict):
            errors.append("[%s] must be a table" % section)
    tables = data.get("source", []) if isinstance(data.get("source"), list) else []
    components = [c for c in (_component(i, t, errors) for i, t in enumerate(tables)) if c]
    if errors:
        raise ParamError("Problems in %s:\n  " % (path or "the parameters") + "\n  ".join(errors))
    coronagraphs = data.get("coronagraph", []) if isinstance(data.get("coronagraph"), list) else []
    return make({k: v for k, v in data.items() if k in SECTIONS}, components, coronagraphs, path, origins, notes)


def _component(i, table, errors):
    """
    A Component from a [[source]] table.

    Args:
        i (int): index of the table (0-based; messages say #i+1).
        table (dict): the [[source]] table as read.
        errors (list of str): problems found are appended to it.

    Returns:
        Component or None: None (with errors) if kind is missing or unknown; otherwise the Component
            (its values may still have errors, then appended to errors).

    Raises:
        ValueError: `module` cannot be imported (see _source_class()).
    """
    kind = table.get("kind")
    where = "[[source]] #%d (%s)" % (i+1, kind or "no kind")
    if not kind:
        errors.append(where + ": kind is missing; one of " + ", ".join(_sources().source_kinds()))
        return None
    cls = _source_class(kind, table.get("module", ""))
    if cls is None:
        close = difflib.get_close_matches(kind, _sources().source_kinds(), 1)
        errors.append("%s: unknown kind %r%s" % (where, kind, (" (did you mean %s?)" % close[0]) if close else
                      "; for your own module also set module = \"package.module\""))
        return None
    values = _check_table(table, SOURCE_COMMON+tuple(cls.params), where, errors)
    common = {p.name for p in SOURCE_COMMON}
    source = cls(**{k: v for k, v in values.items() if k not in common})
    return Component(source, **{k: values[k] for k in ("name", "scale", "enabled", "method") if k in values})


def make(sections, components, coronagraphs=(), path=None, origins=None, notes=()):
    """
    Checked parameters from core sections, Components and [[coronagraph]] tables (used by
    from_data() and by scenes built in Python).

    Components named "auto" are named here (in place): the kind, then kind_2, kind_3, ... by position
    among the components of the same kind. Without coronagraphs, one vortex of charge 2 is used.

    Args:
        sections (dict): {section name: {key: value}} for the core SECTIONS.
        components (list of Component): the scene's components.
        coronagraphs (list of dict): [[coronagraph]] tables, {key: value} with "type" (default "vortex").
        path (str): the parameter file, or None (built in Python).
        origins (dict): {(section, key): file path} from _read_toml(), or None.
        notes (iterable of str): notes from reading (e.g. the old layout).

    Returns:
        Params: after Params.check().

    Raises:
        ParamError: every problem found, one per line (also those of Params.check()).
    """
    errors = []
    checked = {}
    for name, table in sections.items():
        if name not in SECTIONS:
            errors.append("unknown section [%s]" % name)
        else:
            checked[name] = _check_table(table, SECTIONS[name][1], "[%s]" % name, errors)
    for i, c in enumerate(components):  # automatic names: ring, ring_2, ring_3, ... (by position among the same kind)
        # components built in Python were not checked by _component(): check their common settings
        _check_table({"scale": c.scale, "enabled": c.enabled, "method": c.method}, SOURCE_COMMON, "source %s" % c.name, errors)
        if c.name == AUTO:
            k = 1 + sum(d.source.kind == c.source.kind for d in components[:i])
            c.name = c.source.kind if k == 1 else "%s_%d" % (c.source.kind, k)
    names = [c.name for c in components]
    errors += ["two components are named %r: give them different names" % n for n in sorted(set(names)) if names.count(n) > 1]
    errors += ['a component cannot be named %r: it names the sum in plots and saved files' % n
               for n in ("total", "input") if n in names]
    coros = [c for c in (_coronagraph(i, t, errors) for i, t in enumerate(coronagraphs)) if c]
    if len({c.name for c in coros}) != len(coros):
        errors.append("two [[coronagraph]] tables have the same name: give them different names")
    if errors:
        raise ParamError("Problems in %s:\n  " % (path or "the parameters") + "\n  ".join(errors))
    params = Params(path, checked, components, coros or [SimpleNamespace(type="vortex", name="vortex_c2", values={})],
                    origins or {}, notes)
    params.check()
    return params


def _coronagraph(i, table, errors):
    """
    Check a [[coronagraph]] table.

    Args:
        i (int): index of the table (0-based; messages say #i+1).
        table (dict): the table as read; "type" defaults to "vortex".
        errors (list of str): problems found are appended to it.

    Returns:
        SimpleNamespace or None: None (with errors) for an unknown type; otherwise type (str), name (str;
            "auto" = <type>_c<charge> for types with a charge, e.g. vortex_c2 or lazuli_c0, else
            <type>_<i+1>) and values (dict: the checked values set in the table, type and name included;
            unset keys are absent).
    """
    ctype = table.get("type", "vortex")
    where = "[[coronagraph]] #%d (%s)" % (i+1, ctype)
    if ctype not in CORONAGRAPHS:
        errors.append("%s: unknown type; known: %s" % (where, ", ".join(CORONAGRAPHS)))
        return None
    values = _check_table(table, CORONAGRAPH_COMMON+CORONAGRAPHS[ctype][1], where, errors)
    errors += ["%s: %s" % (where, problem) for problem in coronagraph_problems(ctype, values, set(table)-set(values))]
    name = values.get("name", AUTO)
    if name == AUTO:
        charge = values.get("charge", {p.name: p.default for p in CORONAGRAPHS[ctype][1]}.get("charge"))
        name = "%s_%d" % (ctype, i+1) if charge is None else "%s_c%d" % (ctype, charge)
    return SimpleNamespace(type=ctype, name=name, values=values)


class Params:
    """
    Checked parameters: get() / origin() / used_by(); check() for what is missing or unused;
    resolve() for the final shared values, with what the sources provide.

    Keys of the core sections are written "section.name", e.g. "star.distance".

    Attributes:
        path (str): the parameter file, or None (built in Python).
        sections (dict): {section: {key: checked value}}: only what the user set.
        components (list of Component): all components (enabled() for those used).
        coronagraphs (list of SimpleNamespace): type, name, values (see _coronagraph()).
        load_notes (list of str): notes from reading and check().
        notes (list of str): load_notes + the messages of the last resolve().
        provided (list of tuple): (key, value, origin) supplied by the sources (Source.inspect()), set
            before resolve().
        resolved (dict): {key: (value, origin)}: the final shared values, after resolve().
    """

    def __init__(self, path, sections, components, coronagraphs, origins, notes):
        """
        Store checked parameters (use make() to build them).

        Args:
            path (str): the parameter file, or None.
            sections (dict): {section: {key: checked value}}.
            components (list of Component): with their final names.
            coronagraphs (list of SimpleNamespace): see _coronagraph().
            origins (dict): {(section, key): file path}.
            notes (iterable of str): notes from reading.
        """
        self.path = path
        self.sections = sections
        self.components = components
        self.coronagraphs = coronagraphs
        self._origins = origins
        self._read_notes = list(notes)  # from reading the file
        self.load_notes = list(notes)  # from reading and check()
        self.notes = list(notes)       # load_notes + the messages of the last resolve()
        self.provided = []      # (key, value, origin)
        self.resolved = {}      # key -> (value, origin)

    # ---- access
    def user_set(self, key):
        """
        Whether the user set a core key (even to "auto").

        Args:
            key (str): "section.name".

        Returns:
            bool
        """
        section, name = key.split(".")
        return name in self.sections.get(section, {})

    def get(self, key):
        """
        The value of a core key: resolved, else the user's, else the default.

        Args:
            key (str): "section.name".

        Returns:
            object: the value (None if no value and no default).
        """
        if key in self.resolved:
            return self.resolved[key][0]
        section, name = key.split(".")
        return self.sections.get(section, {}).get(name, _param(key).default)

    def origin(self, key):
        """
        Where the value of a core key comes from.

        Args:
            key (str): "section.name".

        Returns:
            str: the resolved origin (e.g. "header PIXSCALE of x.fits"), else "set in <file name>" (or
                "set in parameters" when built in Python), else "default".
        """
        if key in self.resolved:
            return self.resolved[key][1]
        section, name = key.split(".")
        if name in self.sections.get(section, {}):
            return "set in " + os.path.basename(self._origins.get((section, name), self.path or "parameters"))
        return "default"

    def enabled(self):
        """
        The components used.

        Returns:
            list of Component: those with enabled = True.
        """
        return [c for c in self.components if c.enabled]

    def families(self):
        """
        The unit family of each enabled component.

        Returns:
            dict: {component name: "physical", "relative" or None (decided by the input)}.
        """
        return {c.name: c.source.family_of() for c in self.enabled()}

    def star_flux_needed(self):
        """
        Whether the star's flux density is needed to put the components and the output in one unit:
        the components are in different families, or [display] unit asks for another family.

        Returns:
            bool
        """
        fam = set(f for f in self.families().values() if f)
        unit = self.get("display.unit")
        out_family = None if unit == AUTO else ("relative" if unit.startswith("contrast") else "physical")
        return len(fam) > 1 or bool(out_family and fam and out_family not in fam)

    def used_by(self, key):
        """
        Who uses a shared quantity.

        Args:
            key (str): "section.name".

        Returns:
            list of str: names of the enabled components requiring it, plus "the PSFs" (wavelength,
                diameter), "the image grid" (grid.*) or "unit conversion" (the star keys needed for its flux).
        """
        users = [c.name for c in self.enabled() if key in c.source.requires]
        if key in ("observation.wavelength", "observation.diameter"):
            users.append("the PSFs")
        elif key.startswith("grid."):
            users.append("the image grid")
        elif key.startswith("star.") and self.star_flux_needed():
            blackbody = not self.user_set("star.flux_jy")  # the flux from temperature, radius and distance
            # flux_jy is used when set; the other star keys when the flux is a blackbody's
            if (key == "star.flux_jy") != blackbody:
                users.append("unit conversion")
        return users

    # ---- what is missing, what is unused
    def check(self):
        """
        Check that everything a component needs is given; note the [star] keys nobody uses.

        A required shared key is given if it has a value (user or default) or an enabled component
        provides it. When the star's flux is needed, temperature or flux_jy must be set (with
        temperature, radius and distance are needed too). Notes are appended to load_notes; notes is
        reset to load_notes.

        Raises:
            ParamError: everything missing, one per line.
        """
        errors = []
        self.load_notes = list(self._read_notes)  # (so that checking again does not repeat the notes)
        for c in self.enabled():
            errors += ["[[source]] %s needs %s (%s): %s" % (c.name, p.name, p.unit or p.describe_type(), p.doc)
                       for p in c.source.missing()]
        needs = {}  # key -> names of who needs it
        for c in self.enabled():
            for key in c.source.requires:
                needs.setdefault(key, []).append(c.name)
        if self.star_flux_needed():
            if not (self.user_set("star.flux_jy") or self.user_set("star.temperature")):
                errors.append("the components are in different units (%s): the star's flux is needed to combine them; "
                              "set [star] temperature (with radius and distance) or flux_jy"
                              % ", ".join("%s: %s" % (n, "contrast" if f == "relative" else "physical")
                                          for n, f in self.families().items() if f))
            elif not self.user_set("star.flux_jy"):
                for key in ("star.radius", "star.distance"):
                    needs.setdefault(key, []).append("the blackbody flux (temperature)")
        provided = {key for c in self.enabled() for key in c.source.provides}
        for key, who in needs.items():
            if self.get(key) is None and key not in provided:
                section, name = key.split(".")
                errors.append("%s needs [%s] %s (%s)" % (" and ".join(who), section, name,
                                                         _param(key).unit or _param(key).describe_type()))
        if errors:
            raise ParamError("Problems in %s:\n  " % (self.path or "the parameters") + "\n  ".join(errors))
        undecided = None in self.families().values()  # a file decides its unit: the star's flux may be needed
        for key in ("star.distance", "star.radius", "star.temperature", "star.flux_jy"):
            if self.user_set(key) and not self.used_by(key) and not undecided:
                why = " (the output needs no unit conversion)" if key in ("star.temperature", "star.flux_jy") else ""
                self.load_notes.append("Note: [star] %s is not used by any enabled component%s: you can remove it"
                                       % (key.split(".")[1], why))
        self.notes = list(self.load_notes)

    # ---- the final values
    def provide(self, key, value, origin):
        """
        A source supplies a shared quantity (e.g. a header's pixel scale), used by resolve().

        Args:
            key (str): "section.name", e.g. "grid.pixel_scale".
            value (float or int): in the key's canonical unit.
            origin (str): where it comes from, e.g. "header PIXSCALE of x.fits".
        """
        self.provided.append((key, value, origin))

    def resolve(self):
        """
        The final shared values: the user's, else what a source provides, else the default.

        For [observation] and [star]: the user's value (not "auto"), else the first provided, else the
        default (wavelength: 1e-6 m assumed, with a note); the grid: see _resolve_grid(). Every
        disagreement is checked: over-specified settings (grid, temperature + flux_jy) and the provided
        values against the resolved ones (CONFLICT_RULES). Sets resolved and notes (load_notes + the
        messages).

        Returns:
            list of str: the messages (notes and warnings).

        Raises:
            ParamError: conflicts and missing grid settings, one per line ([checks] on_conflict = "warn"
                turns the CONFLICT_RULES errors into warnings).
        """
        self.resolved = {}  # resolving again starts from the user's settings, not the last result
        given = {}  # key -> [(value, origin)] in the order provided
        for key, value, origin in self.provided:
            given.setdefault(key, []).append((value, origin))
        errors, messages, res = [], [], {}
        for key in ("observation.diameter", "observation.wavelength", "star.distance", "star.radius",
                    "star.temperature", "star.flux_jy"):
            if self.user_set(key) and self.get(key) != AUTO:
                res[key] = (self.get(key), self.origin(key))
            elif given.get(key):
                res[key] = given[key][0]
            elif _param(key).default not in (None, AUTO):
                res[key] = (_param(key).default, "default")
        if "observation.wavelength" not in res:
            res["observation.wavelength"] = (1e-6, "assumed (nothing gives it)")
            messages.append("Note: no wavelength in the parameters or the image header: 1e-06 m assumed; "
                            "set [observation] wavelength")
        lambda_d = res["observation.wavelength"][0]/res["observation.diameter"][0]*RAD_TO_ARCSEC  # arcsec
        self._resolve_grid(res, given, lambda_d, errors, messages)
        self._check_star_flux(res, errors)
        policy = self.get("checks.on_conflict")
        for key, (tol, consequence, why) in CONFLICT_RULES.items():
            # every provided value against the resolved one (the same origin: that value itself)
            for other, other_origin in (given.get(key, []) if key in res else []):
                value, origin = res[key]
                if other_origin == origin or abs(other-value) <= tol*abs(value):
                    continue
                text = "%s: %.6g (%s) vs %.6g (%s): %s" % (key, value, origin, other, other_origin, why)
                if consequence == "warning":
                    messages.append("WARNING: " + text)
                elif policy == "error":
                    errors.append("CONFLICT %s; fix one of them, or set [checks] on_conflict = \"warn\"" % text)
                else:
                    messages.append("WARNING (accepted, [checks] on_conflict = \"warn\"): " + text)
        if errors:
            raise ParamError("Problems in %s:\n  " % (self.path or "the parameters") + "\n  ".join(errors))
        self.resolved = res
        self.notes = self.load_notes + messages
        return messages

    def _resolve_grid(self, res, given, lambda_d, errors, messages):
        """
        Resolve the grid: grid.size, grid.pixel_scale, grid.fov (and grid.sampling if set).

        The pixel scale, from the first that applies: pixel_scale or sampling (both set: they must
        agree); fov / size (the user's size, else the image's); the image's field / the user's size
        (size alone changes the resolution of the image's field); the image's pixel scale; else a
        field of 32 lambda/D across the size (with a note).
        The size: the user's (with fov and a user pixel scale too, the three must agree to 1.5 pixels,
        or GROUP_TOL); else fov / pixel scale, made odd; else the image's field of view at this pixel
        scale (keeping the image's parity); else the image's size.
        "The image" is the first source providing grid.size / grid.pixel_scale.

        Args:
            res (dict): {key: (value, origin)}; the grid keys are added to it.
            given (dict): {key: [(value, origin), ...]} provided by the sources.
            lambda_d (float): lambda/D in arcsec (to convert sampling and "N lambda/D").
            errors (list of str): problems found are appended to it (then res is left without the grid).
            messages (list of str): notes are appended to it.
        """
        user = {k: self.get("grid."+k) for k in ("size", "pixel_scale", "sampling", "fov")
                if self.user_set("grid."+k) and self.get("grid."+k) != AUTO}
        # (value, origin) of the image's size and pixel scale; (None, None) if no source provides them
        image_n, image_ps = (given.get("grid."+k, [(None, None)])[0] for k in ("size", "pixel_scale"))
        fov = user.get("fov")
        fov_origin = self.origin("grid.fov") + (" (%g lambda/D)" % fov if isinstance(fov, LambdaD) else "")
        fov = float(fov)*lambda_d if isinstance(fov, LambdaD) else fov  # arcsec, or None
        options = []  # the user's pixel scales: (arcsec, origin)
        if "pixel_scale" in user:
            options.append((user["pixel_scale"], self.origin("grid.pixel_scale")))
        if "sampling" in user:
            options.append((lambda_d/user["sampling"], "from sampling = %g" % user["sampling"]))
        if len(options) == 2:
            if abs(options[1][0]/options[0][0]-1) > GROUP_TOL:
                errors.append("[grid] pixel_scale = %g arcsec but sampling = %g gives %.6g arcsec: set only one"
                              % (options[0][0], user["sampling"], options[1][0]))
            else:
                messages.append("Note: [grid] pixel_scale and sampling agree; one of them is enough")
        # pixel and size are (value, origin) tuples from here on
        pixel = options[0] if options else None
        size = (user["size"], self.origin("grid.size")) if "size" in user else None
        if pixel is None and fov and (size or image_n[0]):
            size = size or image_n
            pixel = (fov/size[0], "from fov / size")
        if pixel is None and size and image_n[0] and image_ps[0]:  # the image's field at another resolution
            pixel = (image_n[0]*image_ps[0]/size[0], "the image's field / size")
        if pixel is None and image_ps[0]:
            pixel = image_ps
        if pixel is None:  # nothing gives the pixel scale: the old default field of view
            size = size or (image_n if image_n[0] else None)
            if size is None:
                errors.append("not enough to define the image: set [grid] size and pixel_scale (or sampling, or fov), "
                              "or add an image source")
                return
            fov, fov_origin = 32*lambda_d, "assumed: 32 lambda/D"
            pixel = (fov/size[0], "assumed: 32 lambda/D across")
            messages.append("Note: no pixel scale in the parameters or the image header: a field of 32 lambda/D "
                            "assumed; set [grid] fov or pixel_scale")
        if size is None:
            if fov:
                n = int(round(fov/pixel[0]))
                size = (n + (n % 2 == 0), "from fov / pixel_scale (odd)")
            elif image_n[0] and image_ps[0] and abs(image_ps[0]/pixel[0]-1) > 1e-9:
                n = int(round(image_n[0]*image_ps[0]/pixel[0]))
                size = (n + (n-image_n[0]) % 2, "the field of view of the image")  # keep the image's parity
            elif image_n[0]:
                size = image_n
            else:
                errors.append("not enough to define the image: set [grid] size or fov, or add an image source")
                return
        elif fov and len(options) and abs(size[0]-fov/pixel[0]) > max(1.5, GROUP_TOL*size[0]):
            # size, pixel scale and fov all set: they must agree in pixels, to 1.5 (fov / pixel scale rounded,
            # then made odd as above) or GROUP_TOL
            errors.append("[grid] size x pixel_scale = %.6g arcsec but fov = %.6g arcsec: set two of the three"
                          % (size[0]*pixel[0], fov))
        res["grid.size"], res["grid.pixel_scale"] = size, pixel
        res["grid.fov"] = (size[0]*pixel[0], fov_origin if fov else "from size x pixel_scale")
        if "sampling" in user:
            res["grid.sampling"] = (user["sampling"], self.origin("grid.sampling"))

    def _check_star_flux(self, res, errors):
        """
        If both temperature and flux_jy are set, check that they agree (to GROUP_TOL).

        Args:
            res (dict): {key: (value, origin)}: the resolved star and wavelength.
            errors (list of str): problems found are appended to it (also when radius or distance is
                missing to compute the blackbody).
        """
        if not (self.user_set("star.temperature") and self.user_set("star.flux_jy")):
            return
        if "star.radius" not in res or "star.distance" not in res:
            errors.append("[star] set temperature OR flux_jy, not both")
            return
        bb = blackbody_flux_jy(res["star.temperature"][0], res["star.radius"][0], res["star.distance"][0],
                               res["observation.wavelength"][0])
        if abs(bb/res["star.flux_jy"][0]-1) > GROUP_TOL:
            errors.append("[star] flux_jy = %g Jy but the blackbody of temperature = %g K gives %.4g Jy: set only one"
                          % (res["star.flux_jy"][0], res["star.temperature"][0], bb))

    # ---- provenance
    def table(self):
        """
        The provenance table of the resolved values (after resolve()).

        Returns:
            str: aligned columns "quantity, value, from, used by", one line per resolved key.
        """
        rows = [("quantity", "value", "from", "used by")]
        for key, (value, origin) in self.resolved.items():
            unit = _param(key).unit if _is_number(value) else None
            text = ("%.6g" % value if _is_number(value) else str(value)) + (" " + unit if unit else "")
            rows.append((key, text, origin, ", ".join(self.used_by(key)) or "-"))
        widths = [max(len(r[i]) for r in rows) for i in range(4)]
        return "\n".join("  ".join(c.ljust(w) for c, w in zip(r, widths)).rstrip() for r in rows)

    def used_toml(self, coronagraphs=None, method=None):
        """
        The values used, as a TOML file to save with the results.

        The shared sections hold the resolved values (with unit and origin as comments), the other
        sections what the user set, then every component (with all its settings) and every coronagraph
        (with all its parameters).

        Args:
            coronagraphs (list of SimpleNamespace): the coronagraphs a run used (type, name, values), or
                None for the parameters'.
            method (dict): [method] settings a run used instead of the parameters' ({key: value}), or None.

        Returns:
            str: the TOML text (ends with a newline).
        """
        stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        lines = ["# Values used by toycoronagraph (%s), from %s." % (stamp, self.path or "Python"),
                 "# Run toycoronagraph with this file to reproduce the run.", ""]
        for section in SECTIONS:
            if section in SHARED_SECTIONS:
                rows = [(key.split(".")[1], value, "%s; %s" % (_param(key).unit or "", origin))
                        for key, (value, origin) in self.resolved.items() if key.startswith(section+".")]
                rows += [(k, v, "") for k, v in self.sections.get(section, {}).items()  # e.g. [grid] resample
                         if "%s.%s" % (section, k) not in self.resolved and v != AUTO]
            else:
                values = dict(self.sections.get(section, {}), **((method or {}) if section == "method" else {}))
                rows = [(k, v, "") for k, v in values.items()]
            if rows:
                lines.append("[%s]" % section)
                lines += [_line(k, v, comment.lstrip("; "), width=34) for k, v, comment in rows]  # no "; " without unit
                lines.append("")
        builtin = _sources().BUILTIN_KINDS
        for c in self.components:
            lines.append("[[source]]")
            values = dict(kind=c.source.kind, name=c.name, enabled=c.enabled, scale=c.scale, method=c.method)
            if c.source.kind not in builtin:
                values["module"] = type(c.source).__module__
            values.update({k: v for k, v in c.source.settings.items() if v is not None})
            lines += [_line(k, v, "", width=34) for k, v in values.items()] + [""]
        for c in coronagraphs or self.coronagraphs:
            values = dict(type=c.type, name=c.name)
            values.update({p.name: c.values.get(p.name, p.default) for p in CORONAGRAPHS[c.type][1]})
            lines += ["[[coronagraph]]"] + [_line(k, v, "", width=34) for k, v in values.items()] + [""]
        return "\n".join(lines).rstrip()+"\n"

    def write_used(self, path):
        """
        Write used_toml() (with the parameters' coronagraphs and [method]) to a file.

        Args:
            path (str): the file (overwritten).
        """
        with open(path, "w") as f:
            f.write(self.used_toml())


def _param(key):
    """
    The declaration of a core key.

    Args:
        key (str): "section.name", e.g. "star.distance".

    Returns:
        Param

    Raises:
        KeyError: unknown section or name.
    """
    section, name = key.split(".")
    return {p.name: p for p in SECTIONS[section][1]}[name]


def _line(key, value, comment, width=24):
    """
    One TOML line "key = value", with an aligned comment.

    Args:
        key (str): the key.
        value (object): the value (written with _toml()).
        comment (str): the comment; "" = none.
        width (int): column where " # comment" starts (longer "key = value" texts push it right).

    Returns:
        str: the line.
    """
    text = "%s = %s" % (key, _toml(value))
    return ("%-*s # %s" % (width, text, comment)) if comment else text


def blackbody_flux_jy(temperature, radius_rsun, distance_pc, wavelength_m):
    """
    Flux density of a blackbody star: pi B_nu(T) (R/d)^2.

    Args:
        temperature (float): effective temperature in K.
        radius_rsun (float): stellar radius in R_sun.
        distance_pc (float): distance in pc.
        wavelength_m (float): wavelength in m.

    Returns:
        float: the flux density in Jy.
    """
    from astropy import units as u, constants as const
    from astropy.modeling.physical_models import BlackBody
    bnu = BlackBody(temperature*u.K)(wavelength_m*u.m)  # per steradian
    return float((bnu*u.sr*(radius_rsun*const.R_sun/(distance_pc*u.pc)).decompose()**2*3.141592653589793).to_value(u.Jy))


# ====================================================================== the old layout

def legacy_to_new(data, base_dir="."):
    """
    Map the old layout ([target], [telescope], [coronagraph], [method], [compute]) to the new one.

    [target] becomes one [[source]] (mcfost if input_format = "mcfost", mcfost_* keys are set or the
    file looks like an MCFOST cube; else fits_image), [telescope] wavelength / diameter go to
    [observation], target_pixel / visual_range (lambda/D) to [grid] size / fov, the [coronagraph]
    table to one [[coronagraph]] (type vortex if not set), [target] unit to [display] label;
    [method] and [compute] are copied. "auto" values are left out.

    Args:
        data (dict): the old-layout file content.
        base_dir (str): directory of the file (for relative FITS file names).

    Returns:
        dict: the content in the new layout (to check with from_data()).
    """
    target = dict(data.get("target", {}))
    telescope = data.get("telescope", {})
    out = {}
    fname = target.get("fits_filename", "")
    fmt = target.get("input_format", AUTO)
    mcfost_keys = [k for k in target if k.startswith("mcfost_")]
    if fmt == AUTO and fname:
        fmt = "mcfost" if (mcfost_keys or _looks_mcfost(fname, base_dir)) else "image"
    src = {"kind": "mcfost" if fmt == "mcfost" else "fits_image"}
    if src["kind"] == "fits_image":
        src["file"] = fname or "example"
        for old, new in [("fits_hdu", "hdu"), ("fits_slice", "slice"), ("input_unit", "input_unit")]:
            if old in target and target[old] != AUTO:
                src[new] = target[old]
    else:
        src["file"] = fname
        for old, new in [("mcfost_component", "component"), ("mcfost_inclination", "inclination"),
                         ("mcfost_azimuth", "azimuth")]:
            if old in target:
                src[new] = target[old]
        if "mcfost_component" not in target:
            src["component"] = "scattered_star"  # the old default
    if target.get("symmetric", AUTO) != AUTO:
        src["symmetric"] = target["symmetric"]
    out["source"] = [src]
    obs = {k: telescope[k] for k in ("wavelength", "diameter") if k in telescope}
    if obs:
        out["observation"] = obs
    grid = {}
    if target.get("target_pixel", AUTO) != AUTO:
        grid["size"] = target["target_pixel"]
    if telescope.get("visual_range", AUTO) != AUTO:
        grid["fov"] = "%r lambda/D" % float(telescope["visual_range"])
    if grid:
        out["grid"] = grid
    coro = dict(data["coronagraph"]) if isinstance(data.get("coronagraph"), dict) else {}
    if coro.get("type") is None:
        coro["type"] = "vortex"
    out["coronagraph"] = [coro]
    for section in ("method", "compute"):
        if data.get(section):
            out[section] = dict(data[section])
    if target.get("unit", "none") != "none":
        out["display"] = {"label": target["unit"]}
    return out


def _looks_mcfost(fname, base_dir):
    """
    Whether a FITS file looks like an MCFOST cube (5-D primary HDU with FLUX_1 in its header).

    Args:
        fname (str): the file name (".fits" added if it has no extension).
        base_dir (str): directory of a relative fname.

    Returns:
        bool: False also if the file cannot be read.
    """
    path = fname if os.path.isabs(fname) else os.path.join(base_dir, fname)
    if not os.path.splitext(path)[1]:
        path += ".fits"
    try:
        from astropy.io import fits
        with fits.open(path) as h:
            return h[0].header.get("NAXIS") == 5 and "FLUX_1" in h[0].header
    except Exception:
        return False
