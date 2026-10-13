"""
Parameters of Target (one image through one coronagraph), read from toycoronagraph_para.toml.

This is the parameter reader of the older, flat configuration used by Target (main.py);
the scene framework (toycoronagraph.observe) uses toycoronagraph.params directly.

Two layouts are read:
  * the scene layout of toycoronagraph.params ([observation], [grid], [star], [[source]],
    [[coronagraph]], ...): checked there, then mapped onto the flat configuration below
    (one image source; the first coronagraph);
  * the older layout ([target], [telescope], [coronagraph], [method], [compute]) in SCHEMA;
    its [coronagraph] type is "vortex" or "lazuli". A lazuli [coronagraph] has the keys of the
    scene layout's "lazuli" type and is checked with its declarations
    (toycoronagraph.params.CORONAGRAPHS); its own keys (LAZULI_KEYS) are refused under a vortex.
Every key is optional. Keys set to "auto" are filled from the input image and its FITS
header (resolve(), called by inputs.read_target). The values actually used are written to
<parameter file>_used.toml (write_used); the parameter file itself is not modified unless
[run] write_back = true. An old-style toycoronagraph_para.py is still read (with a notice)
when no .toml file exists. If there is no parameter file at all, a minimal one is created
(toycoronagraph.params_docs.new_config).
"""
import os
import importlib.util
from types import SimpleNamespace
try:
    import tomllib # Python >= 3.11
except ModuleNotFoundError: # pragma: no cover
    import tomli as tomllib
from toycoronagraph.params import CORONAGRAPHS, Param, convert, coronagraph_problems

PARA_FILE = "toycoronagraph_para.toml"
LEGACY_PARA_FILE = "toycoronagraph_para.py"
RAD_TO_ARCSEC = 206264.806247
LIGHT_SPEED = 299792458 # m/s
JY = 1e26 # W / m^2 / Hz to Jy
MISMATCH_TOL = 0.01 # relative difference between parameters and header that is reported
RESOLVED_MARK = "# ---- Values used in the last run, after checking the input (written by toycoronagraph; not read back) ----"

AUTO = "auto"
UNITS = ("W/m^2/pixel", "Jy/arcsec^2", "Jy/pixel")
# MCFOST image cubes (RT.fits.gz): axis 0 holds these 8 planes (FLUX_1 ... FLUX_8)
MCFOST_COMPONENTS = {"total": 0, "I": 0, "Q": 1, "U": 2, "V": 3, "direct_star": 4, "scattered_star": 5,
                     "direct_thermal": 6, "scattered_thermal": 7}
# "disk": all the light of the disk, i.e. the total minus the direct starlight (the star itself)
MCFOST_CHOICES = ("disk",)+tuple(MCFOST_COMPONENTS)
# [coronagraph] keys used only by type = "lazuli": declared (and checked) as in toycoronagraph.params
LAZULI_KEYS = ("wavelengths", "throughput", "opd", "pupil_samples")

# section -> key -> (default, allowed values / type / "auto|<type>" / "index list" / a params.Param, description)
SCHEMA = {
    "target": {
        "fits_filename": ("", str, "FITS file of the target (\".fits\" is added if no extension is given); \"\" uses the example target"),
        "input_format": (AUTO, (AUTO, "image", "mcfost"), "\"image\": a 2-D image; \"mcfost\": an MCFOST image cube; \"auto\" recognizes MCFOST cubes"),
        "mcfost_component": ("scattered_star", MCFOST_CHOICES, "MCFOST plane: \"disk\" = total minus the star; \"total\" (I) includes the star itself"),
        "mcfost_inclination": (0, int, "MCFOST: index of the inclination"),
        "mcfost_azimuth": (0, int, "MCFOST: index of the azimuth"),
        "fits_hdu": (0, int, "FITS extension (HDU) holding the data"),
        "fits_slice": (AUTO, "index list", "explicit indices into the leading axes, e.g. [5, 0, 0] (overrides the MCFOST keys); \"auto\": from input_format"),
        "input_unit": (AUTO, (AUTO,)+UNITS, "unit of the data; \"auto\" reads the BUNIT header keyword"),
        "target_pixel": (AUTO, "auto|int", "pixels on one edge of the (square) image; \"auto\" uses the image size"),
        "unit": ("none", str, "colorbar label; \"none\" uses Jy/arcsec^2"),
        "symmetric": (AUTO, (AUTO, True, False), "is the target circularly symmetric about the star? \"auto\" measures it"),
    },
    "telescope": {
        "wavelength": (AUTO, "auto|float", "observing wavelength in m; \"auto\" reads WAVELEN (m) or WAVE (um) from the header, else 1e-6"),
        "diameter": (2.4, float, "aperture diameter D, in meters"),
        "visual_range": (AUTO, "auto|float", "width of the image in lambda/D; \"auto\" follows the header's pixel scale (PIXSCALE, or CDELT), else 32"),
    },
    "coronagraph": {
        "type": ("vortex", ("vortex", "lazuli"), "coronagraph model: \"vortex\", or \"lazuli\" (the Lazuli telescope, prysm)"),
        "charge": (2, int, "vortex charge (positive, even); lazuli: 0 = the telescope alone (its default when "
                           "charge is left out) or even; used when a method is called without one"),
        "lyot_mask_size": (0.95, float, "Lyot stop diameter relative to the aperture"),
        "symmetric": (AUTO, (AUTO, True, False), "are the PSFs circularly symmetric? \"auto\" uses the optics model"),
        **{p.name: (p.default, p, "lazuli: " + p.doc) for p in CORONAGRAPHS["lazuli"][1] if p.name in LAZULI_KEYS},
    },
    "method": {
        "method": (AUTO, (AUTO, "circular", "spatial"), "\"auto\": circular if the target and the PSFs are both symmetric, else spatial"),
        "spatial_mode": (AUTO, (AUTO, "rotate", "direct"), "\"auto\": rotate for symmetric PSFs, else direct"),
        "error_target": (0.02, float, "target relative error of the final image (spatial: sets kernel size and grid)"),
        "symmetry_tol": (0.01, float, "measured target asymmetry above which \"auto\" treats it as not symmetric"),
        "rot_number": (0, int, "circular method: number of rotations; 0 uses target_pixel"),
    },
    "compute": {
        "num_cores": (-1, int, "CPU cores for parallel work; -1 uses all"),
        "use_gpu": (False, bool, "use CuPy for the spatial convolutions"),
        "max_kernel_gb": (4.0, float, "largest spatial kernel cube to build, in GB"),
        "psf_cache_dir": ("", str, "folder for computed PSFs; \"\" uses the current folder"),
    },
}

USAGE_NOTE = """\
# Using your own target
#   * An MCFOST image (e.g. data_0.6/RT.fits.gz): set fits_filename to it; the cube is
#     recognized, and mcfost_component / mcfost_inclination / mcfost_azimuth pick the plane.
#     Use "disk" (or "scattered_star", ...) rather than "total": the total includes the star
#     itself, which the coronagraph is meant to remove. The header gives the unit, pixel
#     scale and wavelength.
#   * Any other image: save a FITS file with a square 2-D array (or a cube and set
#     fits_slice), the star at the centre of the array ((N-1)/2 in 0-based pixels), and
#     header keywords BUNIT (W/m^2/pixel, Jy/arcsec^2 or Jy/pixel), PIXSCALE (arcsec per
#     pixel; or CDELT1/CDELT2 in degrees) and WAVELEN (m) - or set input_unit / wavelength
#     / visual_range here instead.
#   The image pixel scale must match the PSF sampling: with visual_range = "auto" it does.
#   The values actually used are saved to <this file's name>_used.toml (this file is not modified).
"""


def _toml_value(value):
    """
    A value as TOML text (lists/tuples as arrays, booleans as true/false, strings quoted).

    Args:
        value: A str, bool, int, float, or a list/tuple of these.

    Returns:
        str: The TOML literal (strings are not escaped).
    """
    if isinstance(value, (list, tuple)):
        return "[%s]" % ", ".join(_toml_value(v) for v in value)
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return '"%s"' % value
    return repr(value)


def example_text():
    """
    The text of a commented example toycoronagraph_para.toml in the older layout (all defaults).

    Returns:
        str: The file content: a header with USAGE_NOTE, then every SCHEMA section and key
            with its default and description (the LAZULI_KEYS commented out: the type is "vortex").
    """
    lines = ["# Toy Coronagraph parameters. Every key is optional; these are the defaults.", "#", USAGE_NOTE]
    for section, keys in SCHEMA.items():
        lines.append("[%s]" % section)
        for key, (default, _, description) in keys.items():
            line = "%-18s = %-16s # %s" % (key, _toml_value(default), description)
            lines.append("# " + line if section == "coronagraph" and key in LAZULI_KEYS else line)
        lines.append("")
    return "\n".join(lines)


def write_example(path=PARA_FILE):
    """
    Create an example parameter file in the older layout (example_text); overwrites path.

    Args:
        path (str): The file to write.
    """
    with open(path, "w") as f:
        f.write(example_text())


def _is_number(value, kind):
    """
    Whether value is a number of the given kind (booleans are not numbers).

    Args:
        value: Any value.
        kind (type): int (only int accepted) or float (int or float accepted).

    Returns:
        bool
    """
    if isinstance(value, bool):
        return False
    return isinstance(value, int) if kind is int else isinstance(value, (int, float))


def _check(section, key, value):
    """
    Check one value against SCHEMA and convert it.

    Args:
        section (str): Section name in SCHEMA.
        key (str): Key name in that section.
        value: The value read from the file (or the default).

    Returns:
        The value: "auto" kept as is; numbers converted to int/float as the schema says; keys
            declared by a params.Param converted by params.convert (e.g. "632 nm" to m).

    Raises:
        ValueError: If the value does not match the schema.
    """
    _, allowed, _ = SCHEMA[section][key]
    if isinstance(allowed, Param):
        return convert(allowed, value, "[%s]" % section)
    if allowed == "index list":
        if value == AUTO or (isinstance(value, list) and all(_is_number(v, int) for v in value)):
            return value
        raise ValueError("[%s] %s = %r: must be \"auto\" or a list of integers, e.g. [5, 0, 0]" % (section, key, value))
    if isinstance(allowed, str) and allowed.startswith("auto|"):
        kind = int if allowed == "auto|int" else float
        if value == AUTO:
            return value
        if _is_number(value, kind):
            return kind(value)
        raise ValueError("[%s] %s = %r: must be \"auto\" or a%s %s" % (section, key, value, "n" if kind is int else "", kind.__name__))
    if isinstance(allowed, tuple):
        # booleans must match booleans exactly (in Python, True == 1 and False == 0)
        if value not in allowed or (isinstance(value, bool) != any(isinstance(a, bool) and a == value for a in allowed)):
            raise ValueError("[%s] %s = %r: must be one of %s" % (section, key, value, ", ".join(map(_toml_value, allowed))))
        return value
    if allowed is float and _is_number(value, float):
        return float(value)
    # (bool is a subclass of int: a boolean is only accepted where bool is expected)
    if not isinstance(value, allowed) or (allowed is not bool and isinstance(value, bool)):
        raise ValueError("[%s] %s = %r: must be a %s" % (section, key, value, allowed.__name__))
    return value


def from_dict(data, source="parameters"):
    """
    Validate a parsed TOML dict and return the flat configuration (not yet resolved).

    Args:
        data (dict): {section: {key: value}}; a [resolved] section is ignored.
        source (str): Name used in messages (the file path).

    Returns:
        SimpleNamespace: every parameter, flat (the two `symmetric` keys become
            image_symmetric and psf_symmetric, [coronagraph] type becomes
            coronagraph_type, fits_filename "" becomes None, diameter is also
            stored as D). "auto" values stay "auto" until resolve(). Also has
            source (str) and legacy (False).

    Raises:
        ValueError: For an unknown section or key, a lazuli-only key (LAZULI_KEYS) under
            another type, a value that does not match SCHEMA (lazuli: the declarations
            of that type), a vortex charge that is not positive and even, a lazuli
            charge that is not even, or a lazuli throughput without one weight per
            wavelength.
    """
    for section, keys in data.items():
        if section == "resolved":
            continue
        if section not in SCHEMA:
            raise ValueError("%s: unknown section [%s]; known: %s" % (source, section, ", ".join(SCHEMA)))
        if not isinstance(keys, dict):
            raise ValueError("%s: [%s] must be a table of keys" % (source, section))
        for key in keys:
            if key not in SCHEMA[section]:
                raise ValueError("%s: unknown key [%s] %s; known: %s" % (source, section, key, ", ".join(SCHEMA[section])))
    coronagraph = data.get("coronagraph", {})
    # a lazuli [coronagraph] is checked with that type's declarations (its own defaults: charge 0, the telescope alone)
    lazuli = {p.name: p for p in CORONAGRAPHS["lazuli"][1]} if coronagraph.get("type") == "lazuli" else {}
    extra = [key for key in LAZULI_KEYS if key in coronagraph and not lazuli]
    if extra:
        raise ValueError("%s: [coronagraph] %s: only used with type = \"lazuli\"" % (source, ", ".join(extra)))
    rename = {("target", "symmetric"): "image_symmetric", ("coronagraph", "symmetric"): "psf_symmetric",
              ("coronagraph", "type"): "coronagraph_type"}
    cfg = SimpleNamespace(source=source, legacy=False)
    for section, keys in SCHEMA.items():
        for key, (default, _, _) in keys.items():
            declared = lazuli.get(key) if section == "coronagraph" else None
            if declared is None:
                value = _check(section, key, data.get(section, {}).get(key, default))
            else:
                value = convert(declared, coronagraph.get(key, declared.default), "[coronagraph]")
            setattr(cfg, rename.get((section, key), key), value)
    cfg.D = cfg.diameter
    cfg.fits_filename = cfg.fits_filename or None
    if cfg.coronagraph_type == "vortex" and (cfg.charge <= 0 or cfg.charge%2):
        raise ValueError("[coronagraph] charge = %r: a vortex charge must be positive and even" % cfg.charge)
    if lazuli:
        problems = coronagraph_problems("lazuli", {key: getattr(cfg, key) for key in ("charge", "wavelengths", "throughput")})
        if problems:
            raise ValueError("[coronagraph] " + "; ".join(problems))
    return cfg


def _header_unit(header):
    """
    The unit named by the header's BUNIT, restricted to UNITS.

    Args:
        header: FITS header (dict-like).

    Returns:
        tuple:
            unit (str or None): One of UNITS, or None for any other (or no) BUNIT.
            raw (str): The BUNIT text as written ("" if absent).
    """
    from toycoronagraph.units import header_unit
    unit, raw = header_unit(header)
    return (unit if unit in UNITS else None), raw


def _header_pixel_scale(header):
    """
    The pixel scale from PIXSCALE (arcsec) or CDELT2/CDELT1 (unit from CUNIT, default deg).

    Args:
        header: FITS header (dict-like).

    Returns:
        tuple:
            scale (float or None): arcsec per pixel (absolute value); None if not found
                or the CUNIT is not understood.
            key (str or None): Where it came from, e.g. "PIXSCALE" or "CDELT2 (deg)".
    """
    if "PIXSCALE" in header:
        return abs(float(header["PIXSCALE"])), "PIXSCALE"
    for axis in ("2", "1"):
        if "CDELT" + axis in header:
            unit = str(header.get("CUNIT" + axis, "deg")).strip().lower()
            factor = {"deg": 3600.0, "arcmin": 60.0, "arcsec": 1.0, "mas": 1e-3, "rad": RAD_TO_ARCSEC}.get(unit)
            if factor is not None:
                return abs(float(header["CDELT" + axis]))*factor, "CDELT%s (%s)" % (axis, unit)
    return None, None


def _header_wavelength(header):
    """
    The wavelength from WAVELEN (m) or MCFOST's WAVE (um).

    Args:
        header: FITS header (dict-like).

    Returns:
        tuple:
            wavelength (float or None): In m (WAVE rounded to 12 significant digits after
                the um -> m conversion); None if neither keyword is present.
            key (str or None): "WAVELEN" or "WAVE (um)".
    """
    if "WAVELEN" in header:
        return float(header["WAVELEN"]), "WAVELEN"
    if "WAVE" in header:
        return float("%.12g" % (float(header["WAVE"])*1e-6)), "WAVE (um)"
    return None, None


def _differs(a, b):
    """
    Whether a differs from the reference b by more than MISMATCH_TOL (relative to b).

    Args:
        a, b (float): The value and the reference.

    Returns:
        bool
    """
    return abs(a-b) > MISMATCH_TOL*abs(b)


def resolve(cfg, image_shape, header=None):
    """
    Fill the "auto" parameters from the input image and its header, check the explicit
    ones against the header, and compute the derived quantities (psf_scale, F_transfer,
    psf_range, rot_number).

    Sets (in cfg): target_pixel, wavelength (m), visual_range (lambda/D across the image),
    psf_scale (arcsec/pixel), F_transfer (factor from nu F_nu in W/m^2/pixel to F_nu in
    Jy/arcsec^2), psf_range (lambda/D, half of visual_range), input_unit, rot_number
    (target_pixel if <= 0), sources ({parameter: where its value came from}) and notes.
    For a legacy (.py) configuration, wavelength, visual_range, psf_scale, F_transfer and
    psf_range are kept as read.

    Args:
        cfg (SimpleNamespace): From load_config / from_dict; modified in place.
        image_shape (tuple): Shape of the (sliced) input image; must be square and 2-D.
        header: The image's FITS header (a dict-like), or None.

    Returns:
        list: Messages: a summary of the input (where each value came from), then any
            mismatch (WARNING lines). Also stored as cfg.notes, with cfg.sources.

    Raises:
        ValueError: If the image is not square and 2-D, target_pixel does not match it,
            or input_unit is "auto" and the header has no known BUNIT.
    """
    header = header or {}
    notes = []
    sources = {}
    if len(image_shape) != 2 or image_shape[0] != image_shape[1]:
        raise ValueError("the input image must be square and 2-D (got shape %s); check fits_hdu, fits_slice "
                         "or the MCFOST keys" % (image_shape,))
    n = image_shape[0]
    if cfg.target_pixel == AUTO:
        cfg.target_pixel, sources["target_pixel"] = n, "image size"
    elif cfg.target_pixel != n:
        raise ValueError("target_pixel = %d but the input image is %d x %d; set target_pixel = \"auto\" or match the image"
                         % (cfg.target_pixel, n, n))
    else:
        sources["target_pixel"] = "parameters (matches the image)"

    if cfg.legacy: # an old .py file: its own derived values are used as they are
        sources.update(wavelength="old .py file", visual_range="old .py file", psf_scale="old .py file")
    else:
        h_wave, h_wave_key = _header_wavelength(header)
        if cfg.wavelength == AUTO:
            cfg.wavelength, sources["wavelength"] = (h_wave, "header " + h_wave_key) if h_wave else (1.0e-6, "default (no wavelength in the header)")
        else:
            sources["wavelength"] = "parameters"
            if h_wave and _differs(cfg.wavelength, h_wave):
                notes.append("WARNING: wavelength %.4g m (parameters) differs from the header's %.4g m (%s)" % (cfg.wavelength, h_wave, h_wave_key))
        lambda_d = cfg.wavelength/cfg.diameter*RAD_TO_ARCSEC # arcsec
        h_scale, h_scale_key = _header_pixel_scale(header)
        if cfg.visual_range == AUTO:
            if h_scale:
                cfg.psf_scale = h_scale
                cfg.visual_range = n*h_scale/lambda_d
                sources["visual_range"] = sources["psf_scale"] = "header " + h_scale_key
            else:
                cfg.visual_range = 32.0
                cfg.psf_scale = lambda_d*cfg.visual_range/n
                sources["visual_range"] = sources["psf_scale"] = "default (no pixel scale in the header)"
        else:
            cfg.psf_scale = lambda_d*cfg.visual_range/n
            sources["visual_range"] = "parameters"
            sources["psf_scale"] = "from visual_range"
            if h_scale and _differs(cfg.psf_scale, h_scale):
                notes.append("WARNING: the PSFs are sampled at %.5g arcsec/pixel (visual_range = %g) but the image header says %.5g (%s): "
                             "the target would be scaled wrongly; use visual_range = \"auto\"" % (cfg.psf_scale, cfg.visual_range, h_scale, h_scale_key))
        cfg.F_transfer = JY*cfg.wavelength/LIGHT_SPEED/cfg.psf_scale**2 # nu F_nu (W/m^2/pixel) to F_nu (Jy/arcsec^2)
        cfg.psf_range = cfg.visual_range/2

    h_unit, h_unit_raw = _header_unit(header)
    if cfg.input_unit == AUTO:
        if h_unit is None:
            raise ValueError("the unit of the input image is unknown (BUNIT = %r); set [target] input_unit to one of %s, "
                             "or the BUNIT header keyword" % (h_unit_raw, ", ".join(UNITS)))
        cfg.input_unit, sources["input_unit"] = h_unit, "header BUNIT"
    else:
        sources["input_unit"] = "parameters"
        if h_unit is not None and h_unit != cfg.input_unit:
            notes.append("WARNING: input_unit %s (parameters) differs from the header's BUNIT = %r" % (cfg.input_unit, h_unit_raw))
    if cfg.rot_number <= 0:
        cfg.rot_number = cfg.target_pixel
    cfg.sources = sources
    notes.insert(0, "Input: %s, %d x %d, %s (%s); %.5g arcsec/pixel (%s) = %.4g lambda/D across at %.4g m (%s)"
                 % (getattr(cfg, "input_description", "image"), n, n, cfg.input_unit, sources["input_unit"],
                    cfg.psf_scale, sources.get("psf_scale"), cfg.visual_range, cfg.wavelength, sources.get("wavelength")))
    cfg.notes = notes
    return notes


def used_path(cfg):
    """
    Where the values used are saved: <parameter file name>_used.toml next to it.

    Args:
        cfg (SimpleNamespace): A configuration (uses cfg.source).

    Returns:
        str: The path; "toycoronagraph_para_used.toml" when cfg.source is not a .toml or .py file.
    """
    source = str(cfg.source)
    stem = os.path.splitext(source)[0] if source.endswith((".toml", ".py")) else "toycoronagraph_para"
    return stem + "_used.toml"


def write_used(cfg):
    """
    Save the values used (after resolve()) to used_path(cfg); with write_back, also into the parameter file.

    Args:
        cfg (SimpleNamespace): A resolved configuration; cfg.write_back (if present and true)
            also writes the [resolved] section into cfg.source.

    Returns:
        str: The path written (used_path(cfg)).
    """
    path = used_path(cfg)
    with open(path, "w") as f:
        f.write("# Values used by toycoronagraph for %s (written after each run; not read back)\n" % cfg.source)
    write_resolved(cfg, path)
    if getattr(cfg, "write_back", False):
        write_resolved(cfg)
    return path


def write_resolved(cfg, path=None):
    """
    Write the values used (after resolve()) into a [resolved] section at the end of a
    .toml file, replacing the previous one; the rest of the file is kept as is.

    Args:
        cfg (SimpleNamespace): A resolved configuration.
        path (str): The .toml file; default cfg.source (skipped if it is not a .toml file).
    """
    path = path or cfg.source
    if not str(path).endswith(".toml") or not os.path.exists(path):
        return
    with open(path) as f:
        text = f.read()
    if RESOLVED_MARK in text:
        text = text[:text.index(RESOLVED_MARK)]
    text = text.rstrip("\n")+"\n\n"
    lines = [RESOLVED_MARK, "[resolved]"]
    used = getattr(cfg, "used_slice", [])
    values = [("input", getattr(cfg, "input_description", "")), ("fits_slice", used if isinstance(used, list) else str(used)),
              ("input_unit", cfg.input_unit), ("target_pixel", cfg.target_pixel), ("wavelength", cfg.wavelength),
              ("visual_range", round(cfg.visual_range, 6)), ("psf_scale", cfg.psf_scale), ("psf_range", round(cfg.psf_range, 6)),
              ("rot_number", cfg.rot_number)]
    for key, value in values:
        source = cfg.sources.get(key, "")
        lines.append("%-13s = %-28s%s" % (key, _toml_value(value), ("# " + source) if source else ""))
    for note in cfg.notes[1:]:
        lines.append("# " + note)
    with open(path, "w") as f:
        f.write(text+"\n".join(lines)+"\n")


def _from_legacy(path):
    """
    Read an old-style toycoronagraph_para.py (its derived values are taken as given).

    The file is executed as a Python module; its variables are mapped onto the older
    layout and checked with from_dict. psf_scale, F_transfer and psf_range are taken
    from the file if it defines them, else derived as in resolve().

    Args:
        path (str): The .py parameter file.

    Returns:
        SimpleNamespace: The configuration, with legacy = True.
    """
    spec = importlib.util.spec_from_file_location("toycoronagraph_para_legacy", path)
    par = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(par)
    data = {"target": {"target_pixel": int(par.target_pixel)},
            "telescope": {"wavelength": float(getattr(par, "wavelength", 1.0e-6)),
                          "diameter": float(getattr(par, "D", 2.4)),
                          "visual_range": float(getattr(par, "visual_range", 32.0))},
            "coronagraph": {"type": getattr(par, "coronagraph_type", "vortex"),
                            "lyot_mask_size": float(getattr(par, "lyot_mask_size", 0.95))},
            "method": {"rot_number": int(getattr(par, "rot_number", 0))},
            "compute": {"num_cores": int(getattr(par, "num_cores", -1))}}
    if getattr(par, "fits_filename", None):
        data["target"]["fits_filename"] = par.fits_filename
    if getattr(par, "unit", "none") != "none":
        data["target"]["unit"] = par.unit
    cfg = from_dict(data, source=path)
    cfg.legacy = True
    lambda_d = cfg.wavelength/cfg.diameter*RAD_TO_ARCSEC
    cfg.psf_scale = getattr(par, "psf_scale", lambda_d*cfg.visual_range/cfg.target_pixel)
    cfg.F_transfer = getattr(par, "F_transfer", JY*cfg.wavelength/LIGHT_SPEED/cfg.psf_scale**2)
    cfg.psf_range = getattr(par, "psf_range", cfg.visual_range/2)
    return cfg


def from_params(params):
    """
    The old flat configuration (toycoronagraph.config) for Target, from new-layout parameters.

    Only files with exactly one enabled component, a fits_image or mcfost source, can be
    mapped; Target uses the first [[coronagraph]] (a note is added to params.notes when
    there are several).

    Args:
        params (toycoronagraph.params.Params): Checked new-layout parameters (params.load).

    Returns:
        SimpleNamespace: The configuration (see from_dict), with params (the Params object)
            and write_back ([run] write_back) added.

    Raises:
        ParamError: When the file needs the scene framework (not exactly one image source).
    """
    from toycoronagraph.params import ParamError
    comps = params.enabled()
    if len(comps) != 1 or comps[0].source.kind not in ("fits_image", "mcfost"):
        raise ParamError("%s: Target simulates one image source (fits_image or mcfost); this file has %s. "
                         "Use toycoronagraph.observe() for scenes with several components."
                         % (params.path, ", ".join(c.name for c in comps) or "none"))
    kind, s = comps[0].source.kind, comps[0].source.settings
    # The image itself is read through the scene (toycoronagraph.scene); these keys only describe it
    target = {"fits_filename": "" if s["file"] == "example" else s["file"],
              "input_format": "image" if kind == "fits_image" else "mcfost"}
    if kind == "fits_image":
        unit = s["input_unit"] if s["input_unit"] in UNITS else AUTO
        target.update(fits_hdu=s["hdu"], fits_slice=s["slice"], input_unit=unit)
    else:
        target.update(mcfost_component=s["component"], mcfost_inclination=s["inclination"], mcfost_azimuth=s["azimuth"])
    target["symmetric"] = s["symmetric"]
    size = params.get("grid.size")
    target["target_pixel"] = size if isinstance(size, int) else AUTO
    if params.get("display.label"):
        target["unit"] = params.get("display.label")
    telescope = {"wavelength": params.get("observation.wavelength"), "diameter": params.get("observation.diameter")}
    c = params.coronagraphs[0]
    coronagraph = {"type": c.type, **{k: v for k, v in c.values.items() if k not in ("type", "name")}}
    method = {k: v for k, v in params.sections.get("method", {}).items() if k not in ("combine", "error_metric")}
    data = {"target": target, "telescope": telescope, "coronagraph": coronagraph, "method": method,
            "compute": dict(params.sections.get("compute", {}))}
    cfg = from_dict(data, source=params.path)
    cfg.params = params
    cfg.write_back = params.get("run.write_back")
    if len(params.coronagraphs) > 1:
        params.notes.append("Note: Target uses the first [[coronagraph]] (%s) as its default; observe() uses all of them"
                            % c.name)
    return cfg


def load_config(path=None, verbose=True):
    """
    Load the parameters (not yet resolved against the input image; see resolve()).

    With no path: toycoronagraph_para.toml in the current directory; if it does not
    exist, an old toycoronagraph_para.py (with a notice); if neither exists, a minimal
    new-layout .toml is created (params_docs.new_config) and used.

    Args:
        path (str): A .toml (or legacy .py) parameter file.
        verbose (bool): Print notices.

    Returns:
        SimpleNamespace: The configuration (see from_dict; from_params for a new-layout file).
    """
    if path is None:
        if os.path.exists(PARA_FILE):
            path = PARA_FILE
        elif os.path.exists(LEGACY_PARA_FILE):
            if verbose:
                print("Reading the old-style %s; toycoronagraph.config.write_example() creates a %s to replace it"
                      % (LEGACY_PARA_FILE, PARA_FILE))
            return _from_legacy(LEGACY_PARA_FILE)
        else:
            from toycoronagraph.params_docs import new_config
            new_config(["fits_image"], ["vortex:2"], path=PARA_FILE)
            if verbose:
                print("There is no %s file, created a minimal one (python -m toycoronagraph new --help: other templates)"
                      % PARA_FILE)
            path = PARA_FILE
    if path.endswith(".py"):
        return _from_legacy(path)
    with open(path, "rb") as f:
        data = tomllib.load(f)
    from toycoronagraph import params
    if params.is_legacy(data) and "include" not in data:
        return from_dict(data, source=path)
    p = params.load(path)
    cfg = from_params(p)
    if verbose:
        for note in p.notes:
            print(note)
    return cfg
