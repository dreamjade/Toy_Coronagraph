"""
Text generated from the parameter declarations (toycoronagraph.params and the source
modules), so it never goes out of date:

    new_config()          a parameter file holding only the chosen modules (+ the JSON schema)
    json_schema()         the JSON schema: editor hover, completion and checks (Even Better TOML)
    write_schema(path)    the JSON schema written to a file
    reference_markdown()  the parameter reference tables
    describe(query)       where a key is set, its unit, default, allowed values and users
    migrate(path)         an old-layout file rewritten in the new layout

Command line: python -m toycoronagraph new | describe | migrate | schema | reference
"""
import os
import json
import difflib
import datetime
try:
    import tomllib # Python >= 3.11
except ModuleNotFoundError: # pragma: no cover
    import tomli as tomllib
from toycoronagraph.params import (AUTO, SCHEMA_FILE, SECTIONS, SOURCE_COMMON, CORONAGRAPHS, CORONAGRAPH_COMMON,
                                   CONFLICT_RULES, LambdaD, ParamError, _QUANTITY, _toml, _line, _sources,
                                   _source_class, from_data, is_legacy, legacy_to_new)

# The [observation] wavelength of templates with a lazuli coronagraph: Lazuli's design wavelength (config_um
# telescope.general.design_wavelength), in the WCC's r' band, where its "auto" wavelengths and throughput have light
LAZULI_WAVELENGTH = 633e-9


def _family(cls, settings):
    """
    The unit family of a module with these settings.

    Args:
        cls (type): a Source subclass.
        settings (dict): {key: value}; keys that are not cls's parameters and None values are ignored.

    Returns:
        str or None: "physical", "relative", or None (decided by the input).
    """
    names = {p.name for p in cls.params}
    return cls(**{k: v for k, v in settings.items() if k in names and v is not None}).family_of()


def _comment(p, used_by=None):
    """
    The comment of a key in generated files: "unit; doc. Used by: ...".

    Args:
        p (Param): the declaration.
        used_by (list of str): who uses it (e.g. source kinds); None or [] = not shown.

    Returns:
        str
    """
    unit = p.unit or ""
    text = (unit + "; " if unit else "") + p.doc
    if used_by:
        text += ". Used by: " + ", ".join(used_by)
    return text


def _example_value(p):
    """
    The value written for a key in generated files: its default, else its example.

    Args:
        p (Param): the declaration.

    Returns:
        object
    """
    return p.default if p.default is not None else p.example


def _odd_image_size(srcs):
    """
    The [grid] size of a template whose grid comes from an image source, for a lazuli coronagraph (which
    needs an odd size): the first image's own size, made odd.

    Args:
        srcs (list): (Source class, {key: value} given) of the template's sources, one at least providing
            grid.size (an image source). The first of those is read, its keys not given taking the values
            the template writes for them.

    Returns:
        int or str: the image's size minus 1 when it is even (the image is resampled onto it); "auto" when
            it is odd, or when the image cannot be read yet (e.g. a file still to be made: the run then names
            the odd sizes to use).
    """
    cls, given = next((cls, s) for cls, s in srcs if "grid.size" in cls.provides)
    try:
        provided = cls(**{p.name: given.get(p.name, _example_value(p)) for p in cls.params}).inspect()
    except (OSError, ValueError):
        return AUTO
    size = next(value for key, value, _ in provided if key == "grid.size")
    return size-1 if size % 2 == 0 else AUTO


def new_config(sources=("fits_image",), coronagraphs=("vortex:2",), path=None, full=False, schema=True, header=True):
    """
    A parameter file holding only what the chosen modules use.

    [observation] and [grid] are always written ("auto" from an image source if one provides the pixel
    scale, else sampling and fov); [star] holds the keys the sources require (and temperature, radius,
    distance if the sources are in different unit families). Then one [[source]] table per source
    (its keys without a default are marked "(required)"), one [[coronagraph]] per coronagraph, and
    [method] if there are several sources. With a lazuli coronagraph, the wavelength is
    LAZULI_WAVELENGTH (in the r' band of the Lazuli camera) and the image's size is made odd
    (_odd_image_size), so the file runs as written.

    Args:
        sources (list): Source kinds (str), or dicts {kind: ..., key: value} whose values are written
            instead of the defaults.
        coronagraphs (list): "type:charge" strings (e.g. "vortex:2"; "" type = vortex), or dicts
            {type: ..., charge: ...}.
        path (str): Write it there (and the JSON schema next to it, if schema); None = only return it.
        full (bool): Also list the optional (advanced) keys and the other sections, commented out.
        schema (bool): Refer to the JSON schema in the header ("#:schema" line) and write it next to path.
        header (bool): Start with the header comment lines.

    Returns:
        str: the TOML text.

    Raises:
        ParamError: an unknown source kind.
    """
    srcs = []  # (Source class, {key: value} given)
    for s in sources:
        s = {"kind": s} if isinstance(s, str) else dict(s)
        cls = _source_class(s["kind"], s.get("module", ""))
        if cls is None:
            raise ParamError("unknown source kind %r; known: %s" % (s["kind"], ", ".join(_sources().source_kinds())))
        srcs.append((cls, s))
    coros = []  # {type: ..., key: value}
    for c in coronagraphs:
        if isinstance(c, str):
            t, _, charge = c.partition(":")
            c = {"type": t or "vortex"}
            if charge:
                c["charge"] = int(charge)
        if c.get("type", "vortex") not in CORONAGRAPHS:
            raise ParamError("unknown coronagraph type %r; known: %s" % (c.get("type"), ", ".join(CORONAGRAPHS)))
        coros.append(c)
    requires = set()
    for cls, _ in srcs:
        requires |= set(cls.requires)
    families = set()
    for cls, s in srcs:
        settings = {p.name: p.default for p in cls.params}
        settings.update(s)
        families.add(_family(cls, settings))
    families.discard(None)
    image_source = any("grid.pixel_scale" in cls.provides for cls, _ in srcs)
    flux_needed = len(families) > 1  # the star's flux puts the components in one unit

    lines = []
    if header:
        if schema:
            lines.append("#:schema ./%s" % SCHEMA_FILE)
        lines += ["# Toy Coronagraph parameters. Hover a key in VS Code (Even Better TOML), or run",
                  "#   python -m toycoronagraph describe KEY     for details. Every key is optional unless noted.",
                  "# Each quantity is set in ONE place; the values used (and where they came from) are saved with the results.",
                  ""]

    def used(key):
        """The kinds of the chosen sources requiring a shared key ("section.name")."""
        return [cls.kind for cls, _ in srcs if key in cls.requires]

    def section(name, keys, note=""):
        """
        Append a core section to lines.

        Args:
            name (str): the section.
            keys (dict): {key: value} written; with full, the section's other keys are written commented out.
            note (str): the comment of the section header; "" = the section's doc.
        """
        doc = SECTIONS[name][0]
        lines.append("%-26s # %s" % ("[%s]" % name, note or doc))
        for p in SECTIONS[name][1]:
            if p.name in keys:
                lines.append(_line(p.name, keys[p.name], _comment(p, used("%s.%s" % (name, p.name))), 26))
            elif full:
                lines.append("# " + _line(p.name, _example_value(p), _comment(p), 24))
        lines.append("")

    lazuli = any(c.get("type") == "lazuli" for c in coros)
    obs = {"wavelength": AUTO if image_source else 1e-6, "diameter": 2.4}
    if lazuli:
        obs["wavelength"] = LAZULI_WAVELENGTH
    section("observation", obs, "the observation; the wavelength: Lazuli's design wavelength (r' band)" if lazuli else "")
    if image_source:
        section("grid", {"size": _odd_image_size(srcs) if lazuli else AUTO, "pixel_scale": AUTO},
                "the simulated image; \"auto\" = the image source's")
    else:
        section("grid", {"sampling": 4.0, "fov": LambdaD(32)}, "the simulated image; star at pixel (N-1)/2")
    star = {}
    for p in SECTIONS["star"][1]:
        key = "star." + p.name
        if key in requires or (flux_needed and p.name in ("temperature", "radius", "distance")):
            star[p.name] = p.example
    if star or full:
        who = sorted({cls.kind for cls, _ in srcs if set(cls.requires) & {"star." + n for n in star}})
        note = ("used by: " + ", ".join(who + (["unit conversion"] if flux_needed else []))) if star else ""
        section("star", star, note)
    for cls, s in srcs:
        lines.append("%-26s # %s" % ("[[source]]", "one table per component; enabled = false switches it off"))
        lines.append(_line("kind", cls.kind, cls.doc, 26))
        for p in SOURCE_COMMON[1:]:  # name, enabled, scale, method, module: when given (module: your own modules)
            if p.name in s:
                lines.append(_line(p.name, s[p.name], _comment(p), 26))
        for p in cls.params:
            if p.name in s:
                lines.append(_line(p.name, s[p.name], _comment(p), 26))
            elif p.default is None:
                lines.append(_line(p.name, p.example, _comment(p) + " (required)", 26))
            elif not p.advanced:
                lines.append(_line(p.name, p.default, _comment(p), 26))
            elif full:
                lines.append("# " + _line(p.name, p.default, _comment(p), 24))
        lines.append("")
    for c in coros:
        doc, params = CORONAGRAPHS[c.get("type", "vortex")]
        lines.append("%-26s # %s" % ("[[coronagraph]]", "one table per coronagraph; each gets its own results"))
        lines.append(_line("type", c.get("type", "vortex"), doc, 26))
        for p in params:
            if p.name in c:
                lines.append(_line(p.name, c[p.name], _comment(p), 26))
            elif not p.advanced:
                lines.append(_line(p.name, p.default, _comment(p), 26))
            elif full:
                lines.append("# " + _line(p.name, p.default, _comment(p), 24))
        lines.append("")
    if len(srcs) > 1 or full:
        section("method", {"combine": "separate"} if len(srcs) > 1 else {})
    if full:
        for name in ("compute", "checks", "run", "display"):
            section(name, {})
    text = "\n".join(lines).rstrip()+"\n"
    if path:
        with open(path, "w") as f:
            f.write(text)
        if schema:
            write_schema(os.path.join(os.path.dirname(os.path.abspath(path)), SCHEMA_FILE))
    return text


def _param_schema(p, prefix=""):
    """
    The JSON schema of one key.

    Args:
        p (Param): the declaration.
        prefix (str): text starting the description, e.g. "(wind) ".

    Returns:
        dict: the schema: the type (or enum) with its bounds, alternatives in "anyOf" ("auto", the
            auto_word, a string with a unit; for a list with a unit, its items), "description" /
            "markdownDescription" (doc, type, default, help) and "default" (unless there is none or
            it is a LambdaD).
    """
    desc = "%s%s\n\n**Type:** %s" % (prefix, p.doc[0].upper()+p.doc[1:] if p.doc else "", p.describe_type())
    if p.default is not None:
        desc += "  \n**Default:** `%s`" % _toml(p.default)
    if p.help:
        desc += "\n\n" + p.help
    if p.choices:
        base = {"enum": list(p.choices)}
    elif p.type == "bool":
        base = {"type": "boolean"}
    elif p.type == "str":
        base = {"type": "string"}
    elif p.type == "table":
        base = {"type": "object"}
    elif p.type.endswith("_list"):
        base = {"type": "array", "items": {"type": {"int_list": "integer", "float_list": "number",
                                                    "str_list": "string"}[p.type]}}
        if p.unit:  # items may also be text with a unit
            base["items"] = {"anyOf": [base["items"], {"type": "string", "pattern": _QUANTITY.pattern}]}
    else:
        base = {"type": "integer" if p.type == "int" else "number"}
        if p.min is not None:
            base["exclusiveMinimum" if p.exclusive_min else "minimum"] = p.min
        if p.max is not None:
            base["maximum"] = p.max
    options = [base]
    if p.auto and not p.choices:
        options.append({"const": AUTO})
    if p.auto_word:
        options.append({"const": p.auto_word})
    if p.unit and p.type == "float":
        options.append({"type": "string", "pattern": _QUANTITY.pattern,
                        "description": "a value with a unit, e.g. \"3.22 pc\""})
    schema = options[0] if len(options) == 1 else {"anyOf": options}
    schema["description"] = desc
    schema["markdownDescription"] = desc
    if p.default is not None and not isinstance(p.default, LambdaD):
        schema["default"] = p.default
    return schema


def _either(first, second):
    """
    The schema of a key declared by several modules (or coronagraph types): any of their schemas.

    Args:
        first (dict): the schema so far (from _param_schema, or an earlier _either).
        second (dict): the next module's schema of the key.

    Returns:
        dict: {"anyOf": [first, second]} with both descriptions.
    """
    return {"anyOf": [first, second], "description": first["description"] + "\n\n" + second["description"],
            "markdownDescription": first.get("markdownDescription", "") + "\n\n" + second["markdownDescription"]}


def json_schema():
    """
    JSON schema of the parameter file (editor hover, completion and checks; e.g. Even Better TOML).

    [[source]] tables accept the keys of every kind (a key shared by several kinds: "anyOf" their
    schemas), then each kind is limited to its own keys; an unknown kind needs `module`.
    [[coronagraph]] tables likewise, each type also checking the values of its own keys (e.g. a
    vortex charge >= 2, a lazuli charge >= 0).

    Returns:
        dict: the schema (draft-07), ready for json.dump.
    """
    props = {"include": {"type": "array", "items": {"type": "string"},
                         "description": "Other parameter files to read first (e.g. a shared star description). "
                                        "A key set in two files is an error."}}
    for name, (doc, params) in SECTIONS.items():
        props[name] = {"type": "object", "description": doc[0].upper()+doc[1:], "additionalProperties": False,
                       "properties": {p.name: _param_schema(p) for p in params}}
    kinds = _sources().source_kinds()
    common = [p.name for p in SOURCE_COMMON]
    union, per_kind = {}, {}  # {key: schema} of all kinds; {kind: its allowed keys}
    for p in SOURCE_COMMON:
        union[p.name] = _param_schema(p)
    kind_doc = ("Which module. Built-in: " + "; ".join("`%s`: %s" % (k, _sources().get_source(k).doc) for k in kinds)
                + ". Your own module: its kind, with `module` set.")
    union["kind"] = {"type": "string", "description": kind_doc, "markdownDescription": kind_doc, "examples": kinds}
    for k in kinds:
        cls = _sources().get_source(k)
        per_kind[k] = common + [p.name for p in cls.params]
        for p in cls.params:
            entry = _param_schema(p, "(%s) " % k)
            union[p.name] = _either(union[p.name], entry) if p.name in union and p.name not in common else entry
    source_item = {"type": "object", "required": ["kind"], "properties": union,
                   "allOf": [{"if": {"properties": {"kind": {"const": k}}, "required": ["kind"]},
                              "then": {"propertyNames": {"enum": per_kind[k]}}} for k in kinds] +
                            [{"if": {"properties": {"kind": {"not": {"enum": kinds}}}, "required": ["kind"]},
                              "then": {"required": ["module"]}}]}
    props["source"] = {"type": "array", "items": source_item,
                       "description": "One table per component of the scene ([[source]]); kind picks the module."}
    cunion, cper = {p.name: _param_schema(p) for p in CORONAGRAPH_COMMON}, {}  # as union, per_kind for the types
    cunion["type"]["enum"] = list(CORONAGRAPHS)
    cown = {}  # {type: {key: schema}}: each type's own value schemas (a key of several types gets the type's own)
    for t, (doc, params) in CORONAGRAPHS.items():
        cper[t] = [p.name for p in CORONAGRAPH_COMMON] + [p.name for p in params]
        cown[t] = {p.name: _param_schema(p, "(%s) " % t) for p in params}
        for name, entry in cown[t].items():
            cunion[name] = _either(cunion[name], entry) if name in cunion else entry
    # a table without type matches every "if": it gets the limits of every type, the vortex's (its default) included
    props["coronagraph"] = {"type": "array", "description": "One table per coronagraph ([[coronagraph]]).",
                            "items": {"type": "object", "properties": cunion,
                                      "allOf": [{"if": {"properties": {"type": {"const": t}}},
                                                 "then": {"propertyNames": {"enum": cper[t]}, "properties": cown[t]}}
                                                for t in CORONAGRAPHS]}}
    return {"$schema": "http://json-schema.org/draft-07/schema#", "title": "Toy Coronagraph parameters",
            "type": "object", "additionalProperties": False, "properties": props}


def write_schema(path):
    """
    Write json_schema() to a file (indented JSON).

    Args:
        path (str): the file (overwritten), normally SCHEMA_FILE next to the parameter file.
    """
    with open(path, "w") as f:
        json.dump(json_schema(), f, indent=1)
        f.write("\n")


def _md_rows(params, prefix=""):
    """
    Markdown table rows "| key | default | unit | meaning |".

    Args:
        params (iterable of Param): the declarations.
        prefix (str): text before each key name.

    Returns:
        list of str: one row per parameter.
    """
    rows = []
    for p in params:
        default = "—" if p.default is None else "`%s`" % _toml(p.default)
        rows.append("| `%s%s` | %s | %s | %s |" % (prefix, p.name, default, p.unit or "", p.doc.replace("|", "\\|")))
    return rows


def reference_markdown():
    """
    Parameter reference tables (README / docs), generated from the declarations.

    Returns:
        str: Markdown: the core sections, [[source]] (common keys, then each kind), [[coronagraph]]
            types, and the rules when values disagree.
    """
    out = ["# Parameter reference", "",
           "Generated from the declarations in `toycoronagraph/params.py` and the source modules. "
           "Every key is optional unless a module needs it; each quantity is set in one place.", ""]
    head = ["| key | default | unit | meaning |", "|---|---|---|---|"]
    for name, (doc, params) in SECTIONS.items():
        out += ["## [%s]: %s" % (name, doc), ""] + head + _md_rows(params) + [""]
    out += ["## [[source]]: the components", "", "Keys every source table may have:", ""] + head + _md_rows(SOURCE_COMMON) + [""]
    for k in _sources().source_kinds():
        cls = _sources().get_source(k)
        extra = []
        if cls.requires:
            extra.append("needs " + ", ".join("`%s`" % r for r in cls.requires))
        if cls.provides:
            extra.append("can supply " + ", ".join("`%s`" % r for r in cls.provides))
        fam = {"physical": "physical units", "relative": "contrast (relative to the star)", None: "from the input"}[cls.family]
        out += ["### kind = \"%s\": %s" % (k, cls.doc), "", "Output: %s. %s" % (fam, "; ".join(extra).capitalize() + "." if extra else ""), ""]
        out += head + _md_rows(cls.params) + [""]
    out += ["## [[coronagraph]]", ""]
    for t, (doc, params) in CORONAGRAPHS.items():
        out += ["### type = \"%s\": %s" % (t, doc), ""] + head + _md_rows(params) + [""]
    out += ["## When values disagree", "", "| quantity | tolerance | consequence | why |", "|---|---|---|---|"]
    for key, (tol, consequence, why) in CONFLICT_RULES.items():
        out.append("| `%s` | %g%% | %s | %s |" % (key, 100*tol, consequence, why))
    out += ["| an image's pixel scale or size vs `[grid]` | — | resampled onto the grid (a note) | flux conserved |",
            "| `grid`: pixel_scale + sampling | 0.1% | error | over-specified |",
            "| `grid`: size x pixel_scale vs fov | 1.5 px (or 0.1%) | error | over-specified |",
            "| `star`: temperature + flux_jy | 0.1% | error | over-specified |", ""]
    return "\n".join(out)


def describe(query):
    """
    Where a key is set, its unit, default and allowed values, and who uses it.

    Args:
        query (str): "section" or "section.key" (e.g. "grid", "star.distance"; brackets are ignored),
            "kind" or "kind.key" of a source module (e.g. "wind", "wind.mdot"), "type" or "type.key" of
            a coronagraph (e.g. "vortex.charge"), or a bare key name searched everywhere (e.g. "charge").

    Returns:
        str: one block per match, separated by blank lines; if nothing matches, "Nothing called ..."
            with the closest names.
    """
    query = query.strip().lstrip("[").rstrip("]")
    head, _, tail = query.partition(".")
    out = []
    kinds = _sources().source_kinds()

    def param_text(p, where, users=None):
        """
        The description block of one key.

        Args:
            p (Param): the declaration.
            where (str): where it is set, e.g. "[star] distance" (also used to find its CONFLICT_RULES).
            users (list of str): source kinds needing it; None or [] = not shown.

        Returns:
            str: lines: where and type, doc, default, help, needed by, conflict rule.
        """
        lines = ["%s  (%s)" % (where, p.describe_type()), "  " + p.doc]
        if p.default is not None:
            lines.append("  default: %s" % _toml(p.default))
        if p.help:
            lines.append("  " + p.help)
        if users:
            lines.append("  needed by: " + ", ".join(users))
        key = where.strip("[]").replace("] ", ".")  # "[star] distance" -> "star.distance"
        if key in CONFLICT_RULES:
            tol, consequence, why = CONFLICT_RULES[key]
            lines.append("  if values disagree (> %g%%): %s; %s" % (100*tol, consequence, why))
        return "\n".join(lines)

    if head in SECTIONS:
        doc, params = SECTIONS[head]
        for p in params:
            if not tail or p.name == tail:
                users = [k for k in kinds if "%s.%s" % (head, p.name) in _sources().get_source(k).requires]
                out.append(param_text(p, "[%s] %s" % (head, p.name), users))
        if not tail:
            out.insert(0, "[%s]: %s" % (head, doc))
    elif head in kinds:
        cls = _sources().get_source(head)
        if not tail:
            block = ['[[source]] kind = "%s": %s' % (head, cls.doc),
                     "  output: %s" % {"physical": "physical units", "relative": "contrast", None: "from the input"}[cls.family]]
            if cls.requires:
                block.append("  needs: " + ", ".join(cls.requires))
            if cls.provides:
                block.append("  can supply: " + ", ".join(cls.provides))
            out.append("\n".join(block))
        # a key: also the common keys (but kind and name); the whole kind: its own keys only
        for p in ((SOURCE_COMMON[2:] + tuple(cls.params)) if tail else cls.params):
            if not tail or p.name == tail:
                out.append(param_text(p, '[[source]] (kind = "%s") %s' % (head, p.name)))
    elif head in CORONAGRAPHS:
        doc, params = CORONAGRAPHS[head]
        if not tail:
            out.append('[[coronagraph]] type = "%s": %s' % (head, doc))
        for p in params:
            if not tail or p.name == tail:
                out.append(param_text(p, '[[coronagraph]] (type = "%s") %s' % (head, p.name)))
    else: # a bare key: search everywhere
        for section, (_, params) in SECTIONS.items():
            for p in params:
                if p.name == query:
                    users = [k for k in kinds if "%s.%s" % (section, p.name) in _sources().get_source(k).requires]
                    out.append(param_text(p, "[%s] %s" % (section, p.name), users))
        for k in kinds:
            for p in _sources().get_source(k).params:
                if p.name == query:
                    out.append(param_text(p, '[[source]] (kind = "%s") %s' % (k, p.name)))
        for t, (_, params) in CORONAGRAPHS.items():
            for p in params:
                if p.name == query:
                    out.append(param_text(p, '[[coronagraph]] (type = "%s") %s' % (t, p.name)))
    if not out:
        names = list(SECTIONS) + kinds + list(CORONAGRAPHS) + sorted(
            {p.name for _, ps in SECTIONS.values() for p in ps} | {p.name for k in kinds for p in _sources().get_source(k).params}
            | {p.name for _, ps in CORONAGRAPHS.values() for p in ps})
        close = difflib.get_close_matches(tail or head, names, 3)
        return "Nothing called %r.%s" % (query, (" Did you mean: " + ", ".join(close) + "?") if close else "")
    return "\n\n".join(out)


def migrate(path, out_path=None):
    """
    Rewrite an old-layout file in the new layout (only the values set there are kept).

    Args:
        path (str): the old file: a .toml file, or a .py parameter file (read by config._from_legacy).
        out_path (str): where to write the new file (and the JSON schema next to it); None = only return it.

    Returns:
        str: the TOML text in the new layout.

    Raises:
        ParamError: the .toml file is already in the new layout, or the mapped parameters are not valid
            (from_data()).
    """
    if path.endswith(".py"):
        from toycoronagraph import config
        cfg = config._from_legacy(path)
        data = {"target": {"target_pixel": cfg.target_pixel}, "telescope": {"wavelength": cfg.wavelength,
                "diameter": cfg.diameter, "visual_range": float(cfg.visual_range)},
                "coronagraph": {"type": cfg.coronagraph_type, "lyot_mask_size": cfg.lyot_mask_size},
                "method": {"rot_number": cfg.rot_number}, "compute": {"num_cores": cfg.num_cores}}
        if cfg.fits_filename:
            data["target"]["fits_filename"] = cfg.fits_filename
    else:
        with open(path, "rb") as f:
            data = tomllib.load(f)
        data.pop("resolved", None)  # values written back by old runs: not parameters
        if not is_legacy(data):
            raise ParamError("%s is already in the new layout" % path)
    new = legacy_to_new(data, os.path.dirname(os.path.abspath(path)))
    from_data(new)  # check it
    lines = ["#:schema ./%s" % SCHEMA_FILE,
             "# Migrated from %s (%s) to the new layout. Only the values set there are kept;"
             % (os.path.basename(path), datetime.date.today()),
             "# everything else uses the defaults (python -m toycoronagraph describe KEY).", ""]
    if "fov" in new.get("grid", {}):  # inserted before the blank line ending the header
        lines[-1:-1] = ["# Note: [grid] fov sets the simulated field. An image whose header gives its pixel scale is",
                        "# resampled onto it; the old layout stretched such an image to fit visual_range (with a warning)."]
    for name in ("observation", "grid", "star"):
        if new.get(name):
            lines.append("[%s]" % name)
            for p in SECTIONS[name][1]:
                if p.name in new[name]:
                    lines.append(_line(p.name, new[name][p.name], _comment(p), 26))
            lines.append("")
    for s in new["source"]:
        cls = _source_class(s["kind"])
        lines.append("[[source]]")
        lines.append(_line("kind", s["kind"], cls.doc, 26))
        for p in SOURCE_COMMON[1:] + tuple(cls.params):
            if p.name in s:
                lines.append(_line(p.name, s[p.name], _comment(p), 26))
        lines.append("")
    for c in new["coronagraph"]:
        ctype = c.get("type", "vortex")
        lines.append("[[coronagraph]]")
        lines.append(_line("type", ctype, CORONAGRAPHS[ctype][0], 26))
        for p in CORONAGRAPHS[ctype][1]:
            if p.name in c:
                lines.append(_line(p.name, c[p.name], _comment(p), 26))
        lines.append("")
    for name in ("method", "compute", "display"):
        if new.get(name):
            lines.append("[%s]" % name)
            for p in SECTIONS[name][1]:
                if p.name in new[name]:
                    lines.append(_line(p.name, new[name][p.name], _comment(p), 26))
            lines.append("")
    text = "\n".join(lines).rstrip()+"\n"
    if out_path:
        with open(out_path, "w") as f:
            f.write(text)
        write_schema(os.path.join(os.path.dirname(os.path.abspath(out_path)), SCHEMA_FILE))
    return text
