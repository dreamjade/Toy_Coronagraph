"""
Source modules: each input of a scene (an image, an MCFOST model, a stellar wind, a CME,
planets, or your own module) is one class following the Source interface.

A module declares, as class attributes:
    kind      the name used in the parameter file ([[source]] kind = "...")
    doc       one line describing it
    family    "physical" (e.g. Jy/arcsec^2), "relative" (contrast: a fraction of the star's
              flux) or None (decided by the input, e.g. a FITS file's BUNIT)
    params    its own parameters (toycoronagraph.params.Param); shared quantities such as
              the star's distance are never repeated here, they are listed in `requires`
    requires  the shared quantities it needs, e.g. ("star.distance", "star.radius")
    provides  the shared quantities it can supply, e.g. ("grid.pixel_scale",) from a header

and implements render(grid, star), returning its image centred on the star (pixel (N-1)/2 of
the grid, or of its own grid with pixel_scale set); a module that wraps another code makes
sure of this itself. Register a module with @register_source; a module outside the package
is found from the parameter file's `module = "my_package.my_module"` (imported, which
registers it). tests/extra_sources.py (CmeBlob) is a documented example.
"""
import importlib

_SOURCES = {}  # registry: kind -> Source subclass (filled by @register_source)


class Source:
    """
    Base class of the source modules (see the module docstring).

    Attributes:
        kind (str): name in parameter files; must be set by every registered subclass.
        doc (str): one-line description.
        family (str): "physical", "relative", or None when the input decides.
        params (tuple of Param): the module's own parameters.
        requires (tuple of str): shared quantities it needs, e.g. "star.distance".
        provides (tuple of str): shared quantities it can supply (reported by inspect()).
        settings (dict): {param name: value} of every parameter, defaults filled in (set by __init__).
    """
    kind = None
    doc = ""
    family = None
    params = ()
    requires = ()
    provides = ()

    def __init__(self, **settings):
        """
        Check and store the module's own parameters.

        Args:
            **settings: {param name: value}; each is converted with toycoronagraph.params.convert
                (units in text, "auto", choices...); parameters not given take their default.

        Raises:
            ParamError: unknown parameter names (with a suggestion) or invalid values, all listed.
        """
        from toycoronagraph.params import convert, ParamError
        import difflib
        known = {p.name: p for p in self.params}
        errors, values = [], {}
        for key, value in settings.items():
            if key not in known:
                close = difflib.get_close_matches(key, list(known), 1, 0.6)
                errors.append("%s: unknown parameter %s%s" % (self.kind, key, (" (did you mean %s?)" % close[0]) if close
                                                              else "; known: " + ", ".join(known)))
                continue
            try:
                values[key] = convert(known[key], value, self.kind)
            except ValueError as err:
                errors.append(str(err))
        for p in self.params:
            if p.name not in values:
                values[p.name] = p.default
        if errors:
            raise ParamError("\n".join(errors))
        self.settings = values

    def missing(self):
        """
        Required parameters (no default) that are not set.

        Returns:
            list of Param: empty when the module can be rendered.
        """
        return [p for p in self.params if self.settings.get(p.name) is None and p.default is None]

    def family_of(self):
        """
        The family of this source's images, given its settings.

        Returns:
            str: "physical", "relative", or None when only the input decides (e.g. a file's BUNIT).
        """
        return self.family

    def inspect(self):
        """
        Shared quantities this source supplies (e.g. from a file header); see `provides`.

        Returns:
            list of tuple: [(key, value, origin)], e.g. ("grid.pixel_scale", 0.01, "header PIXSCALE of x.fits");
                empty by default.
        """
        return []

    def point_sources(self, grid, star):
        """
        Unresolved sources, for modules that are not images (render() is then not called).

        Args:
            grid (toycoronagraph.scene.Grid): the scene grid.
            star (toycoronagraph.scene.Star): the star.

        Returns:
            list or None: [((radius arcsec, angle rad), flux)], the flux in Jy (family "physical")
                or as a fraction of the star's flux ("relative"); None for image modules (default).
        """
        return None

    def render(self, grid, star):
        """
        The image(s) of this source; every image module overrides it.

        Args:
            grid (toycoronagraph.scene.Grid): the scene grid (size, pixel_scale arcsec, wavelength m).
            star (toycoronagraph.scene.Star): the star (distance pc, radius R_sun, temperature K, flux_jy).

        Returns:
            SourceImage or dict: a units.SourceImage centred on the star, or {name: SourceImage}
                for several components.

        Raises:
            NotImplementedError: in the base class.
        """
        raise NotImplementedError("%s: rendering is added in a later stage" % self.kind)

    def symmetric(self):
        """
        Whether the image is circularly symmetric, from the "symmetric" setting if the module has one.

        Returns:
            bool or None: True / False if known (by construction or set by the user); None: measure it.
        """
        value = self.settings.get("symmetric", "auto")
        return None if value == "auto" else bool(value)


def register_source(cls):
    """
    Class decorator: make a Source subclass available under its kind.

    Args:
        cls (type): the Source subclass; registering a kind again replaces the previous class.

    Returns:
        type: cls itself.

    Raises:
        ValueError: no kind, or two parameters with the same name.
    """
    if not cls.kind:
        raise ValueError("%s has no kind" % cls.__name__)
    names = [p.name for p in cls.params]
    if len(set(names)) != len(names):
        raise ValueError("%s: duplicate parameter names" % cls.kind)
    _SOURCES[cls.kind] = cls
    return cls


def get_source(kind, module=None):
    """
    The Source class of a kind.

    Args:
        kind (str): the registered kind, e.g. "wind".
        module (str): Python module to import first if the kind is not registered yet
            (importing it registers its classes); None = none.

    Returns:
        type: the Source subclass.

    Raises:
        KeyError: the kind is still not registered.
    """
    _load_builtins()
    if kind not in _SOURCES and module:
        importlib.import_module(module)
    if kind not in _SOURCES:
        raise KeyError(kind)
    return _SOURCES[kind]


BUILTIN_KINDS = ("fits_image", "mcfost", "scene_fits", "wind", "cme", "planet")


def source_kinds():
    """
    All registered kinds.

    Returns:
        list of str: the built-in ones first (in the order of BUILTIN_KINDS), then yours (in
            registration order).
    """
    _load_builtins()
    return [k for k in BUILTIN_KINDS if k in _SOURCES] + [k for k in _SOURCES if k not in BUILTIN_KINDS]


def _load_builtins():
    """Import the built-in source modules, which registers them (imported on first use, not with the package)."""
    from toycoronagraph.sources import images, thomson, cme, planets  # noqa: F401 (registers them)
