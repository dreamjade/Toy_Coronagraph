"""
The scene: the image before the coronagraph, built from source modules.

    scene = Scene.from_file("eps_eri.toml")                     # from a parameter file
    scene = Scene(wavelength=500e-9, diameter=6.0, sampling=4, fov=2.0,
                  star=Star(distance=3.22, radius=0.73))          # or in Python
    scene.add(Wind(mdot=6e-13))                                  # any Source module, yours too
    images = scene.render()                                      # {name: image}, and scene.total

Each component is rendered by its module (centred on the star, on its own grid, in its own
unit), put per pixel (Jy, or the fraction of the star's flux), resampled onto the scene grid
with the flux conserved, and converted to the scene's internal unit (units.choose_internal):
"Jy/arcsec^2" when the star's flux is known (or every component is physical), else
"contrast_pixel" (the fraction of the star's flux per pixel). Point sources (planets) are
not images: they are kept as positions and fluxes (Scene.points).

Conventions: images are (N, N) float64 arrays with the star at pixel ((N-1)/2, (N-1)/2);
rows are +y (upwards, as in FITS; plotted with origin="lower"), columns +x; angles are
measured from +x towards +y. The shared quantities (grid, wavelength, star) are resolved once
by toycoronagraph.params.

Also here: the plotting helpers shared with toycoronagraph.observe (show_image, mark_points,
point_labels, blank_panel, radial_profile, plot_profiles) and the scene FITS writer
(write_scene, image_hdu).
"""
from dataclasses import dataclass
import numpy as np
from astropy.io import fits
from toycoronagraph import params as P
from toycoronagraph import units as U


@dataclass(frozen=True)
class Grid:
    """
    The simulated image: size x size pixels of pixel_scale arcsec, star at pixel (size-1)/2.

    Attributes:
        size (int): Pixels on one edge of the (square) image.
        pixel_scale (float): Pixel size in arcsec.
        wavelength (float): Observing wavelength in m.
        diameter (float): Telescope aperture diameter D in m.
    """
    size: int
    pixel_scale: float     # arcsec
    wavelength: float      # m
    diameter: float        # m

    @property
    def lambda_d(self):
        """float: lambda/D in arcsec."""
        return self.wavelength/self.diameter*P.RAD_TO_ARCSEC

    @property
    def sampling(self):
        """float: Pixels per lambda/D."""
        return self.lambda_d/self.pixel_scale

    @property
    def visual_range(self):
        """float: Image width in lambda/D (twice the Optics psf_range)."""
        return self.size*self.pixel_scale/self.lambda_d

    @property
    def pixel_area(self):
        """float: Area of one pixel in arcsec^2."""
        return self.pixel_scale**2

    @property
    def peak_norm(self):
        """float: The peak of an ideal circular aperture's PSF (no coronagraph, no Lyot stop) as a fraction
        of the star's light per pixel, pi/4 (pixel_scale / lambda_d)^2: the contrast_peak reference of
        input images (final images use their coronagraph's own, Coronagraph.peak_norm)."""
        return np.pi/4/self.sampling**2

    @property
    def extent(self):
        """list: [left, right, bottom, top] pixel edges in arcsec from the star (for imshow's extent)."""
        h = self.size/2*self.pixel_scale
        return [-h, h, -h, h]

    def coordinates(self):
        """
        Offsets from the star of the pixel centres.

        Returns:
            tuple: (x, y), each a (size, size) float64 array in arcsec; rows are y (upwards, as in
                FITS), columns x (np.meshgrid order: x varies along axis 1, y along axis 0).
        """
        offsets = (np.arange(self.size)-(self.size-1)/2)*self.pixel_scale
        return np.meshgrid(offsets, offsets)

    def pixel(self, radius, angle):
        """
        Pixel position of a point given in polar coordinates around the star.

        Args:
            radius (float or np.ndarray): Distance from the star in arcsec.
            angle (float or np.ndarray): Angle in rad, from +x towards +y.

        Returns:
            tuple: (row, column) as floats (not rounded; arrays if the inputs are), rows = +y.
        """
        c = (self.size-1)/2
        return c + radius*np.sin(angle)/self.pixel_scale, c + radius*np.cos(angle)/self.pixel_scale

    def axis(self, unit="arcsec", distance=None):
        """
        Scale factor and label for plot axes in another unit.

        Args:
            unit (str): "arcsec", "mas", "lambda/D", "au" or "pixel".
            distance (float): The star's distance in pc; needed for "au" (1 arcsec at d pc is d au).

        Returns:
            tuple: (factor, label): multiply arcsec by factor (float) to get `unit`; label is `unit` (str).

        Raises:
            ParamError: "au" without a distance, or an unknown unit.
        """
        if unit == "au" and not distance:
            raise P.ParamError('axis = "au" needs the distance: set [star] distance')
        factors = {"arcsec": 1.0, "mas": 1e3, "lambda/D": 1/self.lambda_d, "au": distance or 0.0, "pixel": 1/self.pixel_scale}
        if unit not in factors:
            raise P.ParamError("axis = %r: use one of %s" % (unit, ", ".join(factors)))
        return factors[unit], unit


@dataclass(frozen=True)
class Star:
    """
    The star / system. Any field may be None (unknown).

    Attributes:
        distance (float): Distance in pc.
        radius (float): Stellar radius in R_sun.
        temperature (float): Effective temperature in K (blackbody flux, with radius and distance).
        flux_jy (float): Flux density in Jy at the scene's wavelength (overrides the blackbody).
    """
    distance: float = None
    radius: float = None
    temperature: float = None
    flux_jy: float = None

    def flux_at(self, wavelength):
        """
        The star's flux density at a wavelength.

        Args:
            wavelength (float): Wavelength in m.

        Returns:
            float or None: flux_jy if set, else the blackbody flux in Jy (needs temperature, radius
                and distance); None if neither is known.
        """
        if self.flux_jy is not None:
            return self.flux_jy
        if None not in (self.temperature, self.radius, self.distance):
            return P.blackbody_flux_jy(self.temperature, self.radius, self.distance, wavelength)
        return None


def resample(values, scale_in, scale_out, size_out):
    """
    Resample per-pixel values onto another grid centred on the same star, conserving the flux.

    The values must be per pixel (fluxes, not surface brightnesses): each input pixel's flux
    is shared among the output pixels it overlaps, in proportion to the overlap area (integer
    factors give block sums). The 2-D overlap is separable, so it is one 1-D weight matrix w
    applied to both axes: w @ values @ w.T. Light outside the output grid is dropped.

    Args:
        values (np.ndarray): (n_in, n_in) per-pixel values, star at pixel (n_in-1)/2.
        scale_in (float): Input pixel size in arcsec.
        scale_out (float): Output pixel size in arcsec.
        size_out (int): Output pixels on one edge.

    Returns:
        np.ndarray: (size_out, size_out) float64 per-pixel values, star at pixel (size_out-1)/2.
    """
    n_in = values.shape[0]
    e_in = (np.arange(n_in+1)-n_in/2)*scale_in        # pixel edges, from the star
    e_out = (np.arange(size_out+1)-size_out/2)*scale_out
    overlap = np.minimum(e_out[1:, None], e_in[None, 1:]) - np.maximum(e_out[:-1, None], e_in[None, :-1])
    w = np.clip(overlap, 0, None)/scale_in          # (size_out, n_in): share of each input pixel
    return w @ values @ w.T


class Scene:
    """
    Components (source modules) on one grid, in one unit; see the module docstring.

    Attributes:
        sections (dict): Parameter sections {section: {key: value}} given by the user.
        components (list): The params.Component of each source module (enabled or not).
        coronagraphs (list): [[coronagraph]] tables of the parameters (used by observe()).
        params (Params or None): Checked parameters; built from the above by build() (None until then,
            and reset by add()).

    Set by build():
        grid (Grid): The scene grid.
        star (Star): The resolved star.

    Set by render():
        images (dict or None): {component name: (N, N) float64 image in internal_unit}; None until
            rendered (and after add()).
        points (dict): Point-source components {name: (positions, fluxes)}: positions a list of
            (radius arcsec, angle rad), fluxes a 1-D array of each source's flux per pixel area, in
            internal_unit (multiply by grid.pixel_area for Jy, or for the fraction of the star's flux).
        total (np.ndarray): (N, N) sum of `images` (zeros if there is none; point sources are not in it).
        internal_unit (str): "Jy/arcsec^2" or "contrast_pixel".
        star_flux (float or None): The star's flux density in Jy at the wavelength (None: unknown).
        families (dict): {name: "physical" or "relative"} of every image and point-source component.
        units (dict): {name: the unit the module gave its image in} (images only).
        origins (dict): {name: where the image came from, e.g. a file name} (images only).
        component_of (dict): {name: the params.Component it came from} (images and point sources).
        notes (list): Notes and warnings (str) from the parameters and the rendering.
    """

    def __init__(self, star=None, **settings):
        """
        Args:
            star (Star or dict): The star (distance pc, radius R_sun, temperature K or flux_jy Jy);
                None fields are left out (so the sources or defaults may provide them).
            settings: [observation] / [grid] keys (wavelength, diameter, size, pixel_scale,
                sampling, fov, ...), and other sections as dicts (method={...}, checks={...}, ...).

        Raises:
            ParamError: A setting that is neither an [observation] / [grid] key nor a section given as a dict.
        """
        sections = {}
        for key, value in settings.items():
            # a plain key belongs to whichever of [observation] / [grid] declares it
            section = next((s for s in ("observation", "grid") if key in {p.name for p in P.SECTIONS[s][1]}), None)
            if section:
                sections.setdefault(section, {})[key] = value
            elif key in P.SECTIONS and isinstance(value, dict):
                sections[key] = dict(value)
            else:
                raise P.ParamError("Scene: unknown setting %s" % key)
        if star is not None:
            star = star.__dict__ if isinstance(star, Star) else dict(star)
            sections["star"] = {k: v for k, v in star.items() if v is not None}
        self.sections = sections
        self.components = []
        self.coronagraphs = []
        self.params = None      # built from the sections and components when needed
        self.images = None
        self.path, self._origins = None, {}  # the parameter file (from_params), for messages and provenance

    @classmethod
    def from_file(cls, path):
        """
        A scene from a parameter file (see toycoronagraph.params).

        Args:
            path (str): The .toml parameter file.

        Returns:
            Scene: Not rendered yet.
        """
        return cls.from_params(P.load(path))

    @classmethod
    def from_params(cls, p):
        """
        A scene from loaded parameters.

        Args:
            p (Params): From toycoronagraph.params.load (its sections, components and coronagraphs are used).

        Returns:
            Scene: Not rendered yet.
        """
        scene = cls()
        scene.sections, scene.components, scene.coronagraphs, scene.params = p.sections, p.components, p.coronagraphs, p
        scene.path, scene._origins = p.path, p._origins  # kept if components are added later (params rebuilt)
        return scene

    def add(self, source, name=P.AUTO, scale=1.0, enabled=True, method=P.AUTO):
        """
        Add a component; the parameters and the images are built again when next needed.

        Args:
            source (Source): A source module instance (toycoronagraph.sources).
            name (str): Component name (FITS extension, plots); "auto": the module's kind (ring, ring_2, ...).
            scale (float): Factor applied to the component's image or fluxes.
            enabled (bool): Whether the component is rendered.
            method (str): "auto", "circular" or "spatial": this component's method in observe().

        Returns:
            Source: `source` (to keep a handle on it).
        """
        self.components.append(P.Component(source, name, scale, enabled, method))
        self.params = self.images = None
        return source

    def build(self):
        """
        Resolve the grid, wavelength and star with what the sources provide (every disagreement is checked).

        Sets self.params (if None), self.grid and self.star.

        Raises:
            ParamError: Conflicting or missing parameters (from Params.resolve).
        """
        if self.params is None:
            tables = [dict(c.values, type=c.type, name=c.name) for c in self.coronagraphs]
            self.params = P.make(self.sections, self.components, tables, path=self.path, origins=self._origins)
        p = self.params
        p.provided = [item for c in p.enabled() for item in c.source.inspect()]
        p.resolve()
        self.grid = Grid(int(p.get("grid.size")), float(p.get("grid.pixel_scale")), float(p.get("observation.wavelength")),
                         float(p.get("observation.diameter")))
        self.star = Star(**{k: p.get("star."+k) for k in ("distance", "radius", "temperature", "flux_jy")})

    def render(self, verbose=True):
        """
        Render every enabled component onto the grid, in the internal unit.

        Each module's image is put per pixel (units.per_pixel), resampled onto the grid if it
        has another pixel size or field ([grid] resample), scaled by the component's scale and
        converted to the internal unit. Point-source modules give positions and fluxes instead.

        Args:
            verbose (bool): Print the notes and a one-line summary of the scene.

        Returns:
            dict: {component name: (N, N) float64 image in self.internal_unit}; also sets self.total,
                self.points, self.internal_unit, self.star_flux, self.families, self.units,
                self.origins, self.component_of and self.notes (see the class docstring).

        Raises:
            ParamError: An image covers another field than the grid while [grid] resample = false.
            UnitError: The components cannot be put in one unit (units.choose_internal).
        """
        self.build()
        g = self.grid
        self.notes = list(self.params.notes)
        per_pixel, self.families, self.units, self.origins, self.component_of = {}, {}, {}, {}, {}
        points = {}  # point-source components: name -> (positions, fluxes per pixel in their family)
        mismatches = []
        for c in self.params.enabled():
            sources = c.source.point_sources(g, self.star)
            if sources is not None:
                points[c.name] = ([pos for pos, _ in sources], np.array([f for _, f in sources])*c.scale)
                self.families[c.name], self.component_of[c.name] = c.source.family_of(), c
                continue
            out = c.source.render(g, self.star)
            for part, img in (out.items() if isinstance(out, dict) else [(None, out)]):
                # several images from one source are named after them (prefixed if you named the source)
                name = c.name if part is None else part if c.name == c.source.kind else "%s_%s" % (c.name, part)
                img.pixel_scale = img.pixel_scale or g.pixel_scale  # None: already on the scene grid
                self.notes += ["Note: %s: %s" % (name, n) for n in img.notes]
                fam, values = U.per_pixel(img, g.wavelength)
                # differences below 0.01 grid pixel (at the edge) are round-off
                n_in, field_in = values.shape[0], values.shape[0]*img.pixel_scale
                same_field = abs(field_in - g.size*g.pixel_scale) < 0.01*g.pixel_scale
                same_pixels = n_in == g.size and same_field
                if not same_pixels:
                    text = "%s: %d px of %.6g arcsec (%.6g arcsec across), the grid is %d px of %.6g arcsec (%.6g across)" % (
                        img.origin or name, n_in, img.pixel_scale, field_in, g.size, g.pixel_scale, g.size*g.pixel_scale)
                    # another resolution of the same field is always rebinned; another field needs [grid] resample
                    if not (same_field or self.params.get("grid.resample")):
                        mismatches.append(text)
                        continue
                    self.notes.append("Note: %s (flux conserved) %s" % ("rebinned to the grid's resolution, same field:"
                                                                         if same_field else "resampled", text))
                    values = resample(values, img.pixel_scale, g.pixel_scale, g.size)
                per_pixel[name], self.families[name] = values*c.scale, fam
                self.units[name], self.origins[name], self.component_of[name] = img.unit, img.origin, c
        if mismatches:
            raise P.ParamError("[grid] resample = false, but these images cover another field than the grid:\n  "
                               + "\n  ".join(mismatches) + "\n  give the grid their field (for another resolution set only "
                               "[grid] size, with pixel_scale = \"auto\"), or set [grid] resample = true")
        self.star_flux = self.star.flux_at(g.wavelength)
        self.internal_unit = U.choose_internal(self.families, self.star_flux)
        self.images = {n: U.to_internal(self.families[n], v, self.internal_unit, g.pixel_area, self.star_flux)
                       for n, v in per_pixel.items()}
        # point sources: positions, and each one's flux per pixel area in the internal unit
        self.points = {n: (pos, U.to_internal(self.families[n], f, self.internal_unit, g.pixel_area, self.star_flux))
                       for n, (pos, f) in points.items()}
        self.total = sum(self.images.values()) if self.images else np.zeros((g.size, g.size))
        if verbose:
            for note in self.notes:
                print(note)
            print("Scene: %d x %d px, %.5g arcsec/pixel (%.3g px per lambda/D at %.4g m); %s; in %s"
                  % (g.size, g.size, g.pixel_scale, g.sampling, g.wavelength, ", ".join(list(self.images)+list(self.points)),
                     self.internal_unit))
        return self.images

    def in_unit(self, image, unit, star_flux=None):
        """
        An image of this scene in another output unit (see units.convert_image).

        Args:
            image (np.ndarray): (N, N) image in self.internal_unit.
            unit (str): "Jy/arcsec^2", "Jy/pixel", "contrast_pixel" or "contrast_peak" (relative to the
                ideal aperture's PSF peak, Grid.peak_norm); None or "auto": unchanged.
            star_flux (float or Star): The star's flux in Jy, or a Star (its flux at the scene's
                wavelength), to convert a contrast scene to a brightness (and back) without
                rendering again; None: self.star_flux.

        Returns:
            np.ndarray: (N, N) image in `unit` (the same array if no conversion is needed).

        Raises:
            UnitError: The conversion needs the star's flux and it is unknown.
        """
        if unit in (None, P.AUTO):
            return image
        flux = U.star_flux_of(star_flux, self.grid.wavelength) if star_flux is not None else self.star_flux
        return U.convert_image(image, self.internal_unit, unit, self.grid.pixel_area, flux, self.grid.peak_norm)

    def plot(self, log=None, unit=None, axis=None, path=None, show=True, dpi=None, log_range=None, star_flux=None,
             **kwargs):
        """
        Plot each component and the total side by side (rendering the scene first if needed).

        The total panel is added when there are several images or point sources; point sources
        are marked on it with their flux.

        Args:
            log (bool): Log colour scale. None here and below: the [display] setting (see display_setting).
            unit (str): Output unit (see in_unit); "auto": the internal unit.
            axis (str): Axis unit (see Grid.axis).
            path (str): Save the figure there; None: not saved.
            show (bool): Show the figure (else it is closed).
            dpi (int): Resolution of the saved file.
            log_range (float): Log scale: decades shown below each panel's peak.
            star_flux (float or Star): See in_unit.
            **kwargs: Passed to matplotlib's imshow for every panel, e.g. cmap="magma".

        Returns:
            matplotlib.figure.Figure: The figure.
        """
        import matplotlib.pyplot as plt
        if self.images is None:
            self.render(verbose=False)
        log, unit, axis, log_range = (display_setting(self.params, key, value) for key, value in
                                      (("log", log), ("unit", unit), ("axis", axis), ("log_range", log_range)))
        unit = self.internal_unit if unit == P.AUTO else unit
        panels = list(self.images.items()) + ([("total", self.total)] if len(self.images) > 1 or self.points else [])
        fig, axes = plt.subplots(1, len(panels), figsize=(4.2*len(panels), 3.6), squeeze=False)
        for ax, (name, img) in zip(axes[0], panels):
            factor = show_image(ax, fig, self.in_unit(img, unit, star_flux), self.grid, axis, self.star.distance, log,
                                log_range, unit, **kwargs)
            if name == "total":
                mark_points(ax, point_labels(self.points, self.internal_unit, self.grid.pixel_area, unit,
                                             U.star_flux_of(star_flux, self.grid.wavelength) if star_flux is not None
                                             else self.star_flux), factor)
            ax.set_title(name)
        fig.tight_layout()
        return finish(fig, path, show, display_setting(self.params, "dpi", dpi))

    def plot_profiles(self, log=None, unit=None, axis=None, path=None, show=True, dpi=None, star_flux=None, extra=None,
                      log_range=None, **kwargs):
        """
        Azimuthal mean of each component (and of the total, if several) against the distance from the star.

        Args:
            log (bool): Log-log axes (else linear, both from 0). None here and below: the [display] setting.
            unit (str): Output unit (see in_unit); "auto": the internal unit.
            axis (str): Distance unit (see Grid.axis).
            path (str): Save the figure there; None: not saved.
            show (bool): Show the figure (else it is closed).
            dpi (int): Resolution of the saved file.
            star_flux (float or Star): See in_unit.
            extra (dict): More curves {label: (radius in arcsec, values in `unit`)}, dashed.
            log_range (float): Log axes show log_range + PROFILE_DECADES decades below the highest point.
            **kwargs: Passed to matplotlib's plot for every curve, e.g. lw=2.

        Returns:
            matplotlib.figure.Figure: The figure.
        """
        import matplotlib.pyplot as plt
        if self.images is None:
            self.render(verbose=False)
        log, unit, axis, log_range = (display_setting(self.params, key, value) for key, value in
                                      (("log", log), ("unit", unit), ("axis", axis), ("log_range", log_range)))
        unit = self.internal_unit if unit == P.AUTO else unit
        curves = {n: self.in_unit(img, unit, star_flux) for n, img in self.images.items()}
        if len(curves) > 1:
            curves["total"] = self.in_unit(self.total, unit, star_flux)
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        plot_profiles(ax, curves, self.grid, axis, self.star.distance, log, unit, extra,
                      log_range=log_range+PROFILE_DECADES, **kwargs)
        fig.tight_layout()
        return finish(fig, path, show, display_setting(self.params, "dpi", dpi))

    def write(self, path, overwrite=False):
        """
        Save the components as a scene FITS file (readable with kind = "scene_fits"; see write_scene).

        Point sources are not written (Result.save adds them as a POINTS table).

        Args:
            path (str): The FITS file to write.
            overwrite (bool): Replace an existing file (else astropy raises OSError).
        """
        if self.images is None:
            self.render(verbose=False)
        write_scene(path, self.images, self.internal_unit, self.grid, self.star, overwrite=overwrite)


# Profiles show this many decades more than images: curves span more than a colour scale can show
PROFILE_DECADES = 3


def display_setting(params, key, value):
    """
    A plot setting: the value given to the plot, else [display] <key> of the parameters.

    Args:
        params (Params or None): The scene's parameters; None: the declared [display] defaults.
        key (str): The [display] key: "unit", "log", "log_range", "axis" or "dpi".
        value: The value given to the plot; None: not given.

    Returns:
        The setting.
    """
    if value is not None:
        return value
    return params.get("display." + key) if params is not None else P._param("display." + key).default


def finish(fig, path, show, dpi):
    """
    The end of every plot: save the figure (if path is given), then show it or close it.

    Args:
        fig (matplotlib.figure.Figure): The figure.
        path (str): Save the figure there; None: not saved.
        show (bool): Show the figure (else it is closed).
        dpi (int): Resolution of the saved file.

    Returns:
        matplotlib.figure.Figure: The figure.
    """
    import matplotlib.pyplot as plt
    if path:
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.show() if show else plt.close(fig)
    return fig


def show_image(ax, fig, data, grid, axis="arcsec", distance=None, log=False, log_range=5, label="", cmap="gnuplot",
               vmin=None, vmax=None, **kwargs):
    """
    Show an image on its grid (origin="lower": rows are +y), with a colour bar.

    Log scale: log10 of the positive pixels (others are blank, drawn in the darkest colour),
    from the 99.99th percentile (a lone extreme pixel, e.g. a wind's star pixel, does not set
    the scale) down `log_range` decades. Linear scale of non-negative data: from 0.

    Args:
        ax (matplotlib.axes.Axes): Where to draw.
        fig (matplotlib.figure.Figure): The figure (for the colour bar).
        data (np.ndarray): (N, N) image on `grid`.
        grid (Grid): The image's grid.
        axis (str): Axis unit (see Grid.axis).
        distance (float): The star's distance in pc (for axis "au").
        log (bool): Log colour scale.
        log_range (float): Log scale: decades shown below the top.
        label (str): Colour bar label (e.g. the unit); "log10 " is put in front on a log scale.
        cmap (str or Colormap): The colour map.
        vmin (float): Lower colour limit (log10 of the value on a log scale); None: automatic.
        vmax (float): Upper colour limit (log10 of the value on a log scale); None: automatic.
        **kwargs: Passed to matplotlib's imshow, e.g. interpolation="nearest".

    Returns:
        float: The axis factor (arcsec -> `axis` units), e.g. for mark_points.
    """
    factor, name = grid.axis(axis, distance)
    if log:
        data = np.log10(np.where(data > 0, data, np.nan))
        if np.isfinite(data).any():
            top = np.nanpercentile(data, 99.99)
            vmax = top if vmax is None else vmax
            vmin = vmax-log_range if vmin is None else vmin
        label = "log10 " + label
    elif vmin is None and np.nanmin(data) >= 0:
        vmin = 0
    import matplotlib
    cmap = matplotlib.colormaps[cmap] if isinstance(cmap, str) else cmap
    cmap = cmap.with_extremes(bad=cmap(0.0))  # zero (log scale: no light) in the darkest colour, not blank
    im = ax.imshow(data, origin="lower", extent=[e*factor for e in grid.extent], cmap=cmap, vmin=vmin, vmax=vmax, **kwargs)
    ax.set(xlabel="x [%s]" % name, ylabel="y [%s]" % name)
    fig.colorbar(im, ax=ax, label=label, fraction=0.046)
    return factor


def point_labels(points, unit, pixel_area, shown=None, star_flux=None):
    """
    Labels of point sources for mark_points: each source's flux in Jy, or as a fraction of the star's.

    Args:
        points (dict): {name: (positions, fluxes)}: positions a list of (radius arcsec, angle rad),
            fluxes each source's flux per pixel area in `unit` (scalar or 1-D array).
        unit (str): The internal unit, "Jy/arcsec^2" or "contrast_pixel".
        pixel_area (float): Pixel area in arcsec^2.
        shown (str): The unit of the plot: a contrast unit labels fractions of the star's flux, a Jy unit
            fluxes in Jy (converted with star_flux when needed and known); None: `unit`'s.
        star_flux (float): The star's flux in Jy (None: unknown; the labels then stay in `unit`'s family).

    Returns:
        dict: {name: (positions, labels)}, labels a list of str such as "planet 5e-05 Jy".
    """
    jy = unit == "Jy/arcsec^2"
    wanted = jy if shown is None or not star_flux else not shown.startswith("contrast")
    factor = 1.0 if wanted == jy else (1/star_flux if jy else star_flux)
    out = {}
    for name, (positions, fluxes) in points.items():
        total = np.atleast_1d(fluxes)*U.flux_per_pixel(unit, pixel_area)*factor  # per source
        text = ["%s %.2g %s" % (name, f, "Jy" if wanted else "of the star") for f in total]
        out[name] = (positions, text)
    return out


def mark_points(ax, points, factor=1.0):
    """
    Mark point sources with a cyan cross and a label.

    Args:
        ax (matplotlib.axes.Axes): Where to draw (axes in arcsec times `factor`, y upwards).
        points (dict): {name: (positions, labels)}: positions a list of (radius arcsec, angle rad
            from +x towards +y), labels a list of str (see point_labels).
        factor (float): Axis factor from arcsec (see Grid.axis / show_image).
    """
    for name, (positions, labels) in points.items():
        for (r, theta), label in zip(positions, labels):
            x, y = r*np.cos(theta)*factor, r*np.sin(theta)*factor
            ax.plot(x, y, "+", color="cyan", ms=14, mew=2)
            right = x > 0  # keep the label inside the frame
            ax.annotate(label, (x, y), color="cyan", fontsize=8, xytext=(-6 if right else 6, 6), textcoords="offset points",
                        ha="right" if right else "left")


def blank_panel(ax, grid, axis="arcsec", distance=None):
    """
    An empty (black) panel on the grid's axes, e.g. for point sources, which have no image.

    Args:
        ax (matplotlib.axes.Axes): Where to draw.
        grid (Grid): The grid giving the axis limits.
        axis (str): Axis unit (see Grid.axis).
        distance (float): The star's distance in pc (for axis "au").

    Returns:
        float: The axis factor (arcsec -> `axis` units).
    """
    factor, name = grid.axis(axis, distance)
    ax.set_facecolor("black")
    ax.set(xlim=[e*factor for e in grid.extent[:2]], ylim=[e*factor for e in grid.extent[2:]], aspect="equal",
           xlabel="x [%s]" % name, ylabel="y [%s]" % name)
    return factor


def radial_profile(image, grid):
    """
    Azimuthal mean in 1-pixel-wide rings around the star.

    Each pixel goes to the ring of its distance from the star rounded to whole pixels; empty
    rings are left out.

    Args:
        image (np.ndarray): (N, N) image on `grid`.
        grid (Grid): The image's grid.

    Returns:
        tuple: (radius, mean): 1-D float64 arrays of the same length; radius is the mean distance
            of each ring's pixels from the star in arcsec (0 for the star pixel of an odd grid),
            mean the ring's mean value (in the image's unit).
    """
    r = np.hypot(*grid.coordinates())/grid.pixel_scale
    rings = np.floor(r + 0.5).astype(int)  # ring 0 is the star pixel(s)
    counts = np.bincount(rings.ravel())
    sums = np.bincount(rings.ravel(), weights=np.asarray(image, dtype=float).ravel())
    keep = counts > 0
    radius = np.bincount(rings.ravel(), weights=r.ravel())[keep]/counts[keep]
    return radius*grid.pixel_scale, sums[keep]/counts[keep]


def plot_profiles(ax, curves, grid, axis="arcsec", distance=None, log=True, label="", extra=None, styles=None,
                  log_range=8, legend=True, **kwargs):
    """
    Plot radial profiles (radial_profile) of images against the distance from the star.

    Curves named "total" or "input" are drawn thicker. A curve's style: these defaults, then its entry in
    `styles`, then kwargs (every curve). Log:
    log-log axes, at most `log_range` decades below the highest point (points at radius 0 or
    with a non-positive mean are left out; linear axes if no point is left, e.g. a run without
    components); linear: both axes start at 0.

    Args:
        ax (matplotlib.axes.Axes): Where to draw.
        curves (dict): {label: (N, N) image on `grid`}.
        grid (Grid): The images' grid.
        axis (str): Distance unit (see Grid.axis).
        distance (float): The star's distance in pc (for axis "au").
        log (bool): Log-log axes.
        label (str): The values' unit, for the y label ("" for none).
        extra (dict): More curves {label: (radius in arcsec, values)}, drawn dashed.
        styles (dict): {curve label: matplotlib line keywords (e.g. color, ls, lw)}; None: the defaults.
        log_range (float): Log axes: decades shown below the highest point.
        legend (bool): Draw the legend.
        **kwargs: Passed to matplotlib's plot for every curve, over `styles`, e.g. lw=2.
    """
    factor, name = grid.axis(axis, distance)
    top = []  # the highest point of each curve (log axes need some positive data: not so for a run without components)
    for curve, image in curves.items():
        r, mean = radial_profile(image, grid)
        keep = (r > 0) & (mean > 0) if log else slice(None)
        style = {"lw": 2 if curve in ("total", "input") else 1.2, **(styles or {}).get(curve, {}), **kwargs}
        ax.plot(r[keep]*factor, mean[keep], label=curve, **style)
        top += [mean[keep].max()] if np.any(keep) else []
    for curve, (r, values) in (extra or {}).items():
        ax.plot(np.asarray(r)*factor, values, **{"ls": "--", **kwargs}, label=curve)
        top += [np.max(values)]
    if log and top:
        ax.set(xscale="log", yscale="log", ylim=(max(top)*10**-log_range, max(top)*2))
    else:
        ax.set(xlim=(0, None), ylim=(0, None))
    ax.set(xlabel="distance from the star [%s]" % name, ylabel=("azimuthal mean [%s]" % label) if label else "azimuthal mean")
    if legend:
        ax.legend(fontsize=8)


def image_hdu(name, image, unit, grid):
    """
    A FITS image extension of a scene or final image.

    Header: EXTNAME (name in capitals), BUNIT, PIXSCALE (arcsec), WAVELEN (m), and a linear
    WCS of arcsec offsets from the star (CRPIX = (N+1)/2: pixel (N-1)/2 counted from 0).

    Args:
        name (str): The extension name (upper-cased).
        image (np.ndarray): (N, N) image, written as float64.
        unit (str): Its unit (BUNIT).
        grid (Grid): Its grid.

    Returns:
        astropy.io.fits.ImageHDU: The extension.
    """
    hdu = fits.ImageHDU(np.asarray(image, dtype=np.float64), name=name.upper())
    e = hdu.header
    e["BUNIT"] = unit
    e["PIXSCALE"] = (grid.pixel_scale, "arcsec per pixel")
    e["WAVELEN"] = (grid.wavelength, "m")
    e["CTYPE1"], e["CTYPE2"] = "OFFSET-X", "OFFSET-Y"
    e["CUNIT1"] = e["CUNIT2"] = "arcsec"
    e["CRPIX1"] = e["CRPIX2"] = ((grid.size+1)/2, "the star (FITS pixels start at 1)")
    e["CRVAL1"] = e["CRVAL2"] = 0.0
    e["CDELT1"] = e["CDELT2"] = grid.pixel_scale
    return hdu


def write_scene(path, images, unit, grid, star=None, overwrite=False, extra=None):
    """
    Write a multi-extension scene FITS file.

    The primary HDU holds the grid (PIXSCALE arcsec, WAVELEN m, DIAMETER m), the star's known
    fields (DISTANCE pc, RADIUS R_sun, TEMPERAT K, FLUX_JY Jy) and the `extra` header cards;
    then one image extension per component (see image_hdu).

    Args:
        path (str): The FITS file to write.
        images (dict): {component name: (N, N) image in `unit`}.
        unit (str): The images' unit (BUNIT).
        grid (Grid): Their grid.
        star (Star): The star; None: no star cards.
        overwrite (bool): Replace an existing file (else astropy raises OSError).
        extra (dict): More primary header cards {key: value or (value, comment)}.
    """
    primary = fits.PrimaryHDU()
    h = primary.header
    h["ORIGIN"] = "toycoronagraph"
    h["SCNSIZE"] = (grid.size, "scene size [pixels]")
    h["PIXSCALE"] = (grid.pixel_scale, "arcsec per pixel")
    h["WAVELEN"] = (grid.wavelength, "m")
    h["DIAMETER"] = (grid.diameter, "telescope diameter [m]")
    for key, comment in [("distance", "pc"), ("radius", "R_sun"), ("temperature", "K"), ("flux_jy", "Jy")]:
        if getattr(star, key, None) is not None:
            h[key.upper()[:8]] = (getattr(star, key), comment)  # FITS keywords have at most 8 characters
    for key, value in (extra or {}).items():
        h[key] = value
    hdus = [primary] + [image_hdu(name, img, unit, grid) for name, img in images.items()]
    fits.HDUList(hdus).writeto(path, overwrite=overwrite)
