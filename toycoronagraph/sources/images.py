"""Image source modules: a FITS image, an MCFOST model, and the scene FITS interchange file."""
import os
from astropy.io import fits
from toycoronagraph import DATADIR
from toycoronagraph.params import Param
from toycoronagraph.sources import Source, register_source
from toycoronagraph.units import SourceImage, header_unit

INPUT_UNITS = ("auto", "W/m^2/pixel", "Jy/arcsec^2", "Jy/pixel", "MJy/sr", "contrast_pixel")
MCFOST_PLANES = ("disk", "total", "I", "Q", "U", "V", "direct_star", "scattered_star", "direct_thermal", "scattered_thermal")


def fits_path(name):
    """
    The path of a FITS file named in a parameter file.

    Args:
        name (str): a path; "example" = the bundled example (DATADIR/example.fits).

    Returns:
        str: the path, with ".fits" added when it has no FITS extension (.fits, .fits.gz, .fit, .fts).
    """
    if name == "example":
        return os.path.join(DATADIR, "example.fits")
    if name.lower().endswith((".fits", ".fits.gz", ".fit", ".fts")):
        return name
    return name + ".fits"


def _header_scale_and_wavelength(header):
    """
    The pixel scale and wavelength written in a FITS header (see config._header_pixel_scale / _header_wavelength).

    Args:
        header (astropy.io.fits.Header): the image's header.

    Returns:
        tuple: (scale, scale_key, lam, lam_key): arcsec per pixel and m, each with the keyword it
            came from (e.g. "PIXSCALE", "CDELT1 (deg)", "WAVE (um)"); (None, None) for a missing one.
    """
    from toycoronagraph.config import _header_pixel_scale, _header_wavelength
    scale, scale_key = _header_pixel_scale(header)
    lam, lam_key = _header_wavelength(header)
    return scale, scale_key, lam, lam_key


class _ImageSource(Source):
    """
    Shared reading of one 2-D image from a FITS file (base of FitsImage and Mcfost).

    Subclasses define _read_args(); the optional settings input_unit, input_pixel_scale and
    input_wavelength ("auto" = from the header) override the header. The file is read once
    and cached in self._cache.
    """
    provides = ("grid.size", "grid.pixel_scale", "observation.wavelength")

    def _read_args(self):
        """
        Keyword arguments of inputs.read_plane selecting the image (defined by each subclass).

        Returns:
            dict: e.g. dict(input_format="image", fits_hdu=0, fits_slice="auto").
        """
        raise NotImplementedError

    def _load(self):
        """
        Read the image and its unit, pixel scale and wavelength (once; cached).

        Returns:
            dict: image ((N, N) float64), unit (a NAMED unit, PHOTONS or an astropy unit text),
                unit_from ("input_unit" or "header BUNIT"), scale (arcsec/pixel or None),
                scale_key, lam (m or None), lam_key (where each came from: a header keyword or
                "input_..."), name (file name), description (file name and selected plane),
                index (the indices of the plane in the file).

        Raises:
            FileNotFoundError: the file does not exist.
            ValueError: the image is not square, or its unit is unknown and input_unit is "auto".
        """
        if getattr(self, "_cache", None) is None:
            from toycoronagraph.inputs import read_plane
            path = fits_path(self.settings["file"])
            if not os.path.exists(path):
                raise FileNotFoundError("%s: file %s not found" % (self.kind, path))
            image, header, fmt, what, index = read_plane(path, **self._read_args())
            if image.shape[0] != image.shape[1]:
                raise ValueError("%s: the image must be square (got %s); the star must be at its centre" % (path, image.shape))
            unit, raw = header_unit(header)
            # the settings below exist only for FitsImage (Mcfost has none: .get with "auto")
            if self.settings.get("input_unit", "auto") != "auto":
                unit, unit_from = self.settings["input_unit"], "input_unit"
            elif unit is None:
                raise ValueError("%s: the unit is unknown (BUNIT = %r); set input_unit in its [[source]] table" % (path, raw))
            else:
                unit_from = "header BUNIT"
            scale, scale_key, lam, lam_key = _header_scale_and_wavelength(header)
            if self.settings.get("input_pixel_scale", "auto") != "auto":
                scale, scale_key = self.settings["input_pixel_scale"], "input_pixel_scale"
            if self.settings.get("input_wavelength", "auto") != "auto":
                lam, lam_key = self.settings["input_wavelength"], "input_wavelength"
            name = os.path.basename(path)
            self._cache = dict(image=image, unit=unit, unit_from=unit_from, scale=scale, scale_key=scale_key, lam=lam,
                               lam_key=lam_key, name=name, description="%s (%s)" % (name, what), index=index)
        return self._cache

    def inspect(self):
        """
        The grid size, and the pixel scale and wavelength when known, from the file.

        Returns:
            list of tuple: [(key, value, origin)] for grid.size, grid.pixel_scale (arcsec) and
                observation.wavelength (m).
        """
        c = self._load()
        out = [("grid.size", c["image"].shape[0], "image size of " + c["name"])]
        if c["scale"]:
            out.append(("grid.pixel_scale", c["scale"], "%s of %s" % (_origin(c["scale_key"]), c["name"])))
        if c["lam"]:
            out.append(("observation.wavelength", c["lam"], "%s of %s" % (_origin(c["lam_key"]), c["name"])))
        return out

    def family_of(self):
        """
        The family given by input_unit; with "auto" (or no input_unit), the class's `family`.

        Returns:
            str: "physical", "relative", or None (decided by the header's BUNIT).
        """
        unit = self.settings.get("input_unit", "auto")
        if unit == "auto":
            return self.family
        return "relative" if unit.startswith("contrast") else "physical"

    def render(self, grid, star):
        """
        The image as read, on its own grid.

        Args:
            grid (toycoronagraph.scene.Grid): the scene grid (not used: the scene resamples).
            star (toycoronagraph.scene.Star): not used.

        Returns:
            SourceImage: (N, N) image in its unit, with the file's pixel scale (arcsec; None = the
                scene's) and wavelength (m; None = the scene's).
        """
        c = self._load()
        return SourceImage(c["image"], c["unit"], pixel_scale=c["scale"], wavelength=c["lam"], origin=c["description"])


def _origin(key):
    """The origin text of a value: a setting ("input_...") as it is, else "header <keyword>"."""
    return key if key.startswith("input_") else "header " + key


@register_source
class FitsImage(_ImageSource):
    """
    A 2-D image (or a plane of a cube) from a FITS file, the star at pixel (N-1)/2.

    family is None: the unit comes from the header's BUNIT or the input_unit setting.
    """
    kind = "fits_image"
    doc = "a 2-D FITS image, star at the centre"
    params = (
        Param("file", "example", doc='FITS file; "example" = the bundled example disk', type="str",
              help="Path of the FITS file (\".fits\" is added if there is no extension). The image must be "
                   "square with the star at pixel (N-1)/2. \"example\" uses toycoronagraph/example_data/example.fits."),
        Param("input_unit", "auto", doc="unit of the data; \"auto\" = header BUNIT", choices=INPUT_UNITS,
              help="\"contrast_pixel\" means the fraction of the star's total flux in each pixel."),
        Param("input_pixel_scale", "auto", unit="arcsec", doc="the file's pixel scale; \"auto\" = header PIXSCALE/CDELT",
              auto=True, min=0, exclusive_min=True, advanced=True,
              help="Only needed when the header has no PIXSCALE or CDELT. This describes the file; the "
                   "simulated grid is set in [grid] (the image is resampled if they differ)."),
        Param("input_wavelength", "auto", unit="m", doc="the file's wavelength; \"auto\" = header WAVELEN/WAVE",
              auto=True, min=0, exclusive_min=True, advanced=True),
        Param("hdu", 0, doc="FITS extension holding the image", type="int", min=0, advanced=True),
        Param("slice", "auto", doc='indices picking a 2-D plane of a cube, e.g. [0]', type="int_list", auto=True,
              advanced=True),
        Param("symmetric", "auto", doc="circularly symmetric? \"auto\" = measure it", choices=("auto", True, False),
              advanced=True),
    )

    def _read_args(self):
        """read_plane arguments: a plain image, from extension `hdu`, plane `slice`."""
        return dict(input_format="image", fits_hdu=self.settings["hdu"], fits_slice=self.settings["slice"])


@register_source
class Mcfost(_ImageSource):
    """One plane of an MCFOST image cube (axes: plane, azimuth, inclination, y, x); its unit is the header's BUNIT."""
    kind = "mcfost"
    doc = "an MCFOST image cube (RT.fits.gz)"
    family = "physical"
    params = (
        Param("file", None, doc="MCFOST image cube, e.g. data_0.6/RT.fits.gz", type="str", example="data_0.6/RT.fits.gz"),
        Param("component", "disk", doc='plane: "disk" = everything but the star', choices=MCFOST_PLANES,
              help="\"disk\" is the total minus the direct starlight. \"total\" includes the star itself, which "
                   "the coronagraph is meant to remove, so it is rarely what you want."),
        Param("inclination", 0, doc="index along the inclination axis", type="int", min=0),
        Param("azimuth", 0, doc="index along the azimuth axis", type="int", min=0, advanced=True),
        Param("symmetric", "auto", doc="circularly symmetric? \"auto\" = measure it", choices=("auto", True, False),
              advanced=True),
    )

    def _read_args(self):
        """read_plane arguments: an MCFOST cube, with the chosen component, inclination and azimuth indices."""
        return dict(input_format="mcfost", mcfost_component=self.settings["component"],
                    mcfost_inclination=self.settings["inclination"], mcfost_azimuth=self.settings["azimuth"])


@register_source
class SceneFits(Source):
    """
    The images of a scene written by write_scene(): one component per 2-D image extension.

    Each extension's unit, pixel scale and wavelength come from its own header (family None);
    the TOTAL and FINAL_TOTAL extensions are skipped. The star's distance can come from the
    primary header (DISTANCE).
    """
    kind = "scene_fits"
    doc = "a multi-extension FITS from write_scene(), one image per component"
    provides = ("grid.size", "grid.pixel_scale", "observation.wavelength", "star.distance")
    params = (
        Param("file", None, doc="the scene FITS file", type="str", example="my_scene.fits"),
        Param("extensions", "all", doc='"all", or a list of EXTNAMEs, e.g. ["CME"]', type="str_list", auto_word="all"),
    )

    def _load(self):
        """
        Read the selected image extensions (once; cached).

        Returns:
            dict: images ({name: (N, N) float64}, name = the EXTNAME in lower case, "image" if
                none), headers ({name: header}), primary (the primary header), name (file name).

        Raises:
            ValueError: no image extension is selected.
        """
        if getattr(self, "_cache", None) is None:
            path = fits_path(self.settings["file"])
            want = self.settings["extensions"]
            images, headers = {}, {}
            with fits.open(path) as hdul:
                for hdu in hdul:
                    name = hdu.header.get("EXTNAME", "")
                    if hdu.data is None or hdu.data.ndim != 2 or name.upper() in ("TOTAL", "FINAL_TOTAL"):
                        continue
                    if want != "all" and name not in want:
                        continue
                    key = name.lower() or "image"
                    images[key] = hdu.data.astype(float)
                    headers[key] = hdu.header
                primary = hdul[0].header
            if not images:
                raise ValueError("%s: no image extensions (selected: %s)" % (path, self.settings["extensions"]))
            self._cache = dict(images=images, headers=headers, primary=primary, name=os.path.basename(path))
        return self._cache

    def inspect(self):
        """
        The grid size, pixel scale and wavelength (from the first image's header) and the star's distance.

        Returns:
            list of tuple: [(key, value, origin)] for grid.size, grid.pixel_scale (arcsec),
                observation.wavelength (m) and star.distance (pc), those present.
        """
        c = self._load()
        first = next(iter(c["headers"].values()))
        out = [("grid.size", first["NAXIS1"], "image size of " + c["name"])]
        scale, scale_key, lam, lam_key = _header_scale_and_wavelength(first)
        if scale:
            out.append(("grid.pixel_scale", scale, "header %s of %s" % (scale_key, c["name"])))
        if lam:
            out.append(("observation.wavelength", lam, "header %s of %s" % (lam_key, c["name"])))
        if "DISTANCE" in c["primary"]:
            out.append(("star.distance", float(c["primary"]["DISTANCE"]), "header DISTANCE of " + c["name"]))
        return out

    def family_of(self):
        """Always None: each extension's BUNIT decides."""
        return None

    def render(self, grid, star):
        """
        One image per extension.

        Args:
            grid (toycoronagraph.scene.Grid): not used (the scene resamples).
            star (toycoronagraph.scene.Star): not used.

        Returns:
            dict: {extension name (lower case): SourceImage} with each header's unit, pixel scale
                (arcsec) and wavelength (m).

        Raises:
            ValueError: an extension's BUNIT is not understood.
        """
        c = self._load()
        out = {}
        for name, data in c["images"].items():
            h = c["headers"][name]
            unit, raw = header_unit(h)
            if unit is None:
                raise ValueError("%s [%s]: unknown BUNIT %r" % (c["name"], name, raw))
            scale, _, lam, _ = _header_scale_and_wavelength(h)
            out[name] = SourceImage(data, unit, pixel_scale=scale, wavelength=lam, origin="%s [%s]" % (c["name"], name.upper()))
        return out
