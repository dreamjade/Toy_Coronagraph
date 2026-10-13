"""
Units of the scene components, and the conversion to one internal unit.

Two families:
    physical   a brightness: any unit astropy converts to Jy (per pixel) or Jy/sr (per solid
               angle) with the spectral_density equivalency, e.g. Jy/arcsec^2, MJy/sr, W/m^2
               per pixel (nu F_nu), erg/s/cm^2/nm per pixel; and photons/s/cm^2/pixel over a band
    relative   a fraction of the star's total flux: "contrast_pixel" (per pixel) or
               "contrast/arcsec^2"
The link between them is the star's flux density F_star (Jy) at the wavelength.

Every component is first put per pixel of its own grid (so that resampling conserves the
flux): Jy per pixel, or the fraction of F_star per pixel. The scene's internal unit is then
Jy/arcsec^2 when F_star is known (or every component is physical), else contrast_pixel.
"""
import numbers
from dataclasses import dataclass
import numpy as np
from astropy import units as u, constants as const

# Names that can be used in parameter files and FITS headers, as astropy units (per pixel unless per area)
NAMED = {
    "Jy/arcsec^2": u.Jy/u.arcsec**2, "Jy/pixel": u.Jy, "MJy/sr": u.MJy/u.sr, "W/m^2/pixel": u.W/u.m**2,
    "contrast_pixel": u.dimensionless_unscaled, "contrast": u.dimensionless_unscaled, "contrast/arcsec^2": u.arcsec**-2,
}
PHOTONS = "photons/s/cm^2/pixel"  # needs a band to become Jy


class UnitError(ValueError):
    """A unit that is unknown, not convertible, or missing information (band, star's flux) to be converted."""


@dataclass
class SourceImage:
    """
    What a source module's render() returns: an image, its unit and how to place it on the scene grid.

    The image must be square and centred on the star (pixel (N-1)/2, rows = +y as in FITS).

    Attributes:
        data (np.ndarray or astropy.units.Quantity): (N, N) image; converted to float64 on creation.
            A Quantity's unit replaces `unit`.
        unit (str or astropy.units.UnitBase): a NAMED unit (e.g. "contrast_pixel", "Jy/arcsec^2"),
            PHOTONS, or anything astropy reads as a unit; a bare astropy unit without "/pixel"
            or a solid angle (e.g. u.Jy, dimensionless) is per pixel.
        pixel_scale (float): arcsec per pixel of data; None = already on the scene grid (the
            scene then sets it to the grid's pixel scale).
        wavelength (float): m, the wavelength of a physical image (used by the spectral_density
            conversions); None = the scene's wavelength.
        band (tuple): (lambda_min, lambda_max) in m, needed only for photon counts (PHOTONS).
        origin (str): text naming the image in messages (e.g. the file name).
        notes (tuple of str): messages for the user (shown with the scene's notes).

    Raises:
        UnitError: when the unit is None and data is not a Quantity.
    """
    data: object
    unit: object = None
    pixel_scale: float = None
    wavelength: float = None
    band: tuple = None
    origin: str = ""
    notes: tuple = ()

    def __post_init__(self):
        """Take the unit from a Quantity, make data a float64 array, and check that a unit is set."""
        if isinstance(self.data, u.Quantity):
            self.unit, self.data = self.data.unit, self.data.value
        self.data = np.asarray(self.data, dtype=np.float64)
        if self.unit is None:
            raise UnitError("%s: the image has no unit" % (self.origin or "a source"))


def per_pixel(img, wavelength):
    """
    The image per pixel of its own grid, in the base unit of its family.

    The unit is tried, in order, as a contrast per pixel, a contrast per solid angle (multiplied
    by the pixel area), a brightness per pixel (-> Jy with the spectral_density equivalency at
    the wavelength), and a brightness per solid angle (x pixel area -> Jy). Photon counts
    (PHOTONS) become Jy as F_nu = N h c / lambda / (c/lambda_min - c/lambda_max), i.e. the
    photons' energy at the wavelength spread evenly over the band in frequency.

    Args:
        img (SourceImage): with pixel_scale set (arcsec).
        wavelength (float): m, used when img.wavelength is None.

    Returns:
        tuple: (family, values): family is "physical" (values in Jy per pixel) or "relative"
            (values as the fraction of the star's flux per pixel); values is a float64 array of
            the shape of img.data.

    Raises:
        UnitError: photon counts without a band, an unknown unit, or a unit of neither family.
    """
    lam = (img.wavelength or wavelength)*u.m
    if img.unit == PHOTONS:
        if not img.band:
            raise UnitError("%s: photons/s/cm^2/pixel needs the band (lambda_min, lambda_max) to become Jy; "
                            "for grey (Thomson) scattering, convert to contrast instead" % (img.origin or "a source"))
        lo, hi = sorted(img.band)
        dnu = const.c/(lo*u.m) - const.c/(hi*u.m)
        return "physical", img.data*(const.h*const.c/lam/u.s/u.cm**2/dnu).to_value(u.Jy)
    try:
        unit = NAMED[img.unit] if isinstance(img.unit, str) and img.unit in NAMED else u.Unit(img.unit)
    except ValueError as err:
        raise UnitError("%s: unknown unit %r (%s)" % (img.origin or "a source", img.unit, err))
    area = img.pixel_scale**2*u.arcsec**2
    eq = u.spectral_density(lam)
    for fam, target, factor in [("relative", u.dimensionless_unscaled, 1), ("relative", u.dimensionless_unscaled, area),
                                ("physical", u.Jy, 1), ("physical", u.Jy, area)]:  # per pixel, then per area
        try:
            return fam, img.data*(1*unit*factor).to_value(target, eq)
        except u.UnitConversionError:
            pass
    raise UnitError("%s: unit %s is neither a brightness (convertible to Jy or Jy/sr) nor a contrast"
                    % (img.origin or "a source", unit))


def family(unit):
    """
    The family of a unit, found by converting a dummy image with per_pixel.

    Args:
        unit (str or astropy.units.UnitBase): a NAMED unit, PHOTONS, or an astropy unit.

    Returns:
        str: "physical" or "relative".

    Raises:
        UnitError: for an unknown unit or a unit of neither family.
    """
    return per_pixel(SourceImage(np.ones(1), unit, pixel_scale=1.0, band=(1e-6, 2e-6)), 1e-6)[0]


def choose_internal(families, star_flux):
    """
    The scene's internal unit.

    Args:
        families (dict): {component name: "physical" or "relative"}.
        star_flux (float): the star's flux density in Jy, or None if unknown.

    Returns:
        str: "Jy/arcsec^2" when the star's flux is known or every component is physical;
            "contrast_pixel" when every component is relative and the flux is unknown.

    Raises:
        UnitError: physical and relative components together without the star's flux.
    """
    fam = set(families.values())
    if star_flux is not None or fam <= {"physical"}:
        return "Jy/arcsec^2"
    if fam <= {"relative"}:
        return "contrast_pixel"
    raise UnitError("the components are in different units (%s): the star's flux is needed to combine them; set "
                    "[star] temperature (with radius and distance) or flux_jy"
                    % ", ".join("%s: %s" % (n, "contrast" if f == "relative" else "physical") for n, f in families.items()))


def to_internal(fam, values, internal, pixel_area, star_flux):
    """
    Per-pixel values on the scene grid -> the internal unit.

    Args:
        fam (str): "physical" (values in Jy per pixel) or "relative" (fraction of the star's flux
            per pixel), as returned by per_pixel.
        values (np.ndarray or float): per-pixel values (any shape; also the fluxes of point sources).
        internal (str): "Jy/arcsec^2" or "contrast_pixel" (choose_internal).
        pixel_area (float): arcsec^2 per pixel of the scene grid.
        star_flux (float): the star's flux density in Jy (may be None when no conversion
            between families is needed).

    Returns:
        np.ndarray or float: values in the internal unit, same shape.
    """
    if internal == "Jy/arcsec^2":
        return (values if fam == "physical" else values*star_flux)/pixel_area
    return values/star_flux if fam == "physical" else values


def star_flux_of(star, wavelength):
    """
    A star's flux density in Jy.

    Args:
        star (float or Star or None): a real number (Jy; a NumPy scalar too) or None, returned as it
            is; or a toycoronagraph.scene.Star, evaluated at the wavelength (Star.flux_at).
        wavelength (float): m.

    Returns:
        float: Jy (None when star is None).

    Raises:
        UnitError: a Star with neither flux_jy nor temperature, radius and distance.
    """
    if star is None or isinstance(star, numbers.Real):
        return star
    flux = star.flux_at(wavelength)
    if flux is None:
        raise UnitError("the star's flux needs temperature (with radius and distance) or flux_jy")
    return flux


OUTPUT_UNITS = ("Jy/arcsec^2", "Jy/pixel", "contrast_pixel", "contrast_peak")


def flux_per_pixel(unit, pixel_area):
    """
    What turns an image in the internal unit into each pixel's flux (summing over pixels gives a total):
    the pixel area for "Jy/arcsec^2" (-> Jy), 1 for "contrast_pixel" (already a fraction of the star's flux).

    Args:
        unit (str): The internal unit, "Jy/arcsec^2" or "contrast_pixel".
        pixel_area (float): arcsec^2 per pixel.

    Returns:
        float: the factor.
    """
    return pixel_area if unit == "Jy/arcsec^2" else 1.0


def convert_image(image, from_unit, to_unit, pixel_area, star_flux=None, peak_norm=None):
    """
    A scene image in another output unit.

    Args:
        image (np.ndarray): (N, N) image in from_unit.
        from_unit (str): the scene's internal unit, "Jy/arcsec^2" or "contrast_pixel".
        to_unit (str): "Jy/arcsec^2", "Jy/pixel", "contrast_pixel", or "contrast_peak" (surface
            brightness / the peak of the star's unocculted PSF).
        pixel_area (float): arcsec^2 per pixel.
        star_flux (float): the star's flux density in Jy; needed to go between the families.
        peak_norm (float): for "contrast_peak": the peak of the star's unocculted PSF as a
            fraction of the star's flux per pixel.

    Returns:
        np.ndarray: (N, N) image in to_unit (the input itself when the units are equal).

    Raises:
        UnitError: an unknown unit, a from_unit that is not an internal unit, or the star's flux or
            peak_norm is needed but None.
    """
    for unit in (from_unit, to_unit):
        if unit not in OUTPUT_UNITS:
            raise UnitError("unknown output unit %r; use one of %s" % (unit, ", ".join(OUTPUT_UNITS)))
    if from_unit not in ("Jy/arcsec^2", "contrast_pixel"):
        raise UnitError("cannot convert from %s: the image must be in an internal unit, Jy/arcsec^2 or contrast_pixel"
                        % from_unit)
    if from_unit == to_unit:
        return image
    physical_in, physical_out = from_unit == "Jy/arcsec^2", to_unit in ("Jy/pixel", "Jy/arcsec^2")
    if physical_in != physical_out and star_flux is None:
        raise UnitError("%s needs the star's flux ([star] temperature or flux_jy)" % to_unit)
    if physical_in:                # Jy per pixel
        jy = image*pixel_area
    elif physical_out:
        jy = image*star_flux
    else:                          # contrast to contrast: not needed
        jy = None
    if physical_out:
        return jy if to_unit == "Jy/pixel" else jy/pixel_area
    frac = jy/star_flux if physical_in else image  # fraction of the star's flux per pixel
    if to_unit == "contrast_pixel":
        return frac
    if peak_norm is None:
        raise UnitError("contrast_peak needs the peak of the star's unocculted PSF (per coronagraph)")
    return frac/peak_norm


# Spellings of BUNIT recognised for each unit (compared in lower case, without spaces, "**" read as "^")
_HEADER_SPELLINGS = {
    "W/m^2/pixel": ("w/m^2/pixel", "w.m-2.pixel-1", "w/m2/pixel", "w/m^2/pix", "w.m-2.pix-1"),
    "Jy/arcsec^2": ("jy/arcsec^2", "jy.arcsec-2", "jy/arcsec2", "jy/sq.arcsec"),
    "Jy/pixel": ("jy/pixel", "jy.pixel-1", "jy/pix", "jy.pix-1"),
    "MJy/sr": ("mjy/sr", "mjy.sr-1"),
    "contrast_pixel": ("contrast", "contrast_pixel", "contrast/pixel", "contrast.pixel-1"),
    "contrast/arcsec^2": ("contrast/arcsec^2", "contrast/arcsec2", "contrast.arcsec-2"),
    PHOTONS: ("photons/s/cm^2/pixel", "photons/s/cm2/pixel", "ph/s/cm2/pix"),
}


def header_unit(header):
    """
    The unit named by a FITS header's BUNIT.

    Args:
        header (astropy.io.fits.Header or dict-like): the image's header.

    Returns:
        tuple: (unit, raw): unit is a NAMED unit, PHOTONS, the raw text when astropy reads it as
            a unit of either family, or None (no BUNIT, or not understood); raw is the BUNIT text
            ("" when absent).
    """
    raw = str(header.get("BUNIT", "")).strip()
    key = raw.lower().replace(" ", "").replace("**", "^")
    for unit, spellings in _HEADER_SPELLINGS.items():
        if key in spellings:
            return unit, raw
    if not raw:  # (astropy would read "" as dimensionless, i.e. a contrast)
        return None, raw
    try:
        family(raw)
        return raw, raw
    except (UnitError, ValueError):
        return None, raw
