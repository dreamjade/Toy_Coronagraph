"""
Unit tests for stage F1 (tests/scene_stage1_core.py): the scene framework core.

Modules covered:
    toycoronagraph/scene.py             Grid, Star, resample, Scene (render, names, units, point sources,
                                        [grid] resample = false, in_unit, write / scene_fits, plot),
                                        radial_profile
    toycoronagraph/units.py             SourceImage, per_pixel, family, choose_internal, to_internal,
                                        convert_image, flux_per_pixel, star_flux_of, header_unit
    toycoronagraph/sources/__init__.py  the registry (register_source, get_source, source_kinds) and Source
    toycoronagraph/sources/images.py    fits_path, FitsImage (header grid, units, errors)
    toycoronagraph/inputs.py            read_plane (planes of cubes, the "N-D after selecting" message)

Images are 16-48 pixels on a side; even and odd sizes are both tested where the centre
convention matters (star at pixel (N-1)/2, rows = +y). Files are written in each test's own
temporary folder (in_tmp). The stage script keeps the figures and the 95/96-pixel runs.
"""
import re

import numpy as np
import pytest
from astropy import units as u, constants as const
from astropy.io import fits

from extra_sources import CmeBlob, QuantityRing, Ring  # (importing it registers the test modules)
from toycoronagraph import params as P, sources as S, units as U
from toycoronagraph.inputs import read_plane
from toycoronagraph.scene import Grid, Scene, Star, radial_profile, resample
from toycoronagraph.sources.images import FitsImage, SceneFits, fits_path
from toycoronagraph.sources.planets import Planet

SIZES = [24, 25]  # an even and an odd size


def _ring(n, radius, width):
    """A Gaussian ring centred on pixel ((n-1)/2, (n-1)/2); radius and width in pixels; (n, n) float64."""
    c = (n-1)/2
    rows, cols = np.indices((n, n))
    return np.exp(-0.5*((np.hypot(rows-c, cols-c)-radius)/width)**2)


def _write_fits(path, data, **cards):
    """Write data as a primary-HDU FITS file with the given header cards; returns the path (str)."""
    fits.PrimaryHDU(np.asarray(data), fits.Header(list(cards.items()))).writeto(str(path))
    return str(path)


def _disk(folder, name, n, scale, wavelength=1.6e-6):
    """
    A ring image file in Jy/arcsec^2 with PIXSCALE and WAVELEN in its header.

    Args:
        folder (pathlib.Path): where to write it.
        name (str): file name without ".fits".
        n (int): pixels on a side.
        scale (float): PIXSCALE, arcsec per pixel.
        wavelength (float): WAVELEN, m.

    Returns:
        tuple: (path (str), image ((n, n) float64)).
    """
    image = _ring(n, 0.25*n, 0.04*n)
    path = _write_fits(folder/(name+".fits"), image, BUNIT="Jy/arcsec^2", PIXSCALE=scale, WAVELEN=wavelength)
    return path, image


def _centroid_offset(image):
    """The (x, y) = (column, row) offset of an image's centroid from pixel ((N-1)/2, (N-1)/2)."""
    n = image.shape[0]
    rows, cols = np.indices(image.shape)
    return (image*cols).sum()/image.sum() - (n-1)/2, (image*rows).sum()/image.sum() - (n-1)/2


def _render_error(scene):
    """The ParamError message of scene.render(), or None if it renders."""
    try:
        scene.render(verbose=False)
    except P.ParamError as err:
        return str(err)
    return None


# ====================================================================== Grid

def test_grid_derived_quantities():
    """lambda_d (arcsec), sampling, visual_range, pixel_area and extent follow from size, scale, lambda, D."""
    g = Grid(size=25, pixel_scale=0.02, wavelength=1.6e-6, diameter=2.4)
    lam_d = 1.6e-6/2.4*P.RAD_TO_ARCSEC
    assert np.isclose(g.lambda_d, lam_d) and np.isclose(g.sampling, lam_d/0.02)
    assert np.isclose(g.visual_range, 25*0.02/lam_d) and np.isclose(g.pixel_area, 4e-4)
    assert np.allclose(g.extent, [-0.25, 0.25, -0.25, 0.25])


@pytest.mark.parametrize("n", SIZES)
def test_grid_coordinates_star_at_centre_rows_are_y(n):
    """coordinates(): (n, n) float64 offsets, zero at pixel (n-1)/2; x grows along columns, y along rows."""
    g = Grid(n, 0.1, 1e-6, 2.0)
    x, y = g.coordinates()
    assert x.shape == y.shape == (n, n) and x.dtype == y.dtype == np.float64
    assert np.allclose(np.diff(x, axis=1), 0.1) and np.allclose(np.diff(x, axis=0), 0)
    assert np.allclose(np.diff(y, axis=0), 0.1) and np.allclose(np.diff(y, axis=1), 0)
    assert np.isclose(x[0, 0], -(n-1)/2*0.1) and np.allclose(x, -x[:, ::-1]) and np.allclose(y, -y[::-1])
    if n % 2:
        assert x[n//2, n//2] == 0 and y[n//2, n//2] == 0


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("angle, d_row, d_col", [(0.0, 0, 1), (np.pi/2, 1, 0), (np.pi, 0, -1), (-np.pi/2, -1, 0)])
def test_grid_pixel_angles_from_x_towards_y(n, angle, d_row, d_col):
    """pixel(r, angle): angle from +x (columns) towards +y (rows), around (n-1)/2."""
    g = Grid(n, 0.1, 1e-6, 2.0)
    row, col = g.pixel(0.5, angle)
    c = (n-1)/2
    assert np.isclose(row, c + 5*d_row) and np.isclose(col, c + 5*d_col)


def test_grid_axis_factors_and_errors():
    """axis(): factors from arcsec to each unit; "au" needs the distance; unknown units are refused."""
    g = Grid(16, 0.5, 1e-6, 2.0)
    assert g.axis("arcsec") == (1.0, "arcsec") and g.axis("mas") == (1e3, "mas") and g.axis("pixel") == (2.0, "pixel")
    assert g.axis("au", 3.0) == (3.0, "au") and np.isclose(g.axis("lambda/D")[0], 1/g.lambda_d)
    with pytest.raises(P.ParamError, match="needs the distance"):
        g.axis("au")
    with pytest.raises(P.ParamError, match="use one of"):
        g.axis("km")


def test_grid_is_frozen_and_compared_by_value():
    """Grid is immutable and two grids with the same values are equal."""
    g = Grid(16, 0.5, 1e-6, 2.0)
    assert g == Grid(16, 0.5, 1e-6, 2.0) and g != Grid(17, 0.5, 1e-6, 2.0)
    with pytest.raises(AttributeError):
        g.size = 3


@pytest.mark.parametrize("n", SIZES)
def test_radial_profile_of_a_constant(n):
    """radial_profile(): a constant image has that mean in every ring; radii increase (0 first for odd n)."""
    g = Grid(n, 0.1, 1e-6, 2.0)
    r, mean = radial_profile(np.full((n, n), 3.0), g)
    assert np.allclose(mean, 3.0) and np.all(np.diff(r) > 0) and r.shape == mean.shape
    assert (r[0] == 0) == bool(n % 2)


# ====================================================================== Star

def test_star_flux_jy_overrides_the_blackbody():
    """flux_jy, when set, is the star's flux whatever the temperature."""
    assert Star(distance=3.0, radius=1.0, temperature=5800, flux_jy=7.0).flux_at(1e-6) == 7.0


def test_star_blackbody_matches_planck():
    """The blackbody flux is pi B_nu(T) (R/d)^2 in Jy."""
    lam, t, r, d = 1.6e-6, 5800.0, 1.0, 10.0
    nu = const.c.value/lam
    bnu = 2*const.h.value*nu**3/const.c.value**2/np.expm1(const.h.value*nu/(const.k_B.value*t))
    planck = np.pi*bnu*(r*const.R_sun.value/(d*u.pc.to(u.m)))**2*1e26
    assert abs(Star(distance=d, radius=r, temperature=t).flux_at(lam)/planck - 1) < 1e-10


@pytest.mark.parametrize("star", [Star(), Star(distance=3.0, radius=1.0), Star(temperature=5000, radius=1.0)])
def test_star_unknown_flux_is_none(star):
    """Without flux_jy or all of temperature, radius and distance, the flux is None."""
    assert star.flux_at(1e-6) is None


# ====================================================================== resample

@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("scale_out, extra", [(2.0, 1), (0.5, 0), (1.37, 3), (0.73, 2), (1.0, 6)])
def test_resample_conserves_flux_and_keeps_the_star_centred(n, scale_out, extra):
    """Output covering the input: total flux kept to 1e-12, centroid at (n_out-1)/2."""
    img = _ring(n, 0.2*n, 0.04*n)
    n_out = int(np.ceil(n/scale_out)) + extra
    out = resample(img, 1.0, scale_out, n_out)
    assert out.shape == (n_out, n_out) and out.dtype == np.float64
    assert abs(out.sum()/img.sum() - 1) < 1e-12
    assert np.allclose(_centroid_offset(out), 0, atol=1e-9)


@pytest.mark.parametrize("n, factor", [(24, 2), (24, 3), (25, 5), (27, 3)])
def test_resample_integer_factor_is_a_block_sum(n, factor):
    """An integer factor (same parity) sums blocks of factor x factor pixels."""
    a = np.random.default_rng(0).random((n, n))
    m = n//factor
    assert np.allclose(resample(a, 1.0, float(factor), m), a.reshape(m, factor, m, factor).sum(axis=(1, 3)))


@pytest.mark.parametrize("n", SIZES)
def test_resample_same_grid_is_identity(n):
    """The same pixel scale and size gives the image back."""
    a = np.random.default_rng(1).random((n, n))
    assert np.allclose(resample(a, 0.3, 0.3, n), a, rtol=1e-12, atol=1e-15)


@pytest.mark.parametrize("n", SIZES)
def test_resample_crop_drops_the_light_outside(n):
    """A smaller grid at the same scale is the central crop: the light outside it is dropped."""
    a = np.random.default_rng(2).random((n, n))
    assert np.allclose(resample(a, 1.0, 1.0, n-4), a[2:-2, 2:-2])


def test_resample_point_on_the_star_stays_on_the_star():
    """An odd image's star pixel lands on the star pixel of an odd output grid, flux kept."""
    a = np.zeros((25, 25))
    a[12, 12] = 1.0
    out = resample(a, 1.0, 1.0/3, 75)
    assert np.isclose(out.sum(), 1.0) and np.isclose(out[37, 37], 1/9) and np.isclose(out[36:39, 36:39].sum(), 1.0)


# ====================================================================== units

def test_source_image_takes_the_unit_of_a_quantity_and_makes_float64():
    """SourceImage: a Quantity's unit replaces `unit`; data becomes float64."""
    img = U.SourceImage(np.ones((3, 3), dtype=np.float32)*u.MJy/u.sr)
    assert img.unit == u.MJy/u.sr and img.data.dtype == np.float64 and not isinstance(img.data, u.Quantity)


def test_source_image_without_unit_is_refused():
    """A SourceImage without a unit raises UnitError naming its origin."""
    with pytest.raises(U.UnitError, match="disk.fits: the image has no unit"):
        U.SourceImage(np.ones((3, 3)), origin="disk.fits")


LAM, SCALE, F_STAR = 1.6e-6, 0.01, 3.7  # m, arcsec/pixel, Jy
BASE = _ring(16, 4, 1.0)*2e-3           # Jy/arcsec^2
JY_PIX = BASE*SCALE**2
BAND = (1.0e-6, 2.2e-6)
DNU = (const.c/(BAND[0]*u.m) - const.c/(BAND[1]*u.m)).to_value(u.Hz)
E_PH = (const.h*const.c/(LAM*u.m)).to_value(u.erg)


@pytest.mark.parametrize("unit, data, band", [
    ("Jy/arcsec^2", BASE, None),
    ("Jy/pixel", JY_PIX, None),
    ("MJy/sr", (BASE*u.Jy/u.arcsec**2).to_value(u.MJy/u.sr), None),
    ("W/m^2/pixel", (JY_PIX*u.Jy).to_value(u.W/u.m**2, u.spectral_density(LAM*u.m)), None),
    ("photons/s/cm^2/pixel", JY_PIX*1e-23*DNU/E_PH, BAND),
    ("contrast_pixel", JY_PIX/F_STAR, None),
    ("contrast/arcsec^2", BASE/F_STAR, None),
    ("erg / (s cm2 Hz sr)", (BASE*u.Jy/u.arcsec**2).to_value(u.erg/u.s/u.cm**2/u.Hz/u.sr), None),
    ("erg / (s cm2 nm)", (JY_PIX*u.Jy).to_value(u.erg/u.s/u.cm**2/u.nm, u.spectral_density(LAM*u.m)), None),
])
def test_every_unit_gives_the_same_internal_image(unit, data, band):
    """The same image in any unit becomes the same Jy/arcsec^2 image (via per_pixel and to_internal)."""
    fam, per = U.per_pixel(U.SourceImage(data, unit, pixel_scale=SCALE, wavelength=LAM, band=band), LAM)
    assert fam == ("relative" if unit.startswith("contrast") else "physical")
    internal = U.to_internal(fam, per, "Jy/arcsec^2", SCALE**2, F_STAR)
    assert np.abs(internal - BASE).max()/BASE.max() < 1e-12


@pytest.mark.parametrize("unit, words", [
    ("not_a_unit_xyz", "unknown unit 'not_a_unit_xyz'"),
    ("kg", "neither a brightness"),
])
def test_per_pixel_refuses_unknown_and_wrong_units(unit, words):
    """An unknown unit, or one of neither family, raises UnitError."""
    with pytest.raises(U.UnitError, match=words):
        U.per_pixel(U.SourceImage(np.ones((3, 3)), unit, pixel_scale=0.1), 1e-6)


def test_photons_need_a_band():
    """Photon counts without a band raise UnitError (no way to Jy)."""
    with pytest.raises(U.UnitError, match="needs the band"):
        U.per_pixel(U.SourceImage(np.ones((3, 3)), U.PHOTONS, pixel_scale=0.1), 1e-6)


@pytest.mark.parametrize("unit, family", [("contrast", "relative"), ("contrast/arcsec^2", "relative"), ("Jy", "physical"),
                                          ("MJy/sr", "physical"), (U.PHOTONS, "physical"), (u.dimensionless_unscaled, "relative")])
def test_family_of_units(unit, family):
    """family(): contrasts are "relative", brightnesses (any astropy form) "physical"."""
    assert U.family(unit) == family


@pytest.mark.parametrize("families, flux, expected", [
    ({"a": "relative"}, None, "contrast_pixel"),
    ({"a": "relative", "b": "physical"}, 2.0, "Jy/arcsec^2"),
    ({"a": "physical"}, None, "Jy/arcsec^2"),
    ({"a": "relative"}, 2.0, "Jy/arcsec^2"),
    ({}, None, "Jy/arcsec^2"),
])
def test_choose_internal(families, flux, expected):
    """Jy/arcsec^2 when the star's flux is known or all is physical, else contrast_pixel."""
    assert U.choose_internal(families, flux) == expected


def test_choose_internal_mixed_without_flux():
    """Physical and relative components without the star's flux: UnitError naming each."""
    with pytest.raises(U.UnitError, match="wind: contrast, disk: physical"):
        U.choose_internal({"wind": "relative", "disk": "physical"}, None)


def test_to_internal():
    """to_internal(): per-pixel values -> Jy/arcsec^2 (/ area, x flux for contrasts) or contrast_pixel (/ flux)."""
    v = np.array([1.0, 2.0])
    assert np.allclose(U.to_internal("physical", v, "Jy/arcsec^2", 0.25, None), v/0.25)
    assert np.allclose(U.to_internal("relative", v, "Jy/arcsec^2", 0.25, 2.0), v*2.0/0.25)
    assert np.allclose(U.to_internal("physical", v, "contrast_pixel", 0.25, 2.0), v/2.0)
    assert U.to_internal("relative", v, "contrast_pixel", 0.25, None) is v


@pytest.mark.parametrize("from_unit, to_unit", [
    ("Jy/arcsec^2", "Jy/sr"), ("MJy/sr", "Jy/pixel"), ("contrast", "Jy/pixel"), ("Jy/arcsec^2", "auto"),
])
def test_convert_image_unknown_units_raise(from_unit, to_unit):
    """Regression: an unknown input or output unit raises UnitError listing the known ones."""
    with pytest.raises(U.UnitError, match="unknown output unit"):
        U.convert_image(np.ones((3, 3)), from_unit, to_unit, 0.25, star_flux=2.0)


def test_convert_image_same_unit_returns_the_input():
    """No conversion needed: the same array is returned."""
    img = np.ones((3, 3))
    assert U.convert_image(img, "contrast_pixel", "contrast_pixel", 0.25) is img


@pytest.mark.parametrize("from_unit, to_unit, factor", [
    ("Jy/arcsec^2", "Jy/pixel", 0.25),
    ("Jy/arcsec^2", "contrast_pixel", 0.25/2.0),
    ("Jy/arcsec^2", "contrast_peak", 0.25/2.0/0.5),
    ("contrast_pixel", "Jy/pixel", 2.0),
    ("contrast_pixel", "Jy/arcsec^2", 2.0/0.25),
    ("contrast_pixel", "contrast_peak", 1/0.5),
])
def test_convert_image_values(from_unit, to_unit, factor):
    """convert_image(): pixel area 0.25 arcsec^2, star 2 Jy, PSF peak 0.5 of the star per pixel."""
    out = U.convert_image(np.full((3, 3), 3.0), from_unit, to_unit, 0.25, star_flux=2.0, peak_norm=0.5)
    assert np.allclose(out, 3.0*factor)


def test_convert_image_needs_the_star_flux_and_the_peak():
    """Between families without the star's flux, or to contrast_peak without peak_norm: UnitError."""
    with pytest.raises(U.UnitError, match="contrast_pixel needs the star's flux"):
        U.convert_image(np.ones((3, 3)), "Jy/arcsec^2", "contrast_pixel", 0.25)
    with pytest.raises(U.UnitError, match="contrast_peak needs the peak"):
        U.convert_image(np.ones((3, 3)), "contrast_pixel", "contrast_peak", 0.25)


def test_convert_image_refuses_or_converts_a_non_internal_input_unit():
    """from_unit must be an internal unit: another one is refused (UnitError) or converted correctly."""
    try:
        out = U.convert_image(np.ones((3, 3)), "Jy/pixel", "Jy/arcsec^2", 0.25, star_flux=2.0)
    except U.UnitError:
        return
    assert np.allclose(out, 1/0.25)


def test_flux_per_pixel():
    """flux_per_pixel(): the pixel area for Jy/arcsec^2, 1 for contrast_pixel."""
    assert U.flux_per_pixel("Jy/arcsec^2", 0.25) == 0.25 and U.flux_per_pixel("contrast_pixel", 0.25) == 1.0


def test_star_flux_of():
    """star_flux_of(): None and numbers as they are, a Star at the wavelength; a Star without flux: UnitError."""
    assert U.star_flux_of(None, 1e-6) is None and U.star_flux_of(3, 1e-6) == 3 and U.star_flux_of(2.5, 1e-6) == 2.5
    assert U.star_flux_of(Star(flux_jy=4.0), 1e-6) == 4.0
    with pytest.raises(U.UnitError, match="needs temperature"):
        U.star_flux_of(Star(distance=3.0), 1e-6)


def test_star_flux_of_numpy_scalars():
    """A NumPy number is a flux in Jy too."""
    assert U.star_flux_of(np.float32(2.0), 1e-6) == 2.0 and U.star_flux_of(np.int64(2), 1e-6) == 2


@pytest.mark.parametrize("bunit, expected", [
    ("JY/PIXEL", "Jy/pixel"), ("W.m-2.pixel-1", "W/m^2/pixel"), ("Jy/arcsec**2", "Jy/arcsec^2"), ("MJy/sr", "MJy/sr"),
    ("contrast", "contrast_pixel"), ("ph/s/cm2/pix", U.PHOTONS), ("erg/(s cm2 Hz sr)", "erg/(s cm2 Hz sr)"),
    ("", None), ("counts", None), ("Jy/beam", None),
])
def test_header_unit(bunit, expected):
    """header_unit(): known spellings map to NAMED units, other astropy units are kept, the rest is None."""
    unit, raw = U.header_unit({"BUNIT": bunit} if bunit else {})
    assert unit == expected and raw == bunit


# ====================================================================== sources registry, Source

def test_source_kinds_builtins_first():
    """source_kinds(): the built-in kinds in their order, then the registered ones."""
    kinds = S.source_kinds()
    assert kinds[:len(S.BUILTIN_KINDS)] == list(S.BUILTIN_KINDS) and {"ring", "cme_blob"} <= set(kinds)


def test_get_source_unknown_kind_is_a_key_error():
    """get_source() of an unregistered kind raises KeyError."""
    with pytest.raises(KeyError):
        S.get_source("no_such_kind_xyz")


def test_get_source_imports_the_module(monkeypatch):
    """get_source(kind, module) imports the module, which registers the kind."""
    import sys
    monkeypatch.setattr(S, "_SOURCES", {k: v for k, v in S._SOURCES.items() if k != "ring"})
    monkeypatch.delitem(sys.modules, "extra_sources")
    assert S.get_source("ring", "extra_sources").kind == "ring"


def test_register_source_checks_and_replaces(monkeypatch):
    """register_source(): a kind is required, parameter names must differ; a kind registered again is replaced."""
    monkeypatch.setattr(S, "_SOURCES", dict(S._SOURCES))

    class NoKind(S.Source):
        pass

    class Twice(S.Source):
        kind = "twice_test"
        params = (P.Param("a"), P.Param("a"))

    class Mine(S.Source):
        kind = "ring"

    with pytest.raises(ValueError, match="has no kind"):
        S.register_source(NoKind)
    with pytest.raises(ValueError, match="duplicate parameter names"):
        S.register_source(Twice)
    assert S.register_source(Mine) is Mine and S.get_source("ring") is Mine


def test_source_settings_defaults_and_conversion():
    """Source(**settings): values converted (units in text), the others at their default."""
    ring = Ring(radius="300 mas")
    assert ring.settings == {"radius": 0.3, "width": 0.05, "peak": 1e-3} and ring.missing() == []


@pytest.mark.parametrize("make, words", [
    (lambda: Ring(radiu=0.3), "ring: unknown parameter radiu (did you mean radius?)"),
    (lambda: Ring(radius="3 kg"), "cannot convert to arcsec"),
    (lambda: CmeBlob(mass=-1), "must be > 0"),
])
def test_source_setting_errors(make, words):
    """Unknown names (with a suggestion) and invalid values raise ParamError."""
    with pytest.raises(P.ParamError, match=re.escape(words)):
        make()


def test_source_missing_required_parameter():
    """missing() lists the parameters without a default that are not set."""
    assert [p.name for p in Planet().missing()] == ["pos", "brightness"]


def test_source_symmetry_and_base_render():
    """symmetric(): from the module or the "symmetric" setting ("auto" = None); the base render raises."""
    assert Ring().symmetric() is True and CmeBlob().symmetric() is False and CmeBlob(r=0).symmetric() is True
    assert FitsImage().symmetric() is None and FitsImage(symmetric=False).symmetric() is False
    with pytest.raises(NotImplementedError):
        S.Source().render(None, None)


# ====================================================================== images.py, read_plane

def test_fits_path():
    """fits_path(): "example" is the bundled file; ".fits" is added only without a FITS extension."""
    assert fits_path("example").endswith("example_data/example.fits")
    assert fits_path("disk") == "disk.fits" and fits_path("a.fits.gz") == "a.fits.gz" and fits_path("b.FTS") == "b.FTS"


@pytest.mark.parametrize("shape, kwargs, words", [
    ((2, 5, 5), {}, "the data are 3-D after selecting image (shape (2, 5, 5))"),
    ((2, 3, 5, 5), {"fits_slice": [0]}, "the data are 3-D after selecting fits_slice [0] (shape (3, 5, 5))"),
    ((24,), {}, "the data are 1-D after selecting image (shape (24,))"),
])
def test_read_plane_not_2d_message(in_tmp, shape, kwargs, words):
    """Regression: data that are not 2-D after the selection: the dimension, the selection and the fix."""
    path = _write_fits(in_tmp/"cube.fits", np.zeros(shape))
    with pytest.raises(ValueError) as err:
        read_plane(path, **kwargs)
    assert words in str(err.value) and "fits_slice to pick a 2-D plane" in str(err.value)


def test_read_plane_selects_a_plane_and_reshapes_flat_squares(in_tmp):
    """fits_slice picks a plane of a cube; a flattened square image is reshaped to (n, n) float64."""
    cube = np.arange(2*5*5, dtype=np.float32).reshape(2, 5, 5)
    image, header, fmt, what, index = read_plane(_write_fits(in_tmp/"cube.fits", cube), fits_slice=[1])
    assert np.array_equal(image, cube[1]) and image.dtype == np.float64 and index == [1] and fmt == "image"
    flat, *_ = read_plane(_write_fits(in_tmp/"flat.fits", np.arange(25.0)))
    assert np.array_equal(flat, np.arange(25.0).reshape(5, 5))


@pytest.mark.parametrize("kwargs, words", [
    ({"fits_slice": [0, 0, 0, 0]}, "4 indices for the 3-D data"),
    ({"fits_hdu": 3}, "fits_hdu = 3 does not exist"),
    ({"input_format": "mcfost"}, "is not an MCFOST image cube (expected 5 axes"),
])
def test_read_plane_selection_errors(in_tmp, kwargs, words):
    """Too many indices, a missing HDU and a non-MCFOST file read as MCFOST raise ValueError."""
    path = _write_fits(in_tmp/"cube.fits", np.zeros((2, 5, 5)))
    with pytest.raises(ValueError, match=re.escape(words)):
        read_plane(path, **kwargs)


@pytest.mark.parametrize("cards, scale, scale_origin, lam, lam_origin", [
    ({"PIXSCALE": 0.01, "WAVELEN": 1.6e-6}, 0.01, "header PIXSCALE of d.fits", 1.6e-6, "header WAVELEN of d.fits"),
    ({"CDELT2": -0.01/3600, "WAVE": 1.6}, 0.01, "header CDELT2 (deg) of d.fits", 1.6e-6, "header WAVE (um) of d.fits"),
])
def test_fits_image_inspect_reads_the_header(in_tmp, cards, scale, scale_origin, lam, lam_origin):
    """inspect(): the image size, pixel scale (PIXSCALE or CDELT, abs) and wavelength (WAVELEN or WAVE um)."""
    path = _write_fits(in_tmp/"d.fits", np.zeros((9, 9)), BUNIT="Jy/pixel", **cards)
    found = {k: (v, o) for k, v, o in FitsImage(file=path).inspect()}
    assert found["grid.size"] == (9, "image size of d.fits")
    assert np.isclose(found["grid.pixel_scale"][0], scale) and found["grid.pixel_scale"][1] == scale_origin
    assert np.isclose(found["observation.wavelength"][0], lam) and found["observation.wavelength"][1] == lam_origin


def test_fits_image_settings_override_the_header(in_tmp):
    """input_unit, input_pixel_scale and input_wavelength replace the header's values."""
    path = _write_fits(in_tmp/"d.fits", np.ones((9, 9)), BUNIT="Jy/pixel", PIXSCALE=0.01, WAVELEN=1.6e-6)
    src = FitsImage(file=path, input_unit="contrast_pixel", input_pixel_scale="20 mas", input_wavelength="1 um")
    found = {k: (v, o) for k, v, o in src.inspect()}
    assert found["grid.pixel_scale"] == (0.02, "input_pixel_scale of d.fits")
    assert np.isclose(found["observation.wavelength"][0], 1e-6) and src.family_of() == "relative"
    assert src.render(None, None).unit == "contrast_pixel"


@pytest.mark.parametrize("unit, family", [("auto", None), ("contrast_pixel", "relative"), ("MJy/sr", "physical")])
def test_fits_image_family_from_input_unit(unit, family):
    """family_of(): from input_unit; "auto" leaves it to the header (None)."""
    assert FitsImage(input_unit=unit).family_of() == family


def test_fits_image_file_errors(in_tmp):
    """A missing file, a non-square image and an unknown BUNIT without input_unit are explained."""
    _write_fits(in_tmp/"rect.fits", np.ones((5, 7)), BUNIT="Jy/pixel")
    _write_fits(in_tmp/"nounit.fits", np.ones((5, 5)))
    with pytest.raises(FileNotFoundError, match="file nothere.fits not found"):
        FitsImage(file="nothere").inspect()
    with pytest.raises(ValueError, match="must be square"):
        FitsImage(file="rect").inspect()
    with pytest.raises(ValueError, match="the unit is unknown .* set input_unit"):
        FitsImage(file="nounit").inspect()
    assert FitsImage(file="nounit", input_unit="Jy/pixel").render(None, None).unit == "Jy/pixel"


# ====================================================================== Scene: building and rendering

def test_scene_unknown_setting_refused():
    """Scene(): a misspelled setting or a section not given as a dict raises ParamError."""
    with pytest.raises(P.ParamError, match="unknown setting wavelenght"):
        Scene(wavelenght=1e-6)
    with pytest.raises(P.ParamError, match="unknown setting method"):
        Scene(method="together")
    assert Scene(method={"combine": "together"}, size=9).sections["method"] == {"combine": "together"}


@pytest.mark.parametrize("n", SIZES)
def test_scene_render_shapes_unit_and_total(n):
    """render(): {name: (n, n) float64} in Jy/arcsec^2 on the grid; total is their sum; the module's image as given."""
    scene = Scene(wavelength=1e-6, diameter=2.0, size=n, pixel_scale=0.05)
    scene.add(Ring(radius=0.3, width=0.05, peak=2e-3))
    scene.add(Ring(radius=0.5, peak=1e-3))
    images = scene.render(verbose=False)
    assert list(images) == ["ring", "ring_2"] and scene.internal_unit == "Jy/arcsec^2" and scene.star_flux is None
    assert all(im.shape == (n, n) and im.dtype == np.float64 for im in images.values())
    assert np.allclose(scene.total, images["ring"] + images["ring_2"])
    x, y = scene.grid.coordinates()
    assert np.allclose(images["ring"], 2e-3*np.exp(-0.5*((np.hypot(x, y)-0.3)/0.05)**2), rtol=1e-12)
    assert scene.grid == Grid(n, 0.05, 1e-6, 2.0) and scene.families == {"ring": "physical", "ring_2": "physical"}


def test_scene_automatic_names_continue_after_rendering():
    """Regression: automatic names are per kind and keep counting when components are added later."""
    scene = Scene(size=9, pixel_scale=0.1, star=Star(distance=3.0, flux_jy=2.0))
    scene.add(Ring())
    scene.add(CmeBlob())
    scene.add(Ring())
    scene.render(verbose=False)
    scene.add(Ring())
    assert list(scene.render(verbose=False)) == ["ring", "cme_blob", "ring_2", "ring_3"]


def test_scene_scale_enabled_and_names():
    """scale multiplies an image, a disabled component is left out (its name still counts), names must differ."""
    one = Scene(size=15, pixel_scale=0.05)
    one.add(Ring())
    one.render(verbose=False)
    scene = Scene(size=15, pixel_scale=0.05)
    scene.add(Ring(), enabled=False)
    scene.add(Ring(), scale=2.0)
    assert list(scene.render(verbose=False)) == ["ring_2"] and np.allclose(scene.total, 2*one.total)
    bad = Scene(size=15, pixel_scale=0.05)
    bad.add(Ring(), name="a")
    bad.add(Ring(), name="a")
    assert "two components are named 'a'" in _render_error(bad)
    negative = Scene(size=15, pixel_scale=0.05)
    negative.add(Ring(), scale=-1.0)
    assert "scale = -1.0: must be >= 0" in _render_error(negative)


def test_scene_without_components_has_a_zero_total():
    """No components: no images and a zero (n, n) total."""
    scene = Scene(size=9, pixel_scale=0.1)
    assert scene.render(verbose=False) == {} and np.array_equal(scene.total, np.zeros((9, 9)))


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("pa", [0, 90, 180, 270])
def test_scene_orientation_angles_from_x_rows_are_y(n, pa):
    """A blob at position angle pa (from +x towards +y) lies along (cos pa, sin pa): columns +x, rows +y."""
    scene = Scene(size=n, pixel_scale=0.05, star=Star(distance=10.0))  # 0.5 au per pixel
    scene.add(CmeBlob(r=2.0, pa=pa, width=0.5))                         # 4 pixels from the star
    scene.render(verbose=False)
    dx, dy = _centroid_offset(scene.total)
    along = dx*np.cos(np.radians(pa)) + dy*np.sin(np.radians(pa))
    across = -dx*np.sin(np.radians(pa)) + dy*np.cos(np.radians(pa))
    assert along > 2.0 and abs(across) < 1e-6


def test_scene_contrast_only_without_star_flux():
    """Contrast components without the star's flux stay in contrast_pixel."""
    scene = Scene(size=15, pixel_scale=0.05, star=Star(distance=3.0))
    scene.add(CmeBlob(r=0.3))
    images = scene.render(verbose=False)
    blob = CmeBlob(r=0.3).render(scene.grid, scene.star).data
    assert scene.internal_unit == "contrast_pixel" and scene.star_flux is None and np.allclose(images["cme_blob"], blob)


def test_scene_mixed_families_through_the_star_flux():
    """Contrast and Jy components are combined in Jy/arcsec^2 through the star's flux."""
    scene = Scene(size=15, pixel_scale=0.05, star=Star(distance=3.0, flux_jy=2.0))
    scene.add(CmeBlob(r=0.3))
    scene.add(Ring())
    images = scene.render(verbose=False)
    blob = CmeBlob(r=0.3).render(scene.grid, scene.star).data
    assert scene.internal_unit == "Jy/arcsec^2" and scene.star_flux == 2.0
    assert np.allclose(images["cme_blob"], blob*2.0/0.05**2) and scene.families == {"cme_blob": "relative", "ring": "physical"}


def test_scene_mixed_families_without_star_flux_refused():
    """Contrast and Jy components without the star's flux: ParamError before rendering."""
    scene = Scene(size=15, pixel_scale=0.05, star=Star(distance=3.0))
    scene.add(CmeBlob())
    scene.add(Ring())
    assert "the star's flux is needed" in _render_error(scene)


def test_scene_module_needs_the_distance():
    """A module requiring [star] distance cannot render without it."""
    scene = Scene(size=15, pixel_scale=0.05)
    scene.add(CmeBlob())
    assert "cme_blob needs [star] distance" in _render_error(scene)


def test_scene_quantity_module_equals_the_plain_one():
    """A module returning an astropy Quantity (erg/s/cm^2/Hz/sr) gives the same scene as Jy/arcsec^2."""
    a, b = Scene(size=16, pixel_scale=0.05), Scene(size=16, pixel_scale=0.05)
    a.add(Ring())
    b.add(QuantityRing())
    a.render(verbose=False)
    b.render(verbose=False)
    assert np.allclose(a.total, b.total, rtol=1e-12)


@pytest.mark.parametrize("unit, star, expected", [("Jy", None, 1e-4/0.05**2), ("contrast", Star(flux_jy=2.0), 1e-4*2.0/0.05**2),
                                                  ("contrast", None, 1e-4)])
def test_scene_point_sources_are_kept_apart(unit, star, expected):
    """Planets are positions and fluxes per pixel area in the internal unit, not in the images or total."""
    scene = Scene(size=15, pixel_scale=0.05, star=star)
    scene.add(Planet(pos=[0.2, 0.1], brightness=1e-4, brightness_unit=unit))
    assert scene.render(verbose=False) == {} and not scene.total.any()
    (radius, angle), = scene.points["planet"][0]
    assert np.isclose(radius, np.hypot(0.2, 0.1)) and np.isclose(angle, np.arctan2(0.1, 0.2))
    assert np.allclose(scene.points["planet"][1], [expected])


def test_scene_in_unit():
    """in_unit(): "auto" / None unchanged (same array); other units with the scene's or a given star flux."""
    scene = Scene(size=15, pixel_scale=0.05, star=Star(distance=3.0))
    scene.add(CmeBlob(r=0.3))
    scene.render(verbose=False)
    img = scene.total
    assert scene.in_unit(img, "auto") is img and scene.in_unit(img, None) is img
    assert np.allclose(scene.in_unit(img, "Jy/pixel", star_flux=2.0), img*2.0)
    assert np.allclose(scene.in_unit(img, "Jy/arcsec^2", star_flux=Star(flux_jy=2.0)), img*2.0/0.05**2)
    with pytest.raises(U.UnitError, match="needs the star's flux"):
        scene.in_unit(img, "Jy/pixel")


def test_scene_from_file_keeps_a_lambda_d_field_when_a_component_is_added():
    """A scene from parameters (fov = "8 lambda/D") keeps its grid when a component is added in Python."""
    data = {"observation": {"wavelength": 1e-6, "diameter": 2.0}, "grid": {"sampling": 2.0, "fov": "8 lambda/D"},
            "source": [{"kind": "ring", "module": "extra_sources", "radius": 0.1}]}
    scene = Scene.from_params(P.from_data(data))
    scene.render(verbose=False)
    before = scene.grid
    scene.add(Ring(radius=0.2))
    scene.render(verbose=False)
    assert before.size == 17 and scene.grid == before


# ====================================================================== Scene: image files

def test_scene_example_disk_equals_the_old_reader():
    """The bundled example through the scene equals read_target's image (Jy/arcsec^2), on its header grid."""
    import contextlib
    import io
    from toycoronagraph import config
    from toycoronagraph.inputs import read_target
    scene = Scene()
    scene.add(FitsImage(file="example"))
    scene.render(verbose=False)
    with contextlib.redirect_stdout(io.StringIO()):
        old = read_target(config.from_dict({}))
    assert scene.internal_unit == "Jy/arcsec^2" and scene.grid.size == 512
    assert np.allclose(scene.total, old, rtol=1e-12, atol=0)


@pytest.mark.parametrize("n", SIZES)
def test_scene_fits_image_is_not_flipped(in_tmp, n):
    """A FITS image keeps its orientation: row r, column c of the file is row r, column c of the scene."""
    data = np.zeros((n, n))
    data[3, n-5] = 1.0  # bottom right (rows = +y, as in FITS)
    _write_fits(in_tmp/"pix.fits", data, BUNIT="Jy/pixel", PIXSCALE=0.01, WAVELEN=1e-6)
    scene = Scene()
    scene.add(FitsImage(file="pix"))
    scene.render(verbose=False)
    assert scene.grid == Grid(n, 0.01, 1e-6, 2.4) and np.allclose(scene.total, data/1e-4)
    assert scene.origins["fits_image"] == "pix.fits (image)" and scene.units["fits_image"] == "Jy/pixel"


def test_scene_wavelength_mismatch_warns(in_tmp):
    """An image made at another wavelength than the scene's: a warning naming both."""
    path, _ = _disk(in_tmp, "disk", 16, 0.01, wavelength=1.6e-6)
    scene = Scene(wavelength=1.0e-6)
    scene.add(FitsImage(file=path))
    scene.render(verbose=False)
    warning = next(n for n in scene.notes if n.startswith("WARNING"))
    assert "observation.wavelength" in warning and "1e-06" in warning and "1.6e-06" in warning


@pytest.mark.parametrize("n", SIZES)
def test_scene_default_resamples_another_field_with_a_note(in_tmp, n):
    """By default an image of another pixel scale is resampled (flux kept when the grid covers it), with a note."""
    path, image = _disk(in_tmp, "disk", n, 0.01)
    scene = Scene(wavelength=1.6e-6, pixel_scale=0.013, size=n)  # 0.013 n arcsec > the image's 0.01 n
    scene.add(FitsImage(file=path))
    scene.render(verbose=False)
    assert any(note.startswith("Note: resampled (flux conserved) disk.fits") for note in scene.notes)
    assert abs(scene.total.sum()*0.013**2/(image.sum()*0.01**2) - 1) < 1e-12
    assert np.allclose(_centroid_offset(scene.total), 0, atol=1e-9)


def test_scene_write_and_read_back(in_tmp):
    """write() then kind "scene_fits": the same components, grid and images (one extension each)."""
    scene = Scene(wavelength=1e-6, diameter=2.0, size=17, pixel_scale=0.05, star=Star(distance=3.0, flux_jy=2.0))
    scene.add(CmeBlob(r=0.3))
    scene.add(Ring(radius=0.3))
    scene.render(verbose=False)
    scene.write(str(in_tmp/"scene.fits"))
    header = fits.getheader(in_tmp/"scene.fits", 1)
    assert header["EXTNAME"] == "CME_BLOB" and header["CRPIX1"] == 9 and header["BUNIT"] == "Jy/arcsec^2"
    back = Scene(diameter=2.0)
    back.add(SceneFits(file=str(in_tmp/"scene.fits")))
    back.render(verbose=False)
    assert list(back.images) == ["cme_blob", "ring"] and back.grid == scene.grid
    assert all(np.allclose(back.images[k], scene.images[k], rtol=1e-12, atol=0) for k in scene.images)
    assert back.params.get("star.distance") == 3.0 and back.params.origin("star.distance") == "header DISTANCE of scene.fits"


def test_scene_plot_has_one_panel_per_component_and_the_total(in_tmp):
    """plot(): one image panel per component plus the total; saved to path, not shown."""
    import matplotlib.pyplot as plt
    scene = Scene(size=16, pixel_scale=0.05, star=Star(distance=3.0, flux_jy=2.0))
    scene.add(Ring())
    scene.add(CmeBlob(r=0.3))
    fig = scene.plot(log=True, path=str(in_tmp/"scene.png"), show=False)
    assert [ax.get_title() for ax in fig.axes if ax.images] == ["ring", "cme_blob", "total"]
    assert (in_tmp/"scene.png").exists()
    plt.close(fig)


# ====================================================================== [grid] resample = false

@pytest.mark.parametrize("n", SIZES)
def test_strict_grid_matching_image_is_used_as_is(in_tmp, n):
    """resample = false and an image on the grid: no error, no note, the image unchanged (to round-off)."""
    path, image = _disk(in_tmp, "a", n, 0.01)
    scene = Scene(wavelength=1.6e-6, resample=False)
    scene.add(FitsImage(file=path))
    assert _render_error(scene) is None and np.allclose(scene.total, image, rtol=1e-12, atol=0)
    assert not [note for note in scene.notes if "resampled" in note or "rebinned" in note]


@pytest.mark.parametrize("n", SIZES)
def test_strict_grid_same_field_other_resolution_is_rebinned(in_tmp, n):
    """resample = false with only [grid] size: the image's field is kept, rebinned, flux conserved, with a note."""
    path, image = _disk(in_tmp, "a", n, 0.01)
    half = (n+1)//2
    scene = Scene(wavelength=1.6e-6, resample=False, size=half)
    scene.add(FitsImage(file=path))
    assert _render_error(scene) is None and np.isclose(scene.grid.pixel_scale*half, 0.01*n)
    assert abs(scene.total.sum()*scene.grid.pixel_area/(image.sum()*1e-4) - 1) < 1e-12
    assert any("rebinned to the grid's resolution" in note for note in scene.notes)


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("grid, words", [
    ({"pixel_scale": 0.013}, "{n} px of 0.01 arcsec"),          # another pixel scale
    ({"size": "crop", "pixel_scale": 0.01}, "the grid is {crop} px"),  # a crop
])
def test_strict_grid_other_field_is_refused_with_the_fix(in_tmp, n, grid, words):
    """resample = false and an image covering another field: ParamError with the image, the grid and the fix."""
    path, _ = _disk(in_tmp, "a", n, 0.01)
    grid = {k: (n-10 if v == "crop" else v) for k, v in grid.items()}
    scene = Scene(wavelength=1.6e-6, resample=False, **grid)
    scene.add(FitsImage(file=path))
    error = _render_error(scene)
    assert error and "resample = false" in error and "set only [grid] size" in error and "a.fits (image)" in error
    assert words.format(n=n, crop=n-10) in error


def test_strict_grid_lists_every_mismatch_at_once(in_tmp):
    """Several images off the grid's field are listed together in one error."""
    a, _ = _disk(in_tmp, "a", 24, 0.01)
    b, _ = _disk(in_tmp, "b", 24, 0.012)
    scene = Scene(wavelength=1.6e-6, resample=False, size=24, pixel_scale=0.02)
    scene.add(FitsImage(file=a))
    scene.add(FitsImage(file=b))
    error = _render_error(scene)
    assert "a.fits (image): 24 px of 0.01 arcsec" in error and "b.fits (image): 24 px of 0.012 arcsec" in error


def test_strict_grid_ignores_header_round_off_and_computed_modules(in_tmp):
    """A pixel scale off by 1e-7 (a rounded header) and modules computed on the grid are fine."""
    path, _ = _disk(in_tmp, "rounded", 25, 0.01*(1+1e-7))
    scene = Scene(wavelength=1.6e-6, resample=False, size=25, pixel_scale=0.01, star=Star(distance=3.22, flux_jy=2.0))
    scene.add(FitsImage(file=path))
    scene.add(Ring(radius=0.05))
    scene.add(CmeBlob(r=0.5))
    assert _render_error(scene) is None and list(scene.images) == ["fits_image", "ring", "cme_blob"]


def test_strict_grid_from_a_parameter_file(in_tmp):
    """[grid] resample = false in a parameter file stops a mismatching image too."""
    path, _ = _disk(in_tmp, "a", 24, 0.01)
    with open(in_tmp/"strict.toml", "w") as f:
        f.write('[grid]\npixel_scale = 0.013\nresample = false\n[[source]]\nkind = "fits_image"\nfile = "%s"\n' % path)
    assert "resample = false" in (_render_error(Scene.from_file(str(in_tmp/"strict.toml"))) or "")
