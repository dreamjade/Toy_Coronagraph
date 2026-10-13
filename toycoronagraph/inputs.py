"""
Reading the target image from a FITS file.

read_target reads it as described by the older parameter layout (toycoronagraph.config,
used by Target); read_plane, the FITS reader itself, is also used by the image sources of
the scene framework (toycoronagraph.sources.images).

Two input methods:
    * "image": a square 2-D FITS image (or a cube with fits_slice), star at the array centre,
      described by the header keywords BUNIT, PIXSCALE (or CDELT1/2) and WAVELEN;
    * "mcfost": an MCFOST image cube (e.g. data_<wavelength>/RT.fits.gz) as written by
      MCFOST, with axes (plane, azimuth, inclination, y, x); the plane is chosen by name
      (mcfost_component; "disk" is the total minus the direct starlight) and the
      inclination/azimuth by index. Its header (BUNIT, CDELT, WAVE) is understood.
"input_format = auto" recognizes MCFOST cubes from their header (FLUX_1 ... FLUX_8).
"""
import os
import numpy as np
from astropy.io import fits
from toycoronagraph import DATADIR
from toycoronagraph.config import AUTO, MCFOST_COMPONENTS, resolve


def target_path(cfg):
    """
    The FITS file of the target: the example, or fits_filename (".fits" added if it has no extension).

    Args:
        cfg (SimpleNamespace): From config.load_config (uses cfg.fits_filename; None for the example).

    Returns:
        str: Path of the FITS file (example_data/example.fits for the example).
    """
    if cfg.fits_filename is None:
        return os.path.join(DATADIR, "example.fits")
    name = cfg.fits_filename
    if name.lower().endswith((".fits", ".fits.gz", ".fit", ".fts")):
        return name
    return name+".fits"


def is_mcfost(data, header):
    """
    Whether the data look like an MCFOST image cube.

    Args:
        data (np.ndarray): The HDU's data.
        header: Its FITS header.

    Returns:
        bool: True for 5 axes with 4 or 8 planes on axis 0 and FLUX_1 in the header.
    """
    return data.ndim == 5 and data.shape[0] in (4, 8) and "FLUX_1" in header


def read_plane(path, input_format=AUTO, fits_hdu=0, fits_slice=AUTO, mcfost_component="scattered_star",
               mcfost_inclination=0, mcfost_azimuth=0):
    """
    Read one 2-D image from a FITS file (a plane of a cube, or of an MCFOST cube).

    The plane is chosen by fits_slice if given, else by the MCFOST keys for an MCFOST
    cube; a plain image is read as it is. A flattened square image (1-D, n*n values) is
    reshaped to (n, n).

    Args:
        path (str): The FITS file.
        input_format (str): "image", "mcfost", or "auto" (MCFOST if is_mcfost, else image).
        fits_hdu (int): Index of the HDU holding the data.
        fits_slice (list or str): Indices into the leading axes, e.g. [5, 0, 0] (overrides the
            MCFOST keys), or "auto".
        mcfost_component (str): MCFOST plane: a key of config.MCFOST_COMPONENTS, or "disk"
            (total minus direct_star).
        mcfost_inclination (int): MCFOST: index along the inclination axis.
        mcfost_azimuth (int): MCFOST: index along the azimuth axis.

    Returns:
        tuple:
            image (np.ndarray): (ny, nx) float64 data, in the file's unit (not checked to be square).
            header (astropy.io.fits.Header): The HDU's header.
            fmt (str): The input format used, "image" or "mcfost".
            what (str): Description of the selection, e.g. "image" or "MCFOST disk, inclination 0, azimuth 0".
            index (list or str): The indices used (list of int); for "disk", a text such as
                "total - direct_star at [0, 0]".

    Raises:
        ValueError: For a missing HDU, an MCFOST format without 5 axes, an out-of-range MCFOST
            index, too many fits_slice indices, or data that are not 2-D after the selection.
    """
    with fits.open(path) as hdul:
        if fits_hdu >= len(hdul):
            raise ValueError("%s has %d HDUs; fits_hdu = %d does not exist" % (path, len(hdul), fits_hdu))
        data = np.asarray(hdul[fits_hdu].data, dtype=np.float64)
        header = hdul[fits_hdu].header
    fmt = input_format
    if fmt == AUTO:
        fmt = "mcfost" if is_mcfost(data, header) else "image"

    if fits_slice != AUTO:
        index = list(fits_slice)
        what = "fits_slice %s" % index
    elif fmt == "mcfost":
        if data.ndim != 5:
            raise ValueError("%s is not an MCFOST image cube (expected 5 axes: plane, azimuth, inclination, y, x; got shape %s)"
                             % (path, data.shape))
        component = MCFOST_COMPONENTS["total" if mcfost_component == "disk" else mcfost_component]
        # Range checks; "disk" also needs plane 4 (direct_star), so that index is checked for it
        for name, value, size in [("mcfost_component", max(component, 4 if mcfost_component == "disk" else 0), data.shape[0]),
                                  ("mcfost_azimuth", mcfost_azimuth, data.shape[1]),
                                  ("mcfost_inclination", mcfost_inclination, data.shape[2])]:
            if not 0 <= value < size:
                raise ValueError("%s = %r is out of range: %s has %d along that axis" % (name, value, path, size))
        index = [component, mcfost_azimuth, mcfost_inclination]
        what = "MCFOST %s, inclination %d, azimuth %d" % (mcfost_component, mcfost_inclination, mcfost_azimuth)
    else:
        index = []
        what = "image"
    if len(index) > data.ndim:
        raise ValueError("%d indices for the %d-D data in %s" % (len(index), data.ndim, path))
    image = data[tuple(index)]
    if fmt == "mcfost" and fits_slice == AUTO and mcfost_component == "disk":
        image = image-data[MCFOST_COMPONENTS["direct_star"], mcfost_azimuth, mcfost_inclination] # remove the star
        index = "total - direct_star at %s" % index[1:]
    if image.ndim == 1 and round(np.sqrt(image.size))**2 == image.size: # a flattened square image
        image = image.reshape(int(round(np.sqrt(image.size))), -1)
    if image.ndim != 2:
        raise ValueError("%s: the data are %d-D after selecting %s (shape %s); set input_format = \"mcfost\" for an MCFOST cube, "
                         "or fits_slice to pick a 2-D plane" % (path, image.ndim, what, image.shape))
    return image, header, fmt, what, index


def read_target(cfg, verbose=True):
    """
    Read the target image, fill the "auto" parameters from it (config.resolve), and convert
    it to Jy/arcsec^2.

    Args:
        cfg (SimpleNamespace): From config.load_config; resolved in place (also sets
            input_description and used_slice).
        verbose (bool): Print the summary of the input and any warning.

    Returns:
        np.ndarray: The (N, N) float64 image in Jy/arcsec^2.

    Raises:
        ValueError: From read_plane or config.resolve (e.g. non-square image, unknown unit).
    """
    path = target_path(cfg)
    image, header, fmt, what, index = read_plane(path, cfg.input_format, cfg.fits_hdu, cfg.fits_slice, cfg.mcfost_component,
                                                 cfg.mcfost_inclination, cfg.mcfost_azimuth)

    cfg.input_description = "%s (%s)" % (os.path.basename(path), what)
    cfg.used_slice = index
    notes = resolve(cfg, image.shape, header)
    if verbose:
        for note in notes:
            print(note)

    # Convert to Jy/arcsec^2
    if cfg.input_unit == "W/m^2/pixel": # nu F_nu per pixel
        return image*cfg.F_transfer
    if cfg.input_unit == "Jy/pixel":
        return image/cfg.psf_scale**2
    return image
