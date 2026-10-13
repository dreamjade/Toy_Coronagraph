"""
Starlight scattered by the free electrons of a stellar wind (Thomson scattering).

A spherical wind of mass-loss rate Mdot and speed v has an electron density
n(r) = Mdot / (4 pi r^2 v mu_e m_u). Each electron scatters sigma_T times the stellar flux
L / (4 pi r^2). Integrated along the line of sight (int dz / r^4 = pi / (2 b^3)), the
scattered light per solid angle, as a fraction of the star's flux, is

    contrast per steradian = K_w / theta^3,   K_w = sigma_T Mdot / (32 pi v mu_e m_u d),

with theta the angle from the star (rad) and d the distance. As in wind_image.ipynb this is
isotropic scattering at constant speed, with mu_e = 1. It is independent of the luminosity,
the band and the wavelength: grey scattering gives a contrast valid in any band.

Each pixel holds the exact integral of K_w / theta^3 over its area, from the corners:
over the quadrant beyond (x, y) >= 0, the integral of r^-3 is 2 / (x + y + r), finite
everywhere but at the star. The wind starts at theta_in = r_in R_star / d; a rectangle
[0, b] x [0, c] with a corner on the star holds pi / (2 theta_in) - sqrt(b^2 + c^2) / (b c).

Options (off by default): mu_e; the Thomson phase function (3/16 pi)(1 + cos^2 chi),
exactly 15/16 of the isotropic result at constant speed; and a beta-law wind
v(r) = v0 + (v - v0)(1 - r_in / r)^beta (line-of-sight integral done numerically).
"""
import numpy as np
from scipy import integrate
from astropy import units as u, constants as const
from toycoronagraph.params import Param
from toycoronagraph.sources import Source, register_source
from toycoronagraph.units import SourceImage

ARCSEC = (1*u.arcsec).to_value(u.rad)


def wind_constant(mdot, v, distance, mu_e=1.0):
    """
    The wind constant K_w = sigma_T Mdot / (32 pi v mu_e m_u d).

    Args:
        mdot (float): mass-loss rate in M_sun/yr.
        v (float): wind speed in km/s.
        distance (float): the star's distance in pc.
        mu_e (float): mass per electron in atomic mass units.

    Returns:
        float: K_w, such that the contrast per steradian is K_w / theta^3 (theta in rad).
    """
    k = const.sigma_T*mdot*u.M_sun/u.yr/(32*np.pi*v*u.km/u.s*mu_e*const.u*distance*u.pc)
    return float(k.decompose().value)


def _beyond(x, y):
    """
    The integral of r^-3 over the quadrant [x, inf) x [y, inf), x, y >= 0 (not both 0).

    Args:
        x, y (float or np.ndarray): corner of the quadrant (any unit; the result is in its inverse).

    Returns:
        float or np.ndarray: 2 / (x + y + hypot(x, y)), broadcast.
    """
    return 2.0/(x + y + np.hypot(x, y))


def pixel_integrals(size, pixel, theta_in):
    """
    The integral of theta^-3 over every pixel of a size x size grid centred on the star.

    Each pixel is split along x = 0 and y = 0 into rectangles [a, b] x [c, d] of the first
    quadrant (by symmetry), integrated exactly from _beyond() at their corners; a rectangle with
    a corner on the star holds pi / (2 theta_in) - hypot(b, d) / (b d) (exact while the disc
    theta < theta_in lies inside it). Exact when theta_in is below half a pixel; otherwise the
    pixels within theta_in of the star are integrated numerically (_numerical_near_star).

    Args:
        size (int): pixels on a side; the star is at pixel (size-1)/2 (between pixels if even).
        pixel (float): pixel size in rad (or any angle unit, the same as theta_in's).
        theta_in (float): angle (rad) below which there is no wind (> 0).

    Returns:
        tuple: (np.ndarray (size, size) float64 in rad^-1 (the integral of theta^-3 d^2theta, in
            the inverse of the unit of pixel), bool: whether the numerical fallback was used).
    """
    edges = (np.arange(size+1) - size/2)*pixel
    lo, hi = edges[:-1], edges[1:]
    # each pixel's [lo, hi] split into its parts at x >= 0 and x <= 0, both mirrored to [0, inf) (empty: lo = hi)
    parts = [(np.maximum(lo, 0), np.maximum(hi, 0)), (np.maximum(-hi, 0), np.maximum(-lo, 0))]
    out = np.zeros((size, size))
    for (y1, y2) in parts:          # rows: y
        for (x1, x2) in parts:      # columns: x
            a, b = x1[None, :], x2[None, :]
            c, d = y1[:, None], y2[:, None]
            empty = (b <= a) | (d <= c)
            corner = (a == 0) & (c == 0) & ~empty   # touches the star
            with np.errstate(divide="ignore", invalid="ignore"):
                rect = _beyond(a, c) - _beyond(a, d) - _beyond(b, c) + _beyond(b, d)
                at_star = np.pi/(2*theta_in) - np.hypot(b, d)/(b*d)
            out += np.where(empty, 0.0, np.where(corner, at_star, rect))
    fallback = theta_in >= pixel/2
    if fallback:
        out = _numerical_near_star(out, size, pixel, theta_in)
    return out, fallback


def _numerical_near_star(out, size, pixel, theta_in, sub=400):
    """
    Replace the pixels near the star by a numerical integral of theta^-3, leaving out theta < theta_in.

    The pixels whose centre is within theta_in + one pixel of the star are summed on a
    sub x sub grid of points (midpoint rule).

    Args:
        out (np.ndarray): (size, size) exact integrals from pixel_integrals; changed in place.
        size (int): pixels on a side (star at pixel (size-1)/2).
        pixel (float): pixel size in rad.
        theta_in (float): inner edge of the wind in rad.
        sub (int): sub-samples per pixel side.

    Returns:
        np.ndarray: out, in rad^-1.
    """
    c = (size-1)/2
    rows, cols = np.indices((size, size))
    near = np.hypot(rows-c, cols-c)*pixel < theta_in + pixel
    t = (np.arange(sub)+0.5)/sub - 0.5
    for row, col in zip(*np.nonzero(near)):
        x = ((col-c) + t)[None, :]*pixel
        y = ((row-c) + t)[:, None]*pixel
        r = np.hypot(x, y)
        out[row, col] = np.sum(np.where(r >= theta_in, r**-3.0, 0.0))*(pixel/sub)**2
    return out


def profile_factor(b, r_in, v, v0, beta, phase):
    """
    The line-of-sight integral at impact parameter b over the isotropic, constant-speed one (pi / (2 b^3)).

    With b = r cos t (t the angle from the plane of the sky; z = b tan t), int dz / r^4 =
    b^-3 int cos^2 t dt. The density scales as 1 / v(r), so a beta law weights the integrand by
    v / v(r); the Thomson phase function (scattering angle chi with cos chi = sin t) by
    (3/4)(1 + sin^2 t) relative to isotropic. The integral is done over t in [0, pi/2] (half the
    line of sight, symmetric), normalised by its constant-speed isotropic value pi/4.

    Args:
        b (float): impact parameter, in the same unit as r_in (rad in Wind).
        r_in (float): where the wind starts (the beta law's r_in; below it v(r) = v0).
        v (float): terminal speed (km/s; only ratios matter).
        v0 (float): speed at r_in (same unit as v).
        beta (float): exponent of the beta law; None = constant speed v.
        phase (str): "isotropic" or "thomson".

    Returns:
        float: the factor (1 for constant speed and isotropic scattering).
    """
    def integrand(t):
        """The integrand at angle t (rad): cos^2 t x phase weight x v / v(r)."""
        r = b/np.cos(t)
        speed = v if beta is None else v0 + (v - v0)*max(1 - r_in/r, 0.0)**beta
        weight = 0.75*(1 + np.sin(t)**2) if phase == "thomson" else 1.0
        return np.cos(t)**2*weight*v/speed
    return 4/np.pi*integrate.quad(integrand, 0, np.pi/2, limit=200)[0]


@register_source
class Wind(Source):
    """
    A spherical stellar wind seen in Thomson-scattered starlight (see the module docstring).

    family "relative": a contrast, valid in any band; requires the star's distance (pc) and
    radius (R_sun, for where the wind starts).
    """
    kind = "wind"
    doc = "starlight scattered by electrons in a stellar wind"
    family = "relative"
    requires = ("star.distance", "star.radius")
    params = (
        Param("mdot", None, unit="M_sun/yr", doc="mass-loss rate", min=0, exclusive_min=True, example=6e-13),
        Param("v", 400.0, unit="km/s", doc="wind speed (constant unless velocity_law = \"beta\")", min=0,
              exclusive_min=True),
        Param("mu_e", 1.0, doc="mass per electron, in atomic mass units (notebooks: 1; solar: ~1.17)", min=0,
              exclusive_min=True, advanced=True),
        Param("r_in", 1.0, doc="where the wind starts, in stellar radii", min=1, advanced=True),
        Param("velocity_law", "constant", doc='"constant", or "beta": v0 + (v - v0)(1 - r_in/r)^beta',
              choices=("constant", "beta"), advanced=True),
        Param("beta", 1.0, doc="exponent of the beta velocity law", min=0, advanced=True),
        Param("v0", 10.0, unit="km/s", doc="speed at r_in (beta law)", min=0, exclusive_min=True, advanced=True),
        Param("phase_function", "isotropic", doc='"isotropic" (notebooks) or "thomson" (1+cos^2)',
              choices=("isotropic", "thomson"), advanced=True),
    )

    def symmetric(self):
        """Always True: a spherical wind is circularly symmetric."""
        return True

    def render(self, grid, star):
        """
        The wind's contrast per pixel: K_w x the integral of theta^-3 over each pixel (x the
        profile factor with a beta law or the Thomson phase function).

        Args:
            grid (toycoronagraph.scene.Grid): the scene grid (size, pixel_scale in arcsec).
            star (toycoronagraph.scene.Star): distance (pc) and radius (R_sun).

        Returns:
            SourceImage: (size, size) contrast_pixel on the scene grid (pixel_scale None), with
                notes for an even grid size or the numerical fallback.
        """
        s = self.settings
        pixel = grid.pixel_scale*ARCSEC
        theta_in = s["r_in"]*(star.radius*u.R_sun/(star.distance*u.pc)).decompose().value
        k = wind_constant(s["mdot"], s["v"], star.distance, s["mu_e"])
        image, fallback = pixel_integrals(grid.size, pixel, theta_in)
        notes = []
        if grid.size % 2 == 0:
            notes.append("even grid size: no pixel is on the star, so the four central pixels hold the wind next to the "
                         "star at 0.7 pixel from it; an odd size is better for winds")
        if fallback:
            notes.append("the wind starts at %.3g pixel from the star (more than half a pixel): the pixels next to it are "
                         "integrated numerically" % (theta_in/pixel))
        if s["velocity_law"] == "beta" or s["phase_function"] == "thomson":
            image *= self._factor_map(grid, pixel, theta_in)
        return SourceImage(k*image, "contrast_pixel", origin="wind (Mdot %.3g M_sun/yr, v %g km/s)" % (s["mdot"], s["v"]),
                           notes=tuple(notes))

    def _factor_map(self, grid, pixel, theta_in, sub=8):
        """
        Each pixel's mean of profile_factor, weighted by theta^-3 over sub x sub points (the
        integral itself stays exact); the factor is tabulated against the impact parameter.

        Args:
            grid (toycoronagraph.scene.Grid): the scene grid.
            pixel (float): pixel size in rad.
            theta_in (float): inner edge of the wind in rad (the beta law's r_in).
            sub (int): sub-samples per pixel side.

        Returns:
            np.ndarray: (size, size) float64 factor; 0 for pixels entirely within theta_in.
        """
        s = self.settings
        beta = s["beta"] if s["velocity_law"] == "beta" else None
        b_max = grid.size*pixel
        table_b = theta_in*np.geomspace(1+1e-6, b_max/theta_in*1.5+2, 400)   # beyond the grid's corners
        table_f = np.array([profile_factor(b, theta_in, s["v"], s["v0"], beta, s["phase_function"]) for b in table_b])
        c = (grid.size-1)/2
        t = (np.arange(sub)+0.5)/sub-0.5
        offsets = (np.arange(grid.size)-c)[:, None] + t[None, :]          # (size, sub) in pixels
        x = offsets.reshape(1, 1, grid.size, sub)*pixel                    # r: (row, sub y, column, sub x)
        y = offsets.reshape(grid.size, sub, 1, 1)*pixel
        r = np.hypot(x, y)
        w = np.where(r >= theta_in, r**-3.0, 0.0)
        f = np.interp(r, table_b, table_f, left=0.0)
        num, den = (w*f).sum(axis=(1, 3)), w.sum(axis=(1, 3))
        with np.errstate(invalid="ignore"):
            return np.where(den > 0, num/den, 0.0)
