"""Planets as point sources: each goes through the coronagraph with its exact PSF."""
import numpy as np
from toycoronagraph.params import Param
from toycoronagraph.planet import planet_position
from toycoronagraph.tool import convert_to_polar
from toycoronagraph.sources import Source, register_source


@register_source
class Planet(Source):
    """
    One planet, an unresolved source (point_sources(); no image is rendered).

    family is "physical" (brightness in Jy) or "relative" (a contrast), set by brightness_unit.
    """
    kind = "planet"
    doc = "a planet (static or on an orbit)"
    family = "physical"
    params = (
        Param("mode", "cartesian", doc='"cartesian" [x, y], "polar" [r, deg] or "moving" orbit',
              choices=("cartesian", "polar", "moving")),
        Param("pos", None, doc="position: [x, y] or [r, angle] in arcsec/deg, or [a, e, PA, inc, t]",
              type="float_list", example=[0.5, 0.0]),
        Param("brightness", None, doc="in Jy (or in contrast, see brightness_unit)", min=0, example=5e-5),
        Param("brightness_unit", "Jy", doc='"Jy", or "contrast" (a fraction of the star\'s flux)',
              choices=("Jy", "contrast"), advanced=True),
    )

    def family_of(self):
        """
        The family given by brightness_unit.

        Returns:
            str: "relative" for "contrast", else "physical".
        """
        return "relative" if self.settings["brightness_unit"] == "contrast" else "physical"

    def position(self):
        """
        The planet's position from the star.

        mode "cartesian": pos = [x, y] in arcsec; "polar": pos = [r arcsec, angle deg]; "moving":
        pos = [a arcsec, e, PA deg, inc deg, t] for planet.planet_position (t: the fraction of
        the orbital period, its integer part dropped).

        Returns:
            tuple: (radius in arcsec, angle in rad from +x towards +y).
        """
        s = self.settings
        pos = s["pos"]
        if s["mode"] == "moving":
            return planet_position(*pos[:5], mode="polar", res=1000)
        if s["mode"] == "polar":
            return pos[0], np.radians(pos[1])
        return tuple(convert_to_polar(pos))

    def point_sources(self, grid, star):
        """
        The planet as one point source.

        Args:
            grid (toycoronagraph.scene.Grid): not used.
            star (toycoronagraph.scene.Star): not used.

        Returns:
            list: [((radius arcsec, angle rad), flux)]: flux in Jy, or the fraction of the star's
                flux (brightness_unit "contrast").
        """
        return [(self.position(), self.settings["brightness"])]
