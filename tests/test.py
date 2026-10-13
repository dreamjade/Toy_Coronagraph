# -------------------- Setup --------------------

# Parameters are read from toycoronagraph_para.toml in the working directory; if it does not exist, a
# minimal one is created. Keys left "auto" are filled from the target's FITS header; the values used are
# saved to toycoronagraph_para_used.toml (the parameter file itself is not modified).
import toycoronagraph.main as toy

# -------------------- The example target --------------------

# Load the example target image (example.fits) in Toy_Coronagraph/toycoronagraph/example_data/ folder
toy_target = toy.Target()

'''
Input: example.fits (image), 512 x 512, W/m^2/pixel; 0.0053715 arcsec/pixel (header PIXSCALE of example.fits) = 32 lambda/D across at 1e-06 m (header WAVELEN of example.fits)
'''

# Plot the preimage on a log colour scale; the image will auto-save to origin.png
toy_target.plot_origin(log=True)

# Which method makes the final image? "auto" (the default) uses the circular method when the target and
# the PSFs are both circularly symmetric about the star, and the spatial method otherwise
method, spatial_mode = toy_target.check_symmetry()

'''
Method: circular (auto; target symmetric: asymmetry 0.13%, light in the corners 0% (tol 1%); PSFs symmetric: vortex model)
'''

# Plot the final image through vortex coronagraph with charge = 2; the image will auto-save to charge2_final.png
# (the radial PSFs are computed once and cached in [compute] psf_cache_dir)
final_image = toy_target.plot_final(charge=2)

# Plot the final image through vortex coronagraph with charge = 6; the image will auto-save to charge6_final.png
final_image = toy_target.plot_final(charge=6)

# -------------------- Planets --------------------

# Now add a static planet, where pos=(x_position, y_position) in arcsec
toy_target.add_planet(pos=[0.5,0], brightness=0.00005, mode="cartesian")

# You could also add a moving planet on an elliptic orbit, where pos=(length of the semi-major axis, eccentricity, position angle, inclination angle, time/period))
toy_target.add_planet(pos=[0.5,0.6,60,0,0.1], brightness=0.00003, mode="moving")

# Or a static planet in polar coordinates, where pos=(r in arcsec, angle in degrees)
toy_target.add_planet(pos=[0.3,45], brightness=0.00004, mode="polar")

# Plot the preimage with planets; the image will auto-save to origin_with_planets.png
toy_target.plot_origin(log=True)

# Plot the final image with planets; the image will auto-save to charge6_with_planets_final.png
final_image = toy_target.plot_final(charge=6)

# Plot the final image by turning off the target inside IWA; the image will auto-save to charge6_with_planets_iwa_ignore_final.png
final_image = toy_target.plot_final(charge=6, iwa_ignore=True)

# List the planets
toy_target.list_planets()

'''
Static Planet 1: (0.5, 0.0) arcsec, brightness: 5.00e-05 Jy
Moving Planet 2: (-0.3751700315402439, 0.08996213698601588) arcsec, brightness: 3.00e-05 Jy
Static Planet 3: (0.21213203435596426, 0.21213203435596423) arcsec, brightness: 4.00e-05 Jy
'''

# Delete a specific planet
toy_target.delete_planet(order=1)

'''
Successfully remove planet #1, here is the latest planet list:
Moving Planet 1: (-0.3751700315402439, 0.08996213698601588) arcsec, brightness: 3.00e-05 Jy
Static Planet 2: (0.21213203435596426, 0.21213203435596423) arcsec, brightness: 4.00e-05 Jy
'''

# Move a moving planet on its orbit by time/period=0.1
toy_target.planet_move(time=0.1, order=1)

'''
Planet has moved to new position
'''

# Plot the orbit of a specific moving planet; the image will auto-save to oribit_planet#.png
toy_target.plot_orbit(order=1)

# Plot the new preimage with planets; the image will auto-save to origin_with_planets.png
toy_target.plot_origin(log=True)

# Plot the new final image with planets; the image will auto-save to charge6_with_planets_final.png
final_image = toy_target.plot_final(charge=6)

# -------------------- Contrast, throughput and videos --------------------

# Print a specific planet brightness, background brightness, and background brightness (ignored dust inside IWA)
# in the final image through vortex coronagraph with charge = 6 (also returned as a tuple)
planet_b, dust_b, dust_b_iwa = toy_target.contrast(charge=6, order=1)

'''
This is a charge-6 vortex coronagraph.
Planet brightness: 1.11e-05 Jy
Dust brightness: 1.68e-05 Jy
Dust brightness but ignore dust inside IWA: 1.64e-05 Jy
'''

# Plot the core throughput and the inner working angle, in lambda/D
import toycoronagraph.psf as psf
psfs = toy_target.radial_psfs(charge=2)  # the radial PSFs (computed once, cached in psf_cache_dir)
lambda_d_per_pixel = 2*toy_target.par.psf_range/toy_target.par.target_pixel  # the image's pixel, in lambda/D
psf.cir_core_throughput_plot(psfs, pixel_scale=lambda_d_per_pixel)

# Make the planet movie; the video will be saved to planet_video.webm
toy_target.planet_video(charge=6)

# Make the planet movie again, but ignore the dust inside IWA; the video will be saved to planet_video_iwa_ignore.webm
toy_target.planet_video(charge=6, iwa_ignore=True)

# -------------------- The spatial method --------------------

# The spatial method works for any target (it can also be set with method = "spatial" in the parameter
# file). It chooses the kernel size for [method] error_target, prints an error analysis per stage and the
# estimated error of each image; the image will auto-save to charge6_with_planets_spatial_final.png
final_image = toy_target.plot_final(charge=6, method="spatial")

'''
Method: spatial, rotate mode (set in the parameters; target symmetric: asymmetry 0.13%, ...)
Kernel size K=213 (+-6.6 lambda/D): estimated crop error 1.6% of the image light (80% of the 2% error target; 43 probe pixels)
rotate mode: 9495 PSF samples (5088 per-pixel within 2.52 lambda/D, 44 rings), K=213, kernel cube 1.72 GB
Spatial PSF error analysis vs exact PSFs (...):  [a table per stage and distance from the star]
Estimated error of this image (flux-weighted, from the PSF check): typical 1.8%, up to 2.9%; largest part: kernel crop (K=213) (1.7%)
'''

# Contrast with the spatial method (the planet PSF is computed exactly at the planet's position)
planet_b, dust_b, dust_b_iwa = toy_target.contrast(charge=6, order=1, method="spatial")

# Planet movie with the spatial method; the video will be saved to planet_video_spatial.webm
toy_target.planet_video(charge=6, method="spatial")

# -------------------- Your own target: an MCFOST disk --------------------

# An MCFOST image cube can be used directly. Here: the example disk at 45 degrees inclination, 575 nm
# (not circularly symmetric). The unit, pixel scale and wavelength come from its FITS header; the
# plane is chosen with component (and inclination / azimuth). Use "disk" (the total minus the direct
# starlight) or "scattered_star": the "total" plane includes the star itself, which outshines the disk,
# and the coronagraph is there to remove it.
from toycoronagraph import DATADIR
with open("inclined_disk.toml", "w") as f:
    f.write('[[source]]\nkind = "mcfost"\nfile = "%s45_575"\ncomponent = "disk"\n' % DATADIR)
disk = toy.Target(para_file="inclined_disk.toml")

'''
Input: 45_575.fits (MCFOST disk, inclination 0, azimuth 0), 455 x 455, W/m^2/pixel; 0.006044 arcsec/pixel (header CDELT2 (deg) of 45_575.fits) = 55.65 lambda/D across at 5.75e-07 m (header WAVE (um) of 45_575.fits)
'''

# The target is not circularly symmetric, so "auto" uses the spatial method
method, spatial_mode = disk.check_symmetry()

'''
Method: spatial, rotate mode (auto; target not symmetric: asymmetry 70%, light in the corners 0% (tol 1%); PSFs symmetric: vortex model)
'''

# Images of other targets start with the target's name: 45_575_origin.png and 45_575_charge6_spatial_final.png
disk.plot_origin(log=True)
final_image = disk.plot_final(charge=6)

# -------------------- Your own target: any FITS image --------------------

# Any image works: a square 2-D array with the star at its centre, and the header keywords BUNIT
# (W/m^2/pixel, Jy/arcsec^2 or Jy/pixel), PIXSCALE (arcsec per pixel) and WAVELEN (m) - or set
# input_unit / wavelength / fov in the parameter file instead. Here: an inclined ring.
import numpy as np
from astropy.io import fits
n = 256
y, x = np.indices((n, n))-(n-1)/2
ring = np.exp(-0.5*((np.hypot(x, y/np.cos(np.radians(60)))-60)/8)**2)  # a ring inclined by 60 degrees
header = fits.Header([("BUNIT", "Jy/arcsec^2"), ("PIXSCALE", 0.01), ("WAVELEN", 0.8e-6)])
fits.PrimaryHDU(ring, header).writeto("my_ring.fits", overwrite=True)
with open("my_ring.toml", "w") as f:
    f.write('[observation]\ndiameter = 2.4\n[[source]]\nkind = "fits_image"\nfile = "my_ring"\n[method]\nerror_target = 0.01\n')
my_target = toy.Target(para_file="my_ring.toml")
method, spatial_mode = my_target.check_symmetry()
final_image = my_target.plot_final(charge=2)  # saved to my_ring_charge2_spatial_final.png

'''
Input: my_ring.fits (image), 256 x 256, Jy/arcsec^2; 0.01 arcsec/pixel (header PIXSCALE of my_ring.fits) = 37.23 lambda/D across at 8e-07 m (header WAVELEN of my_ring.fits)
Method: spatial, rotate mode (auto; target not symmetric: asymmetry 94%, ...)
Kernel size K=117 (+-8.4 lambda/D): estimated crop error 0.77% of the image light (80% of the 1% error target; 48 probe pixels)
...
Estimated error of this image (flux-weighted, from the PSF check): typical 0.89%, up to 0.94%; largest part: kernel crop (K=117) (0.79%)
'''

# -------------------- Scenes: a real system built from modules --------------------

import numpy as np
from toycoronagraph.params import Param, ParamError
from toycoronagraph.sources import Source, register_source
from toycoronagraph.units import SourceImage
from toycoronagraph.scene import Scene, Star, write_scene
from toycoronagraph.sources.thomson import Wind
from toycoronagraph.sources.cme import Cme
from toycoronagraph.sources.planets import Planet
from toycoronagraph.sources.images import SceneFits, FitsImage
from toycoronagraph.observe import observe, video, Result


# Your own module: a class declaring its parameters and drawing its image (template: tests/extra_sources.py)
@register_source
class Exozodi(Source):
    """A face-on exozodiacal dust disk: surface brightness ~ r^-alpha between r_in and r_out, holding
    `total` of the star's light within the image (contrast, so no stellar spectrum is needed)."""
    kind = "exozodi"
    doc = "a face-on exozodiacal dust disk"
    family = "relative"                      # a fraction of the star's light
    requires = ("star.distance",)            # au <-> arcsec
    params = (
        Param("total", 2e-7, doc="the disk's light in the image / the star's light", min=0),
        Param("r_in", 0.1, unit="AU", doc="inner edge", min=0),
        Param("r_out", 10.0, unit="AU", doc="outer edge (10 au: beyond the image, so no edge shows)", min=0),
        Param("alpha", 2.34, doc="surface brightness ~ r^-alpha (zodiacal dust: 2.34)", min=0),
    )

    def render(self, grid, star):
        """The disk's contrast per pixel on the scene grid (SourceImage, contrast_pixel)."""
        s = self.settings
        x, y = grid.coordinates()                      # arcsec from the star
        r = np.hypot(x, y)*star.distance                # au (1 arcsec at d pc = d au)
        inside = (r >= s["r_in"]) & (r <= s["r_out"])
        brightness = np.where(inside, np.maximum(r, 1e-9)**-s["alpha"], 0.0)
        return SourceImage(s["total"]*brightness/brightness.sum(), "contrast_pixel", origin="exozodi")

    def symmetric(self):
        """Face-on: circularly symmetric."""
        return True

# ---------- First, the exozodi alone, through two sets of PSFs ----------

# The exozodi alone, through two sets of PSFs: hcipy's ideal vortex and Lazuli's (prysm), both charge 6
eps_eri = Star(distance=3.22, radius=0.73)                     # pc, R_sun (temperature=5080 K: for outputs in Jy)
dust = Scene(star=eps_eri, wavelength=632e-9, diameter=3.0,    # m: Lazuli's design wavelength and aperture
             sampling=5, fov=3.0)                              # 5 px per lambda/D over the WCC's 3.0" field: 345 px
dust.add(Exozodi(total=2e-7))                                  # ~10 zodis
dust_result = observe(dust, ["vortex:6", {"type": "lazuli", "charge": 6}])
fig = dust_result.plot_compare("vortex_c6", "lazuli_c6", unit="contrast_peak", axis="lambda/D")

# ---------- The whole system ----------

# epsilon Eridani at 632 nm through Lazuli's 3.0 m aperture: the wind, a CME, the exozodi and the planet
scene = Scene(star=eps_eri, wavelength=632e-9, diameter=3.0,   # m: Lazuli's design wavelength and aperture
              sampling=5, fov=3.0)                              # 5 px per lambda/D over the WCC's 3.0" field: 345 px
scene.add(Wind(mdot=6e-13, v=400))                               # M_sun/yr (30x the Sun's), km/s
scene.add(Cme(mass=1e17))                                        # g: as heavy as the largest solar CMEs
scene.add(Exozodi(total=2e-7))                                   # ~10 zodis
scene.add(Planet(mode="polar", pos=[1.09, 30], brightness=3e-9, brightness_unit="contrast"),
          name="eps_eri_b")                                      # 3.5 au; reflected light
images = scene.render()
fig = scene.plot(log=True, unit="contrast_pixel", axis="au")
for name, image in images.items():                                # each input's light / the star's
    print("%-8s %.2e" % (name, scene.in_unit(image, "contrast_pixel").sum()))

'''
wind     4.66e-07
cme      9.23e-12
exozodi  2.00e-07
'''

# -------------------- PSF options: the coronagraphs --------------------

# vortex: an ideal vortex (hcipy) on the scene's aperture; lazuli: Lazuli's PSFs (prysm) - the telescope
# alone (no vortex: charge = 0), with an ideal vortex (charge 6), or alone over part of the r' band (weights
# from the mirrors, the filter and the detector). Without an OPD map the optics are ideal.
coronagraphs = ["vortex:2", "vortex:6",
                {"type": "lazuli", "charge": 0, "name": "lazuli_alone"},
                {"type": "lazuli", "charge": 6},
                {"type": "lazuli", "charge": 0, "wavelengths": ["600 nm", "632 nm", "665 nm"], "name": "lazuli_band"}]
result = observe(scene, coronagraphs)   # combine = "separate": each component through each coronagraph
fig = result.plot(unit="contrast_peak", components=["total"], log_range=6)   # the totals, side by side

# The input and two coronagraphs, one column per component
fig = result.plot(unit="contrast_peak", coronagraphs=["vortex_c6", "lazuli_c6"])

# An illustrative OPD map: 20 nm RMS of astigmatism and coma, in nm, on 512 x 512 pupil samples (enough for this
# field; a coarser map would set fewer pupil samples, and a SamplingWarning would say so)
n = 512
v, u = (np.indices((n, n))-(n-1)/2)/(n/2)                         # pupil coordinates, -1..1
rho, phi = np.hypot(u, v), np.arctan2(v, u)
opd = np.where(rho <= 1, rho**2*np.cos(2*phi) + (3*rho**3-2*rho)*np.sin(phi), 0.0)
opd *= 20/opd[rho <= 1].std()
fits.PrimaryHDU(opd, fits.Header([("BUNIT", "nm")])).writeto("lazuli_example_opd.fits", overwrite=True)
aberrated = observe(scene, [{"type": "lazuli", "charge": 0, "opd": "lazuli_example_opd.fits", "name": "lazuli_opd"}],
                    components=["eps_eri_b"])                    # the planet only: its PSF is computed exactly
ideal = result.final("lazuli_alone", "eps_eri_b")
print("Strehl ratio of the planet's image: %.3f" % (aberrated.final("lazuli_opd", "eps_eri_b").max()/ideal.max()))

'''
Strehl ratio of the planet's image: 0.961
'''

# -------------------- How much of the final image is each input? --------------------

# Each component's share of the image, pixel by pixel: the input and every coronagraph
fig = result.plot_contributions(axis="arcsec")

# Inside the planet's core (its final image above half its peak): the planet against everything else
for coronagraph in result.finals:
    core = result.contrast("eps_eri_b", coronagraph)
    print("%-12s planet / rest of the light in its core: %.2g" % (coronagraph, core["eps_eri_b"]/(sum(core.values())-core["eps_eri_b"])))

'''
vortex_c2    planet / rest of the light in its core: 2.1e+02
vortex_c6    planet / rest of the light in its core: 2.1e+02
lazuli_alone planet / rest of the light in its core: 2.1e+02
lazuli_c6    planet / rest of the light in its core: 2.1e+02
lazuli_band  planet / rest of the light in its core: 2.1e+02
'''

fig = result.plot_profiles(components=True, unit="contrast_peak", axis="au")
for coronagraph in result.finals:
    iwa = result.iwa(coronagraph)                                   # arcsec, None without a coronagraph
    print("%-12s IWA: %s" % (coronagraph, "none" if iwa is None else "%.3f arcsec = %.2f au" % (iwa, iwa*eps_eri.distance)))

'''
vortex_c2    IWA: 0.038 arcsec = 0.12 au
vortex_c6    IWA: 0.103 arcsec = 0.33 au
lazuli_alone IWA: none
lazuli_c6    IWA: 0.103 arcsec = 0.33 au
lazuli_band  IWA: none
'''

# ---------- Choosing what a plot shows ----------

# The same settings in every plot; other keywords go to matplotlib
fig = result.plot(unit="contrast_peak", components=["total"], log_range=4, cmap="magma")

# A figure of your own, from the data behind the plots
import matplotlib.pyplot as plt
from toycoronagraph.scene import radial_profile
fig, ax = plt.subplots(figsize=(6, 4))
for coronagraph in ["vortex_c2", "vortex_c6"]:
    r, mean = radial_profile(result.final(coronagraph, "wind", unit="contrast_peak"), result.grid)
    line, = ax.loglog(r[1:]*eps_eri.distance, mean[1:], label=coronagraph)   # arcsec -> au; [1:]: not the star pixel
    ax.axvline(result.iwa(coronagraph)*eps_eri.distance, ls=":", color=line.get_color())
ax.set(xlabel="distance from the star [au]", ylabel="azimuthal mean [contrast_peak]",
       title="The wind through two vortices (dotted: their IWA)")
ax.legend()
plt.show()

# -------------------- A time series: a CME moving outwards --------------------

# The CME moving outwards, through the charge-6 vortex: the frames minus the wind, the exozodi and the planet
frames = [{"horizontal_height_factor": h} for h in np.linspace(0, 3, 20)]
path, cme_frames = video(scene, "vortex:6", "cme", frames, difference=True, log_range=3)
print("Saved the video to", path)

# -------------------- Save, load, and scene files --------------------

# A run folder: FITS images (input and final, one extension per component), the parameters used, quick looks
folder = result.save("results/eps_eri", overwrite=True)
again = Result.load(folder)                                    # re-plot or change units without computing again

# The scene's images as one multi-extension FITS file, read back as an input module (scene_fits)
write_scene("eps_eri_scene.fits", images, scene.internal_unit, scene.grid, eps_eri, overwrite=True)
from_file = Scene(star=eps_eri, wavelength=632e-9, diameter=3.0, sampling=5, fov=3.0)
from_file.add(SceneFits(file="eps_eri_scene.fits"))
print("scene_fits round trip identical:", all(np.allclose(from_file.render(verbose=False)[n.lower()], images[n]) for n in images))

'''
scene_fits round trip identical: True
'''

# The same kind of scene from a parameter file (built-in modules only; your own need module = "...")
with open("eps_eri.toml", "w") as f:
    f.write("""[observation]
wavelength = "632 nm"
diameter = 3.0
[grid]
size = 163
pixel_scale = 0.018445   # arcsec: Lazuli's WCC detector pixel
[star]
distance = 3.22
radius = 0.73
[[source]]
kind = "wind"
mdot = 6e-13
v = 400
[[source]]
kind = "planet"
name = "eps_eri_b"
mode = "moving"
pos = [1.09, 0.07, 30, 79, 0.0]   # a (arcsec), e, PA (deg), inclination (deg), time (periods)
brightness = 3e-9
brightness_unit = "contrast"
[[coronagraph]]
type = "lazuli"
charge = 6
""")
file_result = observe(Scene.from_file("eps_eri.toml"))   # the file's own [[coronagraph]] tables
fig = file_result.plot(unit="contrast_peak")

# Images from files work in scenes too: the MCFOST disk of the section above (its own grid and units)
disk_result = observe(Scene.from_file("inclined_disk.toml"), ["vortex:6"])
fig = disk_result.plot()

# -------------------- The strict grid --------------------

# [grid] resample = false: images must cover the grid's field; the scene file written above does
strict = Scene(star=eps_eri, wavelength=632e-9, diameter=3.0, sampling=5, fov=3.0, resample=False)
strict.add(SceneFits(file="eps_eri_scene.fits"))
strict.add(Planet(mode="polar", pos=[1.09, 30], brightness=3e-9, brightness_unit="contrast"), name="eps_eri_b")
images_strict = strict.render()

# An image of another field (2.0" across instead of 3.0") stops the run, and the message says how to fix it
n = 200
y, x = (np.indices((n, n))-(n-1)/2)*0.01                          # arcsec
ring = np.exp(-0.5*((np.hypot(x, y)-0.6)/0.05)**2)
header = fits.Header([("BUNIT", "Jy/arcsec^2"), ("PIXSCALE", 0.01), ("WAVELEN", 632e-9)])
fits.PrimaryHDU(ring, header).writeto("ring_2arcsec.fits", overwrite=True)
strict.add(FitsImage(file="ring_2arcsec.fits"))
try:
    strict.render()
except ParamError as error:
    print("Problems:", error)

'''
Scene: 345 x 345 px, 0.0086906 arcsec/pixel (5 px per lambda/D at 6.32e-07 m); wind, cme, exozodi, eps_eri_b; in contrast_pixel
Problems: [grid] resample = false, but these images cover another field than the grid:
  ring_2arcsec.fits (image): 200 px of 0.01 arcsec (2 arcsec across), the grid is 345 px of 0.00869062 arcsec (2.99827 across)
  give the grid their field (for another resolution set only [grid] size, with pixel_scale = "auto"), or set [grid] resample = true
'''
