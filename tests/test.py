import toycoronagraph.main as toy

# Load the example target image (example.fits) in Toy_Coronagraph/toycoronagraph/example_data/ folder
toy_target = toy.Target()

# Plot the preimage; the image will auto-save to origin.png
toy_target.plot_origin()

# Plot the final image through vortex coronagraph with charge = 2; the image will auto-save to charge2_final.png
toy_target.plot_final(charge=2)

# Plot the final image through vortex coronagraph with charge = 6; the image will auto-save to charge6_final.png
toy_target.plot_final(charge=6)

# Now add a static planet, where pos=(x_position, y_position)
toy_target.add_planet(pos=[0.5,0], brightness=0.00005, mode="cartesian")

# You could also add a moving planet on an elliptic orbit, where pos=(length of the semi-major axis, eccentricity, position angle, inclination angle, time/period))
toy_target.add_planet(pos=[0.5,0.6,60,0,0.1], brightness=0.00003, mode="moving")

# Plot the preimage with planets; the image will auto-save to origin_with_planets.png
toy_target.plot_origin()

# Plot the final image with planets; the image will auto-save to charge6_with_planets_final.png
toy_target.plot_final(charge=6)

# Plot the final image by turning off the target inside IWA; the image will auto-save to charge6_with_planets_iwa_ignore_final.png
toy_target.plot_final(charge=6, iwa_ignore=True)

# List the planets
toy_target.list_planets()
'''
Static Planet 1: (0.5, 0.0) arcsec, brightness: 5.00e-05 Jy
Moving Planet 2: (-0.3818058424930897, 0.08394161212137384) arcsec, brightness: 3.00e-05 Jy
'''

# Delete a specific planet
toy_target.delete_planet(order=1)
'''
Successfully remove planet #1, here is the latest planet list:
Moving Planet 1: (-0.3818058424930897, 0.08394161212137384) arcsec, brightness: 3.00e-05 Jy
'''

# Move a moving planet on its orbit by time/period=0.1
toy_target.planet_move(time=0.1, order=1)
'''
Planet has moved to new position
'''

# Plot the orbit of a specific moving planet; the image will auto-save to oribit_planet#.png
toy_target.plot_orbit(order=1)

# Plot the new preimage with planets; the image will auto-save to origin_with_planets.png
toy_target.plot_origin()

# Plot the new final image with planets; the image will auto-save to charge6_with_planets_final.png
toy_target.plot_final(charge=6)

# Print a specific planet brightness, background brightness, and background brightness (ignored dust inside IWA) in the final image through vortex coronagraph with charge = 6
toy_target.contrast(charge=6, order=1)
'''
This is a charge-6 vortex coronagraph.
'''

# Plot the core_throughput and inner working angle (in pixel)
import toycoronagraph.psf as psf
import numpy as np
psfs = np.load("psfs_c2.npy")
psf.cir_core_throughput_plot(psfs)

# Make the planet movie
toy_target.planet_video(charge=6)
# The video will be save to planet_video.mp4

# Make the planet movie again, but ignore the dust inside IWA
toy_target.planet_video(charge=6, iwa_ignore=True)
# The video will be save to planet_video_iwa_ignore.mp4