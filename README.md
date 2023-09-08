[![Version](https://img.shields.io/badge/Version-v1.5.3-red.svg?style=flat-square)](https://github.com/dreamjade/Toy_Coronagraph/blob/main/toycoronagraph/__init__.py)
<a href="https://pypi.org/project/toycoronagraph/"><img src="https://img.shields.io/pypi/v/lrgs.svg" alt="PyPI" /></a>
[![Test](https://img.shields.io/badge/Tests-v1.5.3-yellow.svg?style=flat-square)](https://github.com/dreamjade/Toy_Coronagraph/tree/main/tests)
[![Documentation Status](https://img.shields.io/badge/Docs-v1.5.2-green.svg?style=flat-square)](https://dreamjade.github.io/Toy_Coronagraph/)
<a href="./LICENSE"><img src="https://img.shields.io/cran/l/lrgs.svg" alt="MIT License" /></a>
[![DOI](https://zenodo.org/badge/665310914.svg)](https://zenodo.org/badge/latestdoi/665310914)
[![Test](https://img.shields.io/badge/Dependencies-v1.5.3-purple.svg?style=flat-square)](https://github.com/dreamjade/Toy_Coronagraph/tree/main/requirements.txt)

# Toy Coronagraph

The toycoronagraph package is implemented in Python

Language | Release | Note
---------- | -------- | ------
Python | [On PyPI](https://pypi.python.org/pypi/toycoronagraph) | though not always the most recent [version](./toycoronagraph/__init__.py)

for simulating coronagraphs. It is designed to be simple and easy to use, while still providing a powerful and flexible framework for simulating a variety of coronagraph designs. 

The package includes a number of pre-defined coronagraph designs, as well as a library of functions for creating custom coronagraphs. It also includes a number of tools for visualizing the results of simulations, such as ray tracing plots and intensity maps.

The toycoronagraph package is open source and available on GitHub. It is a valuable tool for anyone interested in learning about coronagraphs or simulating the performance of coronagraph designs.

### Here are some of the features of this package:

* Simple and easy to use
* Powerful and flexible framework
* Pre-defined coronagraph designs
* Library of functions for creating custom coronagraphs
* Tools for visualizing the results of simulations
* Open source and available on GitHub

### Here are some of the applications of this package:

* Learning about coronagraphs
* Simulating the performance of coronagraph designs
* Designing new coronagraph designs
* Testing the performance of coronagraph hardware
* Studying the physics of light scattering

### Example

example target           |  Final image (Charge=2)  |  Final image (Charge=6)
:-------------------------:|:-------------------------:|:-------------------------:
![origin](./static/origin.png) | ![charge2_final](./static/charge2_final.png) | ![charge6_final](./static/charge6_final.png)

#### with planets

Pre-image  |  Final image (Charge=6)
:-------------------------:|:-------------------------:
![Pre_image](./static/origin_with_planets.png) | ![charge6_final](./static/charge6_with_planets_final.png)

#### Turn on/off dust inside IWA

Final image (Charge=6)  |  Final image but ignored dust inside IWA (Charge=6)
:-------------------------:|:-------------------------:
![Pre_image](./static/origin_with_planets.png) | ![charge6_final](./static/charge6_with_planets_iwa_ignore_final.png)

#### Moving planet on an elliptic orbit
Frame definition
:-------------------------:
![frame](./static/frame.svg)

Plot the orbit
:-------------------------:
![frame](./static/oribit_planet1.png)


### Example usage (Python)
This [py file](./tests/test.py) illustrates how the Python package is used.
```Python
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
(1.1044442178760744e-05, 1.7104037862959554e-05, 1.6625126181290586e-05)
'''
```
More instructions could be found in [docs](https://dreamjade.github.io/Toy_Coronagraph/).
