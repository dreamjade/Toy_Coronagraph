import numpy as np
import matplotlib.pyplot as plt
# plt.rcParams['text.usetex'] = True
import hcipy
from toycoronagraph.psf import psf_calculation, cir_psf
from toycoronagraph.planet import planet_position, orbit_position, orbit_plot
from toycoronagraph.tool import convert_to_polar, is_positive_even_integer, is_planet_pos_allowed
from toycoronagraph.para import example_para
import importlib.util
from astropy.io import fits
from toycoronagraph import DATADIR
import os

# Find the parameters of the toycoronagraph at toycoronagraph_para.py
para_spec = importlib.util.find_spec("toycoronagraph_para")
if para_spec is None:
    # If the module doesn't exist, provide an example parameter setup
    example_para()
    print("There is no toycoronagraph_para.py file, create an example one")
import toycoronagraph_para as par
import importlib

class Target(object):
    """A circular symmetric target
    
    A circular symmetric target, which could be a dust ring, a debris disk, or ice remains

    Attributes:
        file_name (str): The name of the FITS file containing the target data.
        px (int): Number of pixels along the x-axis.
        py (int): Number of pixels along the y-axis.
        psf_scale (float): The scale of the PSF in arcseconds per pixel.
        xpix (np.ndarray): The x-coordinates of the pixels in arcseconds.
        ypix (np.ndarray): The y-coordinates of the pixels in arcseconds.
        pre_img (np.ndarray): The pre-processed image.
        planets (list): A list of the planets in the system.
        orbits (list): A list of the orbits of the planets.
        planets_brightness (list): A list of the brightnesses of the planets.       

    Methods:
        __init__(): Initialize the Target object.
    """
    def __init__(self):
        """ Initialization
        
        Initializes the instance based on the target fits file

        Args:
            file_name (str, optional): The name of the FITS file containing the target data. 
            If None, the default file `DATADIR/example.fit` will be used.
            px (int, optional): The number of pixels in the x-direction. Defaults to 512.
            py (int, optional): The number of pixels in the y-direction. Defaults to 512.
        """
        # Reload the module 'par' to get the latest parameters
        importlib.reload(par)

        # Load the fits file
        if par.fits_filename == None:
            self.file_name = DATADIR+"example"
        else:
            self.file_name = par.fits_filename
        fits_read = fits.open(self.file_name+".fits")
        
        # Extract data from the FITS file and perform necessary calculations
        origin = np.reshape(fits_read[0].data[5,0],(par.px,par.py)) # vF_v(W/m^2/pixel)
        F_v = origin*par.F_transfer # F_v(Jy/arcsec^2)
        
        self.xpix = (np.arange (-par.px/2, par.px/2, 1))*par.psf_scale
        self.ypix = (np.arange (-par.py/2, par.py/2, 1))*par.psf_scale
        
        self.pre_img = F_v.astype(np.float64)
        self.planets = []
        self.orbits = []
        self.planets_brightness = []
        
    def plot_origin(self, plot_planets=True, plot_dpi=300, boundary=True):
        """Plot original target
        
        Plot the original image with or without planets.

        Args:
            plot_planets (bool): Whether to plot planets or not.
            plot_dpi (int): Dots per inch (DPI) for the plot.
            boundary (bool): Whether to show the boundary limits.

        Returns:
            origin_*.png
        """
        # Create a new figure and axes for plotting
        fig = plt.figure(dpi=plot_dpi)
        axs = plt.gca()

        # Plot the pre-processed image data
        img = axs.imshow(self.pre_img,
                       cmap='gnuplot',extent=[np.min(self.ypix),np.max(self.ypix),np.min(self.xpix),np.max(self.xpix)])

        # Calculate boundary limits based on pixel scale
        x_range = par.px*par.psf_scale/2
        y_range = par.py*par.psf_scale/2
        
        if boundary:
            # Set axis limits
            axs.set_xlim(-x_range, x_range)
            axs.set_ylim(y_range, -y_range)
        else:
            axs.invert_yaxis()

        # Set axis labels
        axs.set_ylabel('y [arcsec]')
        axs.set_xlabel('x [arcsec]')

        # Add colorbar
        colorbar=plt.colorbar(img,orientation='vertical')
        colorbar.set_label(r"Jy/arcsec^2")

        # Iterate through planets, their brightness, orbits, and orders to add them into the plot
        if plot_planets and self.planets !=[]:
            circle_size = 0.02*min(x_range, y_range)
            colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#bcbd22', u'#17becf']
            for p, b, orbit, order in zip(self.planets, self.planets_brightness, self.orbits, range(len(self.planets))):
                x,y = p[0]*np.cos(p[1]),p[0]*np.sin(p[1])
                color = colors[order%9]
                
                if abs(x) > x_range or abs(y) > y_range:
                    print("Planet #"+str(order+1)+" is outside the range of the plot")
                    color = u'#7f7f7f'
                    
                    if not boundary:
                        circle = plt.Circle((x,y), radius=circle_size, color=color)
                        axs.add_artist(circle)
                        axs.annotate(str(b), (x,y), color=color)  
                else:
                    circle = plt.Circle((x,y), radius=circle_size, color=color)
                    axs.add_artist(circle)
                    axs.annotate(str(b), (x,y), color=color)
                    
                if orbit is not None:
                    orbit_x, orbit_y = orbit_position(orbit[0], orbit[1], orbit[2], orbit[3])
                    plt.plot(orbit_x, orbit_y, color=color, linestyle = '--', alpha = 0.5) # Plot the orbit.
                    
            # Save the figure
            fig.savefig("origin_with_planets.png", format='png', bbox_inches='tight')
        else:
            fig.savefig("origin.png", format='png', bbox_inches='tight')

        # Show the plot
        plt.show()
        
    def add_planet(self, pos, brightness, mode="moving"):
        """Add a planet
        
        Add a planet to the Target object.

        Args:
            pos (array-like): Planet position information (depends on mode). For "moving" mode, pos = [a, e, pa, inc, t]; for "polar" mode, pos =[r, theta(in degree)]; for "cartesian" mode, pos =[x,y].
            brightness (float): Planet brightness.
            mode (str): Mode for specifying planet position coordinates ("moving", "polar", "cartesian").

        Returns:
            None
        """
        pos = np.array(pos)

        if mode == "moving":
            # Check if planet position and brightness are valid
            if not is_planet_pos_allowed(pos, mode):
                print("Planet position is not allowed")
            elif not isinstance(brightness, (float, int)) or brightness<=0:
                print("Brightness is invalid")
            else:
                self.planets.append(planet_position(pos[0],pos[1],pos[2],pos[3],pos[4],mode="polar"))
                self.planets_brightness.append(brightness)
                self.orbits.append(pos) #pos = np.array(a, e, pa, inc, t)
                
        elif mode == "polar":
            # Check if planet position and brightness are valid
            if not is_planet_pos_allowed(pos, mode):
                print("Planet position is not allowed")
            elif not isinstance(brightness, (float, int)) or brightness<=0:
                print("Brightness is invalid")
            else:
                # Convert polar mode angle from degree to radian
                pos[1] = pos[1]*np.pi/180.0
                self.planets.append(pos)
                self.planets_brightness.append(brightness)
                self.orbits.append(None)
            
        elif mode == "cartesian":
            # Check if planet position and brightness are valid
            if not is_planet_pos_allowed(pos, mode):
                print("planet position is not allowed")
            elif not isinstance(brightness, (float, int)) or brightness<=0:
                print("Brightness is invalid")
            else:
                self.planets.append(convert_to_polar(pos))
                self.planets_brightness.append(brightness)
                self.orbits.append(None)
        else:
            print("no "+str(mode)+" mode")

    def list_planets(self):
        """List planets
        
        List the planets present in the Target object with their order numbers.

        Returns:
            None
        """
        # Check if the input order is valid
        if len(self.planets)==0:
            print("There is no planet, but you can add planet via Target.add_planet(pos, brightness)")
        order = 1
        for p, b, o in zip(self.planets, self.planets_brightness, self.orbits):
            if o is None:
                print("Static Planet {:d}: ({}, {}), brightness: {:.2e}".format(order, p[0], p[1], b))
            else:
                print("Moving Planet {:d}: ({}, {}), brightness: {:.2e}".format(order, p[0], p[1], b))
            order += 1

    def delete_planet(self, order):
        """Delete planet
        
        Deletes a planet from self.planet, self.orbit and self.planet_brightness with given order number.

        Args:
            order (int): The order number of the planet to be deleted.

        Returns:
            None
        """
        # Check if the input order is valid
        if not isinstance(order, int) or order < 1 or order > len(self.planets):
            print("Input order number is not existed, please try again or use Target.list_planets() to list all planets")
        else:
            # Delete planet, brightness, and orbit information at the given order
            del self.planets[order-1]
            del self.planets_brightness[order-1]
            del self.orbits[order-1]
            print("Successfully remove planet #"+str(order)+", here is the latest planet list:")
            self.list_planets()

    def planet_move(self, time, order=1, mode="culmulated", plot_pos=False, plot_dpi=300):
        """Planet moving
        
        Move a planet's position.

        Args:
            time (float): Time by which the planet should move.
            order (int): The order number of the planet.
            mode (str): Mode for specifying movement ("culmulated", "specific"). The former will be added to the original time, while the latter will set up a new time regardless of the original time.
            plot_pos (bool): Whether to plot the new planet position on its orbit.
            plot_dpi (int): Dots per inch (DPI) for the plot.

        Returns:
            None
        """
        if mode == "culmulated":
            # Update the time of the planet's orbit by adding to the current time
            self.orbits[order-1][4] += time
        elif mode == "specific":
            # Set a new time for the planet's orbit
            self.orbits[order-1][4] = time

        # Get the updated position of the planet based on the updated time
        pos = self.orbits[order-1]
        self.planets[order-1] = planet_position(pos[0],pos[1],pos[2],pos[3],pos[4],mode="polar")
        if plot_pos:
            # Plot the new planet position on its orbit
            orbit_plot(pos[0],pos[1],pos[2],pos[3], self.planets[order-1], mode="polar", name='_planet'+str(order), plot_dpi=300)
        else:
            print("Planet has moved to new position")

    def plot_orbit(self, order, plot_dpi=300):
        """Orbit plot
        
        Plot the orbit of a planet.

        Args:
            order (int): The order number of the planet, using self.list_planets() to look up the order.
            plot_dpi (int): Dots per inch (DPI) for the plot.

        Returns:
            orbit_*.png
        """
        pos = self.orbits[order-1]
        orbit_plot(pos[0],pos[1],pos[2],pos[3], self.planets[order-1], mode="polar", name='_planet'+str(order), plot_dpi=300)
        
    def plot_final(self, charge, add_planet=True, rot_number=par.px, plot_dpi=300):
        """Final image
        
        Plots the target image after processing with the vortex coronagraph.

        Args:
            charge (int): Charge number (for vortex coronagraph).
            add_planet (bool): Whether to include planets in the plot.
            rot_number (int): Rotation number for creating the circular PSF.
            plot_dpi (int): Dots per inch (DPI) for the plot.

        Returns:
            *_final.png
        """
        # Check and adjust the charge number if using a vortex coronagraph
        if par.coronagraph_type=='vortex':
            if not is_positive_even_integer(charge):
                print("charge number is not compatible with coronagraph type, auto set to 2")
                charge = 2
                
        # Define the path to the PSF file based on charge number
        psf_filename = DATADIR+"psfs_c"+str(charge)+".npy"

        # If the PSF file doesn't exist, calculate it
        if not os.path.exists(psf_filename):
            psf_filename = "psfs_c"+str(charge)+".npy"
            if not os.path.exists(psf_filename):
                psf_calculation(charge, par.px, par.psf_range)
        
        # Generate the final image of the disk using cir_psf function
        final_img_charge = cir_psf(self.pre_img, self.planets, self.planets_brightness, par.psf_scale, add_planet, par.px, par.psf_range, rot_number, psf_filename)
        
        # Create a new figure and axes for plotting
        fig = plt.figure(dpi=plot_dpi)
        axs = plt.gca()

        # Plot the final image
        img = axs.imshow(final_img_charge, cmap='gnuplot', 
                         extent=[np.min(self.ypix), np.max(self.ypix), np.min(self.xpix), np.max(self.xpix)])
        axs.invert_yaxis()
        axs.set_ylabel('y [arcsec]')
        axs.set_xlabel('x [arcsec]')

        # Add colorbar
        colorbar=plt.colorbar(img,orientation='vertical')
        colorbar.set_label(r"Jy/arcsec^2")
        
        # Create a name for the final image
        final_image_name = "charge"+str(charge)
        if add_planet and self.planets != []:
            final_image_name += "_with_planets"
        
        # Save the final image
        fig.savefig(final_image_name+"_final.png", format='png', bbox_inches='tight')

        # Show the plot
        plt.show()