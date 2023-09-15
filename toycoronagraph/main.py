import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
# plt.rcParams['text.usetex'] = True
import hcipy
from toycoronagraph.psf import psf_calculation, cir_psf, cir_psf_planets, cir_psf_contrast
from toycoronagraph.planet import planet_position, orbit_position, orbit_plot
from toycoronagraph.tool import convert_to_polar, is_positive_even_integer, is_planet_pos_allowed
from toycoronagraph.para import example_para
import importlib.util
from astropy.io import fits
from toycoronagraph import DATADIR
import os
import importlib
# Find the parameters of the toycoronagraph at toycoronagraph_para.py
para_spec = importlib.util.find_spec("toycoronagraph_para")
if para_spec is None:
    # If the module doesn't exist, provide an example parameter setup
    example_para()
    print("There is no toycoronagraph_para.py file, create an example one")
import toycoronagraph_para as par
import cv2

class Target(object):
    """A circular symmetric target
    
    A circular symmetric target, which could be a dust ring, a debris disk, or ice remains.

    Attributes:
        file_name (str): The name of the FITS file containing the target data.
        taregt_pixel (int): Number of pixels along the edge of the target image.
        psf_scale (float): The scale of the PSF in arcseconds per pixel.
        xpix (np.ndarray): The x-coordinates of the pixels in arcseconds.
        ypix (np.ndarray): The y-coordinates of the pixels in arcseconds.
        pre_img (np.ndarray): The pre-processed image.
        planets (list): A list of the planets in the system.
        orbits (list): A list of the orbits of the planets.
        planets_brightness (list): A list of the brightnesses of the planets in Jy/arcsec^2.

    Methods:
        __init__(): Initialize the Target object.
    """
    def __init__(self):
        """ Initialization
        
        Initializes the instance based on the target fits file.

        Args:
            None
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
        origin = np.reshape(fits_read[0].data[5,0],(par.taregt_pixel,par.taregt_pixel)) # vF_v(W/m^2/pixel)
        F_v = origin*par.F_transfer # F_v(Jy/arcsec^2)
        
        self.xpix = (np.arange (-par.taregt_pixel/2, par.taregt_pixel/2, 1))*par.psf_scale
        self.ypix = (np.arange (-par.taregt_pixel/2, par.taregt_pixel/2, 1))*par.psf_scale
        
        self.pre_img = F_v.astype(np.float64)
        self.planets = []
        self.orbits = []
        self.planets_brightness = []
        
    def plot_origin(self, plot_planets=True, res=1000, plot_dpi=300, boundary=True, flip=True):
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
        x_range = par.taregt_pixel*par.psf_scale/2
        y_range = par.taregt_pixel*par.psf_scale/2
        
        if boundary and flip:
            # Set axis limits
            axs.set_xlim(-x_range, x_range)
            axs.set_ylim(y_range, -y_range)
        elif flip:
            axs.invert_yaxis()
        elif boundary:
            axs.set_xlim(-x_range, x_range)
            axs.set_ylim(-y_range, y_range)

        # Set axis labels
        axs.set_ylabel('y [arcsec]')
        axs.set_xlabel('x [arcsec]')

        # Add color bar
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
                        axs.annotate(str(b*par.psf_scale**2), (x,y), color=color)  
                else:
                    circle = plt.Circle((x,y), radius=circle_size, color=color)
                    axs.add_artist(circle)
                    axs.annotate(str(b*par.psf_scale**2), (x,y), color=color)
                    
                if orbit is not None:
                    orbit_x, orbit_y = orbit_position(orbit[0], orbit[1], orbit[2], orbit[3], res)
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
            brightness (float): Planet brightness in Jy.
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
                self.planets.append(planet_position(pos[0],pos[1],pos[2],pos[3],pos[4],mode="polar",res=1000))
                self.planets_brightness.append(brightness/par.psf_scale**2)
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
                self.planets_brightness.append(brightness/par.psf_scale**2)
                self.orbits.append(None)
            
        elif mode == "cartesian":
            # Check if planet position and brightness are valid
            if not is_planet_pos_allowed(pos, mode):
                print("planet position is not allowed")
            elif not isinstance(brightness, (float, int)) or brightness<=0:
                print("Brightness is invalid")
            else:
                self.planets.append(convert_to_polar(pos))
                self.planets_brightness.append(brightness/par.psf_scale**2)
                self.orbits.append(None)
        else:
            print("no "+str(mode)+" mode")

    def list_planets(self):
        """List planets
        
        List the planets present in the Target object with their order numbers.

        Returns:
            None
        """
        # Check if the planet list is empty
        if len(self.planets)==0:
            print("There is no planet, but you can add planet via Target.add_planet(pos, brightness)")
        else:
            order = 1
            for p, b, o in zip(self.planets, self.planets_brightness, self.orbits):
                if o is None:
                    print("Static Planet {:d}: ({}, {}) arcsec, brightness: {:.2e} Jy".format(order, p[0]*np.cos(p[1]), p[0]*np.sin(p[1]), b*par.psf_scale**2))
                else:
                    print("Moving Planet {:d}: ({}, {}) arcsec, brightness: {:.2e} Jy".format(order, p[0]*np.cos(p[1]), p[0]*np.sin(p[1]), b*par.psf_scale**2))
                order += 1

    def delete_planet(self, order=1):
        """Delete planet
        
        Deletes a planet from self.planet, self.orbit, and self.planet_brightness with a given order number.

        Args:
            order (int): The order number of the planet to be deleted. Defaults to 1.

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

    def planet_move(self, time, order=1, mode="cumulate", plot_pos=False, plot_dpi=300, res=1000, flip=True, message=True):
        """Planet moving
        
        Move a planet's position.

        Args:
            time (float): The time by which the planet should move.
            order (int): The order number of the planet.
            mode (str): Mode for specifying movement ("culmulated", "specific"). The former will be added to the original time, while the latter will set up a new time regardless of the original time.
            plot_pos (bool): Whether to plot the new planet's position in its orbit.
            plot_dpi (int): Dots per inch (DPI) for the plot.

        Returns:
            None
        """
        # Check if the input order is valid
        if not isinstance(order, int) or order < 1 or order > len(self.planets):
            print("Input order number is not existed, please try again or use Target.list_planets() to list all planets")
        elif len(self.orbits[order-1]) != 5:
            print("A static planet is unable to move")
        else:
            if mode == "cumulate":
                # Update the time of the planet's orbit by adding to the current time
                self.orbits[order-1][4] += time
            elif mode == "specific":
                # Set a new time for the planet's orbit
                self.orbits[order-1][4] = time
    
            # Get the updated position of the planet based on the updated time
            pos = self.orbits[order-1]
            self.planets[order-1] = planet_position(pos[0],pos[1],pos[2],pos[3],pos[4],mode="polar",res=1000)
            if plot_pos:
                # Plot the new planet position on its orbit
                orbit_plot(pos[0],pos[1],pos[2],pos[3], self.planets[order-1], plot_dpi, res, flip, "polar", '_planet'+str(order))
            elif message:
                print("Planet has moved to new position")

    def planet_video(self, charge, order=1, iwa_ignore=False, plot_dpi=300, length=3, fps=20, flip=True):
        """Planet Video
           
        Create a movie of the planet moving in its orbit.
        
        Args:
            charge (int): Charge number (for vortex coronagraph).
            order (int): The order number of the planet, using self.list_planets() to look up the order.
            iwa_ignore (bool): Whether to ignore the dust inside IWA.
            plot_dpi (int): Dots per inch (DPI) for the plot.
            length (int): The length of the movie in seconds.
            fps (int): The frames per second of the movie.
            flip (bool): Whether to flip the x-axis of the plot.

        Returns:
            planet_video_*.mp4        
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
                psf_calculation(charge, par.taregt_pixel, par.psf_range, par.lyot_mask_size, par.num_cores)
                
        # Check if the input order is valid
        if not isinstance(order, int) or order < 1 or order > len(self.planets):
            print("Input order number is not existed, please try again or use Target.list_planets() to list all planets")
        elif len(self.orbits[order-1]) != 5:
            print("A static planet is unable to move")
        else:
            total_frames = length*fps
            target_img = cir_psf(self.pre_img, [], [], par.psf_scale, iwa_ignore, False, par.taregt_pixel, par.rot_number, psf_filename)
            fig,[ax,cax] = plt.subplots(1,2, width_ratios=[50,1], figsize=(10, 8))
            # Set the colormap and norm to correspond to the data for which
            # the colorbar will be used.
            cmap = matplotlib.cm.gnuplot #matplotlib.cm.winter
            planet_brightness = self.planets_brightness[order-1]
            t0 = self.orbits[order-1][-1]
            norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(target_img)+planet_brightness*np.max(np.load(psf_filename)))
            colorbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,norm=norm, orientation='vertical')
            colorbar.set_label(r"Jy/arcsec^2") 
            art = ax.scatter([],[],c=[])
            initial_img = cir_psf_planets([self.planets[order-1]], [planet_brightness], par.psf_scale, par.taregt_pixel, psf_filename)
            img = ax.imshow(target_img+initial_img, extent=[np.min(self.ypix), np.max(self.ypix), np.min(self.xpix), np.max(self.xpix)], cmap=cmap, norm=norm)
            ax.set_ylabel('y [arcsec]')
            ax.set_xlabel('x [arcsec]')
            def animate(num):
                if num % fps == 0:
                    print("{:.0%}".format(num/total_frames))
                self.planet_move(t0+num/total_frames, order=order, mode="specific", res=total_frames*10, flip=True, message=False)
                planet_pos = self.planets[order-1]
                currrent_img = target_img+cir_psf_planets([planet_pos], [planet_brightness], par.psf_scale, par.taregt_pixel, psf_filename)
                img.set_array(currrent_img)
                #art.set_color(cmap(norm(currrent_img)))
                return [img]
            anim= animation.FuncAnimation(fig, animate, interval=1000/fps, frames=total_frames)
            anim_name = 'planet_video_iwa_ignore.mp4' if iwa_ignore else 'planet_video.mp4'
            anim.save(anim_name, fps=fps, extra_args=['-vcodec', 'libx264'])
            
    def plot_orbit(self, order=1, plot_dpi=300, res=1000, flip=True):
        """Orbit plot
        
        Plot the orbit of a planet.

        Args:
            order (int): The order number of the planet, using self.list_planets() to look up the order.
            plot_dpi (int): Dots per inch (DPI) for the plot.

        Returns:
            orbit_*.png
        """
        pos = self.orbits[order-1]
        orbit_plot(pos[0],pos[1],pos[2],pos[3], self.planets[order-1], plot_dpi, res, flip, "polar", '_planet'+str(order))

    def contrast(self, charge=2, order=1, plot_contrast=True, plot_dpi=300, exozodi_limit=10, res=100):
        """Contrast
        
        Find the planet brightness and background brightness in the final image.

        Args:
            charge (int): Charge number (for vortex coronagraph).
            order (int): The order number of the planet, using self.list_planets() to look up the order.
            plot_contrast (bool): Whether to plot the dust to planet brightness ratio.
            plot_dpi (int): Dots per inch (DPI) for the plot.
            exozodi_limit (float): The maximum exozodi to target dust ratio.
            res (int): The number of points in the exozodi ratio range.

        Returns:
            plot or print planet brightness, background brightness, and background brightness (ignored dust inside IWA) in Jy.
        """
        # Check if the input order is valid
        if not isinstance(order, int) or order < 1 or order > len(self.planets):
            print("Input order number is not existed, please try again or use Target.list_planets() to list all planets")
        else:
            # Check if the planet inside the plot
            planet_pos = self.planets[order-1]
            planet_psfs_number = int(planet_pos[0]/par.psf_scale)
            if planet_psfs_number>=par.taregt_pixel/2:
                print("Planet is outside the range of the plot")
            else:
                # Check and adjust the charge number if using a vortex coronagraph
                if par.coronagraph_type=='vortex':
                    if not is_positive_even_integer(charge):
                        print("Charge number is not compatible with vortex coronagraph type, auto set to 2")
                        charge = 2
                    print("This is a charge-" + str(charge) + " vortex coronagraph.")
                        
                # Define the path to the PSF file based on charge number
                psf_filename = DATADIR+"psfs_c"+str(charge)+".npy"
        
                # If the PSF file doesn't exist, calculate it
                if not os.path.exists(psf_filename):
                    psf_filename = "psfs_c"+str(charge)+".npy"
                    if not os.path.exists(psf_filename):
                        psf_calculation(charge, par.taregt_pixel, par.psf_range, par.lyot_mask_size, par.num_cores)

                planet_b, dust_b, dust_b_iwa = cir_psf_contrast(self.pre_img, planet_psfs_number, planet_pos[1], self.planets_brightness[order-1], par.psf_scale, par.taregt_pixel, par.rot_number, psf_filename)
                if plot_contrast:
                    # Create a new figure and axes for plotting
                    fig = plt.figure(dpi=plot_dpi)
                    axs = plt.gca()                       
                    axs.set_ylabel('dust to planet brightness ratio')
                    axs.set_xlabel('exozodi to target dust ratio')
                    exozodi_ratio = np.linspace(0, exozodi_limit, res)
                    dp_ratio = (dust_b / planet_b) * exozodi_ratio
                    dp_ratio_iwa = (dust_b_iwa / planet_b) * exozodi_ratio
                    plt.plot(exozodi_ratio, dp_ratio, label="Complete")
                    plt.plot(exozodi_ratio, dp_ratio_iwa, label="IWA ignored")
                    plt.legend()
                    plt.show()
                else:
                    print (planet_b, dust_b, dust_b_iwa)
    
    def plot_final(self, charge, iwa_ignore=False, add_planet=True, plot_dpi=300, flip=True):
        """Final image
        
        Plots the target image after processing with the vortex coronagraph.

        Args:
            charge (int): Charge number (for vortex coronagraph).
            add_planet (bool): Whether to include planets in the plot.
            rot_number (int): Rotation number for creating the circular PSF.
            plot_dpi (int): Dots per inch (DPI) for the final plot.

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
                psf_calculation(charge, par.taregt_pixel, par.psf_range, par.lyot_mask_size, par.num_cores)
        
        # Generate the final image of the disk using the cir_psf function
        final_img = cir_psf(self.pre_img, self.planets, self.planets_brightness, par.psf_scale, iwa_ignore, add_planet, par.taregt_pixel, par.rot_number, psf_filename)
        
        # Create a new figure and axes for plotting
        fig = plt.figure(dpi=plot_dpi)
        axs = plt.gca()

        # Plot the final image
        img = axs.imshow(final_img, cmap='gnuplot', 
                         extent=[np.min(self.ypix), np.max(self.ypix), np.min(self.xpix), np.max(self.xpix)])
        if flip:
            axs.invert_yaxis()
            
        axs.set_ylabel('y [arcsec]')
        axs.set_xlabel('x [arcsec]')

        # Add color bar
        colorbar=plt.colorbar(img,orientation='vertical')
        colorbar.set_label(r"Jy/arcsec^2")
        
        # Create a name for the final image
        final_image_name = "charge"+str(charge)
        if add_planet and self.planets != []:
            final_image_name += "_with_planets"

        if iwa_ignore:
            final_image_name += "_iwa_ignore"
        
        # Save the final image
        fig.savefig(final_image_name+"_final.png", format='png', bbox_inches='tight')

        # Show the plot
        plt.show()