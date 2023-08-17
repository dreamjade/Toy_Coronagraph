import numpy as np
import matplotlib.pyplot as plt
import hcipy
#from .psf import psf_calculation, cir_psf
from toycoronagraph.psf import psf_calculation, cir_psf
from toycoronagraph.planet import planet_position
from toycoronagraph.tool import convert_to_polar, is_positive_even_integer, is_planet_pos_allowed
from astropy.io import fits
from toycoronagraph import DATADIR
import os

class Target(object):
    """
    A circular symmetric target, which could be a dust ring, a debris disk, or ice remains

    Args:
        file_name: The name of the FITS file containing the target image.

    Attributes:
        px: The number of pixels in the x-direction.
        py: The number of pixels in the y-direction.
        psf_scale: The scale of the PSF in arcseconds per pixel.
        xpix: The x-coordinates of the pixels in arcseconds.
        ypix: The y-coordinates of the pixels in arcseconds.
        data_jy: The target image in units of Jy.
        pre_img: The target image in units of float64.
    """
    def __init__(self, file_name=None, px=512, py=512):
        """
        Constructor for the Target class.
        """
        if file_name == None:
            self.file_name = DATADIR+"example"
        else:
            self.file_name = file_name
        inc = fits.open(self.file_name+".fits")
        self.px = px
        self.py = py
        self.psf_scale = 1e-6/2.4*206264.806*32/512 ##arcsecs/pixel
        self.xpix = (np.arange (-self.px/2, self.px/2, 1))*self.psf_scale
        self.ypix = (np.arange (-self.px/2, self.px/2, 1))*self.psf_scale
        c=2.99792*10**14
        wave_length = 1.0 #in microns#
        jy=10**26
        sst=np.reshape(inc[0].data[5,0],(self.px,self.py))
        self.data_jy = (sst/c)*(wave_length**2)*jy
        self.pre_img = self.data_jy.astype(np.float64)
        self.planets = []
        self.orbits = []
        self.planets_brightness = []
        
    def plot_origin(self):
        """
        Plots the original target image.
        """
        fig=plt.figure(dpi=300)
        ax=plt.subplot(111)
        im=ax.imshow(self.pre_img,
                       cmap='gnuplot',extent=[np.min(self.ypix),np.max(self.ypix),np.min(self.xpix),np.max(self.xpix)])
        ax.invert_yaxis()
        ax.set_ylabel('y [arcsec]')
        ax.set_xlabel('x [arcsec]')
        cb=plt.colorbar(im,orientation='vertical')
        cb.set_label("$Jy$")
        fig.savefig("origin.png", format='png', bbox_inches='tight')
        plt.show()
        
    def plot_planets(self):
        """
        Plots all planets location on the target image.

        Args:
            target: The target object.
        """
        fig=plt.figure(dpi=300)
        ax=plt.subplot(111)
        im=ax.imshow(self.pre_img,
                       cmap='gnuplot',extent=[np.min(self.ypix),np.max(self.ypix),np.min(self.xpix),np.max(self.xpix)])
        ax.invert_yaxis()
        ax.set_ylabel('y [arcsec]')
        ax.set_xlabel('x [arcsec]')
        cb=plt.colorbar(im,orientation='vertical')
        cb.set_label("$Jy$")
        for p, b in zip(self.planets, self.planets_brightness):
            circle = plt.Circle((p[0], p[1]), 0.2, color='red', fill=False)
            ax.add_artist(circle)
            ax.annotate(str(b), (p[0], p[1]), color='red')
        fig.savefig("planets.png", format='png', bbox_inches='tight')
        plt.show()
        
    def add_planet(self, pos, brightness, mode="moving"):
        pos = np.array(pos)

        if mode == "moving":
            if not is_planet_pos_allowed(pos, mode):
                print("Planet position is not allowed")
            elif not isinstance(brightness, (float, int)) or brightness<=0:
                print("Brightness is invalid")
            else:
                self.planets.append(planet_position(pos[0],pos[1],pos[2],pos[3],pos[4],mode="polar"))
                self.planets_brightness.append(brightness)
                self.orbits.append(pos) #pos = np.array(a, e, pa, inc, t)
                
        elif mode == "polar":
            if not is_planet_pos_allowed(pos, mode):
                print("Planet position is not allowed")
            elif not isinstance(brightness, (float, int)) or brightness<=0:
                print("Brightness is invalid")
            else:
                self.planets.append(pos)
                self.planets_brightness.append(brightness)
                self.orbits.append(None)
            
        elif mode == "cartesian":
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
        """
        Lists all planets in order and show the order number.

        Args:
            target: The target object.
        """
        if len(self.planets)==0:
            print("There is no planet, but you can add planet via Target.add_planet(pos, brightness)")
        order = 1
        for p, b, o in zip(self.planets, self.planets_brightness, self.orbits):
            if o==None:
                print("Static Planet {:d}: ({}, {}), brightness: {:.2e}".format(order, p[0], p[1], b))
            else:
                print("Moving Planet {:d}: ({}, {}), brightness: {:.2e}".format(order, p[0], p[1], b))
            order += 1

    def delete_planet(self, order):
        """
        Deletes the planets from self.planet and self.planet_brightness with given order number.

        Args:
            target: The target object.
            order: The order number of the planet to be deleted.
        """
        if not isinstance(order, int) or order < 1 or order > len(self.planets):
            print("Input order number is not existed, please try again or use Target.list_planets() to list all planets")
        else:
            del self.planets[order-1]
            del self.planets_brightness[order-1]
            del self.orbits[order-1]

    def planet_move(self, time, order==1, mode=="culmulated", plot_pos==False):
        if mode == "culmulated":
            self.orbits[order-1][4] += time
        elif mode == "specific":
            self.orbits[order-1][4] = time
        pos = self.orbits[order-1]
        self.planets[order-1] = planet_position(pos[0],pos[1],pos[2],pos[3],pos[4],mode="polar")
        if plot_pos:
            orbit_plot(pos[0],pos[1],pos[2],pos[3], self.planets[order-1], mode="polar", name='_planet'+str(order))
        else:
            print("Planet has moved to new position")

    def plot_orbit(self, order):
        pos = self.orbits[order-1]
        orbit_plot(pos[0],pos[1],pos[2],pos[3], self.planets[order-1], mode="polar", name='_planet'+str(order))
        
    def plot_final(self, charge, coronagraph_type='vortex', add_planet=True, img_pixel=512, psf_range=16, rot_number=360, plot_dpi=300):
        """
        Plots the target image after processing with the vortex coronagraph.

        Args:
            charge: The vortex charge.
            img_pixel: The number of pixels in the image.
            psf_range: The range of the PSF in pixels.
            rot_number: The number of rotations.
            plot_dpi: The DPI of the plot.

        Returns:
            The final image.
        """
        #check charge number
        if coronagraph_type=='vortex':
            if not is_positive_even_integer(charge):
                print("charge number is not compatible with coronagraph type, auto set to 2")
                charge = 2
                
        #find psf files
        psf_filename = DATADIR+"psfs_c"+str(charge)+".npy"
        if not os.path.exists(psf_filename):
            psf_filename = "psfs_c"+str(charge)+".npy"
            if not os.path.exists(psf_filename):
                psf_calculation(charge, img_pixel, psf_range)
        
        #draw final image of the disk
        final_img_charge = cir_psf(self.pre_img, self.planets, self.planets_brightness, add_planet, img_pixel, psf_range, img_pixel, psf_filename)
        
        #show the final results
        fig=plt.figure(dpi=plot_dpi)
        ax2=plt.subplot(111)
        im2=ax2.imshow(final_img_charge,
                   cmap='gnuplot',extent=[np.min(self.ypix),np.max(self.ypix),np.min(self.xpix),np.max(self.xpix)])
        ax2.invert_yaxis()
        ax2.set_ylabel('y [arcsec]')
        ax2.set_xlabel('x [arcsec]')
        cb=plt.colorbar(im2,orientation='vertical')
        cb.set_label("$Jy$")
        
        #save the final image
        final_image_name = "charge"+str(charge)
        if add_planet and self.planets != []:
            final_image_name += "_with_planets"
        fig.savefig(final_image_name+"_final.png", format='png', bbox_inches='tight')
        plt.show()