import numpy as np
import matplotlib.pyplot as plt
import hcipy
#from .psf import psf_calculation, cir_psf
from toycoronagraph.psf import psf_calculation, cir_psf
from astropy.io import fits
from toycoronagraph import DATADIR
import os
def convert_to_polar(pos_car):
    """Converts cartesian coordinates to polar coordinates.
    Args:
        A numpy array of the x-coordinate and y-coordinate of the point.

    Returns:
        A numpy array of the radius and angle of the point.
    """
    r = np.sqrt(pos_car[0]**2 + pos_car[1]**2)
    theta = np.arctan2(pos_car[1], pos_car[0])
    return np.array([r, theta])

def is_positive_even_integer(number):
    """
    Checks if the input number is a positive even integer.

    Args:
        number: The number to check.

    Returns:
        True if the number is a positive even integer, False otherwise.
    """

    if not isinstance(number, int) or number <= 0:
        return False
    if number % 2 != 0:
        return False
    return True

def is_planet_pos_allowed(pos, mode):
    """
    Checks if the input planet position is on the plane.
    """
    if mode == "cartesian":
        if len(pos)!=2:
            return False
        if not isinstance(pos[0], (float, int)) or not isinstance(pos[1], (float, int)):
            return False
        return True
    elif mode == "polar":
        if len(pos)!=2:
            return False
        if not isinstance(pos[0], (float, int)) or not isinstance(pos[1], (float, int)) or pos[0]<=0:
            return False
        return True 
    
    else:
        print("unknown add planet mode")

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
        
    def add_planet(self, pos, brightness, mode="polar"):
        pos = np.array(pos)
        
        if mode=="polar":
            if not is_planet_pos_allowed(pos, mode):
                print("Planet position is not allowed")
            elif not isinstance(brightness, (float, int)) or brightness<=0:
                print("Brightness is invalid")
            else:
                self.planets.append(pos)
                self.planets_brightness.append(brightness)
            
        elif mode=="cartesian":
            if not is_planet_pos_allowed(pos, mode):
                print("planet position is not allowed")
            elif not isinstance(brightness, (float, int)) or brightness<=0:
                print("Brightness is invalid")
            else:
                self.planets.append(convert_to_polar(pos))
                self.planets_brightness.append(brightness)
        else:
            print("no "+str(mode)+" mode")
        
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
            psf_calculation(charge, img_pixel, psf_range)
            psf_filename = "psfs_c"+str(charge)+".npy"
        
        #draw final image of the disk
        final_img_charge = cir_psf(self.pre_img, self.planets, add_planet, img_pixel, psf_range, img_pixel, psf_filename)
        
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
        fig.savefig("charge"+str(charge)+"_final.png", format='png', bbox_inches='tight')
        plt.show()
