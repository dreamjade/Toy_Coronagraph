import numpy as np
import matplotlib.pyplot as plt
import hcipy
from .psf import psf_calculation, cir_psf
from astropy.io import fits
import os

class Target(object):
    """
    A circular symmetric target, which could be a dust ring, a debris disk, or ice remains
    """
    def __init__(self, file_name="example"):
        self.file_name = file_name
        inc = fits.open(file_name+".fits")
        self.px=512
        self.py=512
        self.psf_scale=1e-6/2.4*206264.806*32/512 ##arcsecs/pixel
        self.xpix=(np.arange (-self.px/2, self.px/2, 1))*self.psf_scale
        self.ypix=(np.arange (-self.px/2, self.px/2, 1))*self.psf_scale
        c=2.99792*10**14
        wave_length=1.0 #in microns#
        jy=10**26
        sst=np.reshape(inc[0].data[5,0],(self.px,self.py))
        self.data_jy=(sst/c)*(wave_length**2)*jy
        self.pre_img = self.data_jy.astype(np.float64)
        
    def plot_origin(self):
        fig=plt.figure(dpi=300)
        ax=plt.subplot(111)
        im=ax.imshow(self.pre_img,
                       cmap='gnuplot',extent=[np.min(self.ypix),np.max(self.ypix),np.min(self.xpix),np.max(self.xpix)])
        ax.invert_yaxis()
        ax.set_ylabel('y [arcsec]')
        ax.set_xlabel('x [arcsec]')
        cb=plt.colorbar(im,orientation='vertical')
        cb.set_label("$Jy$")
        #fig.savefig(self.file_name+".png", format='png', bbox_inches='tight')
        plt.show()
        
    def plot_final(self, charge=6, img_pixel = 512, psf_range = 16, rot_number = 360, plot_dpi=300):
        psf_filename = "psfs_c"+str(charge)+".npy"
        if not os.path.exists(psf_filename):
            psf_calculation(charge, img_pixel, psf_range)
        final_img_charge = cir_psf(self.pre_img, img_pixel, psf_range, img_pixel, psf_filename)
        fig=plt.figure(dpi=plot_dpi)
        ax2=plt.subplot(111)
        im2=ax2.imshow(final_img_charge,
                   cmap='gnuplot',extent=[np.min(self.ypix),np.max(self.ypix),np.min(self.xpix),np.max(self.xpix)])
        ax2.invert_yaxis()
        ax2.set_ylabel('y [arcsec]')
        ax2.set_xlabel('x [arcsec]')
        cb=plt.colorbar(im2,orientation='vertical')
        cb.set_label("$Jy$")
        #fig.savefig(self.file_name+"_charge"+str(charge)+"_final.png", format='png', bbox_inches='tight')
        plt.show()

if __name__ == "__main__":
    main()