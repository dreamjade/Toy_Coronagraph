import numpy as np
import matplotlib.pyplot as plt
import hcipy
from .psf import psf_calculation, cir_psf
from astropy.io import fits

class Target(object):
    """
    A circular symmetric target, which could be a dust ring, a debris disk, or ice remains
    """
    def __init__(self, file_name="example"):
        self.file_name = file_name
        inc = fits.open(file_name+".fits")
        px=512
        py=512
        psf_scale=1e-6/2.4*206264.806*32/512 ##arcsecs/pixel
        xpix=(np.arange (-px/2, px/2, 1))*psf_scale
        ypix=(np.arange (-px/2, px/2, 1))*psf_scale
        sq_as_per_pix=psf_scale**2
        c=2.99792*10**14
        wave_length=1.0 #in microns#
        jy=10**26
        sst=np.reshape(inc[0].data[5,0],(px,py))
        self.data_jy=(sst/c)*(wave_length**2)*jy
        
    def plot_origin(self):
        fig=plt.figure(dpi=300)
        ax=plt.subplot(111)
        im=ax.imshow(sst_jy.astype(np.float64),
                       cmap='gnuplot',extent=[np.min(ypix),np.max(ypix),np.min(xpix),np.max(xpix)])
        ax.invert_yaxis()
        ax.set_ylabel('y [arcsec]')
        ax.set_xlabel('x [arcsec]')
        cb=plt.colorbar(im,orientation='vertical')
        cb.set_label("$Jy$")
        #fig.savefig(self.file_name+".png", format='png', bbox_inches='tight')
        plt.show()
        
    def plot_final(self, charge, plot_dpi=300):
        img_pixel = 512
        psf_range = 16
        rot_number = 360

        pre_img = sst_jy.astype(np.float64)
        final_img_charge = cir_psf(pre_img, img_pixel, psf_range, img_pixel, 
                                    "psfs_c"+str(charge)+"10.npy")
        fig=plt.figure(dpi=plot_dpi)
        ax2=plt.subplot(111)
        im2=ax2.imshow(final_img_charge,
                   cmap='gnuplot',extent=[np.min(ypix),np.max(ypix),np.min(xpix),np.max(xpix)])
        ax2.invert_yaxis()
        ax2.set_ylabel('y [arcsec]')
        ax2.set_xlabel('x [arcsec]')
        cb=plt.colorbar(im2,orientation='vertical')
        cb.set_label("$Jy$")
        #fig.savefig(self.file_name+"_charge"+str(charge)+"_final.png", format='png', bbox_inches='tight')
        plt.show()

if __name__ == "__main__":
    main()