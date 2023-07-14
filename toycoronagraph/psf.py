import numpy as np
import matplotlib.pyplot as plt
import hcipy

def psf_calculation(charge, pixel_size):
    """
    Calculates the PSF for a given charge and pixel size.

    Args:
        charge: The number of vortices in the coronagraph.
        pixel_size: The size of each pixel in the image.

    Returns:
        A NumPy array that contains the PSF.
    """

    psfs = np.empty((pixel_size//2+1, pixel_size, pixel_size))
    for i in range(pixel_size//2+1):
        x = 2*i*pixel_size / pixel_size
        wf = Wavefront(Wavefront_pos(x, 0))
        img = prop(lyot_stop(coro(wf))).intensity
        psfs[i] = img.to_dict()["values"].reshape(pixel_size, pixel_size)
    return psfs

def cir_psf(pre_img, pixel_size, psf_range, rot_number):
    """
    Simulates the PSF by summing the PSFs at different positions.

    Args:
        pre_img: The image that the PSF will be simulated on.
        pixel_size: The size of each pixel in the image.
        psf_range: The range of positions that will be used to calculate the PSF.
        rot_number: The number of times the PSF will be rotated.

    Returns:
        A NumPy array that contains the simulated PSF.
    """

    chunk_img = np.zeros([pixel_size, pixel_size])
    for i in range(pixel_size//2+1):
        weight = pre_img[255+i][255]
        if weight != 0:
            chunk_img += 2*np.pi*i*weight*psfs[i]/rot_number
    final_img = np.zeros([pixel_size, pixel_size])
    for i in range(rot_number):
        final_img += rotate(chunk_img, angle=360*i/rot_number)
    return final_img