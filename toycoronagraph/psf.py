import numpy as np
import multiprocessing as mp
from hcipy import *
from skimage.transform import rotate

def Wavefront_pos(x,y,pupil_grid):
    """
    Calculates the wavefront at a given position.

    Args:
        x: The x-coordinate of the position.
        y: The y-coordinate of the position.
        pupil_grid: The pupil grid.

    Returns:
        The wavefront at the given position.
    """
    aperture = evaluate_supersampled(make_circular_aperture(1), pupil_grid, 4)
    return aperture * np.exp(2j * np.pi * (pupil_grid.x * x + pupil_grid.y * y))

def psf_chunk(i, img_pixel, psf_range, pupil_grid, prop, lyot_stop, coro):
    """
    Calculates a single PSF chunk.

    Args:
        i: The index of the chunk.
        img_pixel: The number of pixels in the image.
        psf_range: The range of the PSF in pixels.
        pupil_grid: The pupil grid.
        prop: The FraunhoferPropagator.
        lyot_stop: The Apodizer.
        coro: The VortexCoronagraph.

    Returns:
        The index of the chunk and the PSF chunk.
    """
    x = 2*i*psf_range / img_pixel
    wf = Wavefront(Wavefront_pos(x, 0,pupil_grid))
    img = prop(lyot_stop(coro(wf))).intensity
    return i, img.to_dict()["values"].reshape(img_pixel, img_pixel)

def psf_calculation(charge, img_pixel=512, psf_range=16, num_cores = 16):
    """
    Calculates the PSFs of a vortex coronagraph for a given charge.

    Args:
        charge: The vortex charge.
        img_pixel: The number of pixels in the image.
        psf_range: The range of the PSF in pixels.
        num_cores: The number of cores to use.

    Returns:
        The PSFs.
    """
    pupil_grid = make_pupil_grid(1024, 1.5)
    focal_grid = make_focal_grid(16, 16)
    prop = FraunhoferPropagator(pupil_grid, focal_grid)
    lyot_mask = evaluate_supersampled(make_circular_aperture(0.5), pupil_grid, 4)
    coro = VortexCoronagraph(pupil_grid, charge)
    lyot_stop = Apodizer(lyot_mask)
    
    chunk_size = img_pixel // (2*num_cores)
    psfs = np.empty((img_pixel//2+1, img_pixel, img_pixel))
    pool = mp.Pool(processes=num_cores)
    results = [pool.apply_async(psf_chunk, args=(i, img_pixel, psf_range, pupil_grid, prop, lyot_stop, coro)) for i in range(img_pixel//2+1)]
    pool.close()
    pool.join()
    for result in results:
        i, psf = result.get()
        psfs[i] = psf
    
    #save psfs into npy file
    np.save('psfs_c'+str(charge)+'.npy', psfs)
    return psfs

def cir_psf(pre_img, planets_pos, planet_brightness, add_planet, img_pixel=512, psf_range=16, rot_number=360, psfs_name="psfs_c2.npy"):
    """
    Calculates the final image of a circular symmetric pre-image through circular symmetric PSF

    Args:
        pre_img: The pre-image.
        img_pixel: The number of pixels in the image.
        psf_range: The range of the PSF in pixels.
        rot_number: The number of rotations.
        psfs_name: The name of the PSF file.

    Returns:
        The final image through the coronagraph.
    """
    chunk_img = np.zeros([img_pixel, img_pixel])
    psfs = np.load(psfs_name)
    for i in range(img_pixel//2+1):
        weight = pre_img[255+i][255]
        if weight != 0:
            chunk_img += 2*np.pi*i*weight*psfs[i]/rot_number
    final_img = np.zeros([img_pixel, img_pixel])
    for i in range(rot_number):
        final_img += rotate(chunk_img, angle=360*i/rot_number)
        
    #add planets
    if add_planet:
        for i in range(len(planets_pos)):
            final_img += rotate(planet_brightness[i]*psfs[int(planets_pos[i][0])], angle=360*planets_pos[i][1]/2/np.pi)
            
    return final_img