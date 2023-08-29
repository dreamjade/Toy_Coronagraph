import numpy as np
import multiprocessing as mp
from hcipy import *
from skimage.transform import rotate

def Wavefront_pos(x,y,pupil_grid):
    """
    Compute the wavefront with a given position.

    Args:
        x (float): x-coordinate of the wavefront position.
        y (float): y-coordinate of the wavefront position.
        pupil_grid: Grid representing the pupil plane.

    Returns:
        Field: Wavefront at the specified position.
    """
    # Create an aperture using the hcipy library
    aperture = evaluate_supersampled(make_circular_aperture(1), pupil_grid, 4)

    # Calculate the wavefront with the specified position
    return aperture * np.exp(2j * np.pi * (pupil_grid.x * x + pupil_grid.y * y))

def psf_chunk(i, img_pixel, psf_range, pupil_grid, prop, lyot_stop, coro):
    """
    Calculate a chunk of PSF.

    Args:
        i (int): Index of the chunk.
        img_pixel (int): Number of pixels in the image.
        psf_range (float): Range of the PSF.
        pupil_grid: Grid representing the pupil plane.
        prop: Propagator for the PSF calculation.
        lyot_stop: Lyot stop for the coronagraph.
        coro: Coronagraph for the PSF calculation.

    Returns:
        tuple: Index and the calculated PSF chunk as a 2D numpy array.
    """
    # Calculate the x-coordinate based on the chunk index and PSF range
    x = 2*i*psf_range / img_pixel

    # Calculate the wavefront at the specified position
    wf = Wavefront(Wavefront_pos(x, 0,pupil_grid))

    # Propagate the wavefront through the system and calculate the intensity
    img = prop(lyot_stop(coro(wf))).intensity

    # Reshape and return the calculated PSF chunk
    return i, img.to_dict()["values"].reshape(img_pixel, img_pixel)

def psf_calculation(charge, img_pixel=512, psf_range=16, num_cores = 16):
    """Calculate PSFs
    Calculates the PSFs of a vortex coronagraph for a given charge along the positive x-axis, and saves them to a file.

    Args:
        charge (int): Charge of the vortex coronagraph.
        img_pixel (int): Number of pixels in the image.
        psf_range (float): Range of the PSF.
        num_cores (int): Number of CPU cores to use for multiprocessing.

    Returns:
        np.ndarray: Array of calculated PSFs.
    """
    # Create grids and propagator
    pupil_grid = make_pupil_grid(1024, 1.5)
    focal_grid = make_focal_grid(16, 16)
    prop = FraunhoferPropagator(pupil_grid, focal_grid)

    # Create Lyot mask and coronagraph
    lyot_mask = evaluate_supersampled(make_circular_aperture(0.5), pupil_grid, 4)
    coro = VortexCoronagraph(pupil_grid, charge)
    lyot_stop = Apodizer(lyot_mask)

    # Calculate PSFs in parallel using multiprocessing
    chunk_size = img_pixel // (2*num_cores)
    psfs = np.empty((img_pixel//2+1, img_pixel, img_pixel))
    pool = mp.Pool(processes=num_cores)
    results = [pool.apply_async(psf_chunk, args=(i, img_pixel, psf_range, pupil_grid, prop, lyot_stop, coro)) for i in range(img_pixel//2+1)]
    pool.close()
    pool.join()
    for result in results:
        i, psf = result.get()
        psfs[i] = psf
    
    # Gather results and save PSFs to a file
    np.save('psfs_c'+str(charge)+'.npy', psfs)
    return psfs

def cir_psf(pre_img, planets_pos, planet_brightness, psf_scale, add_planet, img_pixel=512, psf_range=16, rot_number=360, psfs_name="psfs_c2.npy"):
    """Circular symmetric PSF processing calculation
    Calculates the final image of a circular symmetric pre-image through circular symmetric PSF with added planets.

    Args:
        pre_img (np.ndarray): The input image.
        planets_pos (list): List of planet positions.
        planet_brightness (list): List of planet brightness values.
        psf_scale (float): Scale of the PSF.
        add_planet (bool): Whether to add planets to the PSF.
        img_pixel (int): Number of pixels in the image.
        psf_range (float): Range of the PSF.
        rot_number (int): Number of rotations for generating the circular PSF.
        psfs_name (str): Name of the file containing PSFs.

    Returns:
        np.ndarray: Final image.
    """
    # Initialize an empty image chunk
    chunk_img = np.zeros([img_pixel, img_pixel])

    # Load PSFs from the specified file
    psfs = np.load(psfs_name)
    for i in range(img_pixel//2+1):
        weight = pre_img[255+i][255]
        if weight != 0:
            chunk_img += 2*np.pi*i*weight*psfs[i]/rot_number

    # Initialize the final image
    final_img = np.zeros([img_pixel, img_pixel])

    # Rotate and accumulate the chunk image to create the final image
    for i in range(rot_number):
        final_img += rotate(chunk_img, angle=360*i/rot_number)
        
    # Add planets if requested
    if add_planet:
        for i in range(len(planets_pos)):
            psfs_number = int(planets_pos[i][0]/psf_scale)
            if psfs_number<img_pixel/2:
                final_img += rotate(planet_brightness[i]*psfs[psfs_number], angle=-planets_pos[i][1]*180.0/np.pi)
            else:
                print("Planet #"+str(i+1)+" is outside the range of the plot")
            
    return final_img