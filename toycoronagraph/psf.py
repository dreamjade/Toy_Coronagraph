import numpy as np
import importlib.util
mp_spec = importlib.util.find_spec("multiprocessing")
if mp_spec is not None:
    import multiprocessing as mp
from hcipy import *
from skimage.transform import rotate
import matplotlib.pyplot as plt

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

def psf_calculation(charge, img_pixel, psf_range, lyot_mask_size, num_cores):
    """Calculate PSFs
    Calculates the PSFs of a vortex coronagraph for a given charge along the positive x-axis, and saves them to a file.

    Args:
        charge (int): Charge of the vortex coronagraph.
        img_pixel (int): Number of pixels on one edge of the image.
        psf_range (float): Range of the PSF in lambda/D.
        num_cores (int): Number of CPU cores to use for multiprocessing.

    Returns:
        np.ndarray: Array of calculated PSFs.
    """
    # Create grids and propagator
    pupil_grid = make_pupil_grid(img_pixel*2, 1.5)
    focal_grid_resolution = img_pixel/2/psf_range
    focal_grid = make_focal_grid(focal_grid_resolution, psf_range)
    prop = FraunhoferPropagator(pupil_grid, focal_grid)

    # Create a Lyot mask and coronagraph
    lyot_mask = evaluate_supersampled(make_circular_aperture(0.95), pupil_grid, 4)
    coro = VortexCoronagraph(pupil_grid, charge)
    lyot_stop = Apodizer(lyot_mask)

    # Normalization factor
    normal_factor = np.pi*(focal_grid_resolution/2)**2
    
    psfs = np.empty((img_pixel//2+1, img_pixel, img_pixel))
    # Check the existence of the multiprocessing module
    
    if mp_spec is not None:
        # Calculate PSFs in parallel using multiprocessing
        pool = mp.Pool(processes=num_cores)
        results = [pool.apply_async(psf_chunk, args=(i, img_pixel, psf_range, pupil_grid, prop, lyot_stop, coro)) for i in range(img_pixel//2+1)]
        pool.close()
        pool.join()
        for result in results:
            i, psf = result.get()
            psfs[i] = psf/normal_factor
    else:
        for i in range(img_pixel//2+1):
            # Calculate the x-coordinate based on the chunk index and PSF range
            x = 2*i*psf_range / img_pixel
        
            # Calculate the wavefront at the specified position
            wf = Wavefront(Wavefront_pos(x, 0,pupil_grid))
            
            # Propagate the wavefront through the system and calculate the intensity
            img = prop(lyot_stop(coro(wf))).intensity

            # Save results to psfs
            psfs[i] = img.to_dict()["values"].reshape(img_pixel, img_pixel)/normal_factor
        
    # Gather results and save PSFs to a file
    np.save('psfs_c'+str(charge)+'.npy', psfs)
    return psfs
    
def cir_psf_contrast(pre_img, planet_psfs_number, planet_angle, planet_brightness, psf_scale, img_pixel, rot_number, psfs_name):
    """Contrast
        
    Find the planet brightness and background brightness in the final image.

    Args:
        pre_img (np.ndarray): The input image.
        planet_psfs_number (int): According to the distance between the planet and the origin, choose the nearest psf.
        planet_angle (float): The angle between the x-axis and the line of the planet and the origin.
        planet_brightness (list): List of planet brightness values.
        psf_scale (float): Scale of the PSF.
        img_pixel (int): Number of pixels in the image.
        rot_number (int): Number of rotations for generating the circular PSF.
        psfs_name (str): Name of the file containing PSFs.

    Returns:
        brightness (tuple): planet brightness, background brightness, background brightness (ignored dust inside IWA) in Jy.
    """
    # Initialize an empty image chunk
    chunk_img = np.zeros([img_pixel, img_pixel])
    chunk_img_iwa = np.zeros([img_pixel, img_pixel]) # Ignore dust inside IWA

    # Load PSFs from the specified file
    psfs = np.load(psfs_name)
    
    # IWA
    iwa = cir_iwa(psfs)
    
    # Generate the chunk image along the x-axis
    for i in range(iwa):
        weight = pre_img[255+i][255]
        if weight != 0:
            chunk_img += 2*np.pi*i*weight*psfs[i]/rot_number
                
    for i in range(iwa, img_pixel//2+1):
        weight = pre_img[255+i][255]
        if weight != 0:
            part_img = 2*np.pi*i*weight*psfs[i]/rot_number
            chunk_img += part_img
            chunk_img_iwa += part_img
                
    # Initialize the disk image
    disk_img = np.zeros([img_pixel, img_pixel])
    disk_img_iwa = np.zeros([img_pixel, img_pixel])

    # Rotate and accumulate the chunk image to create the final image
    for i in range(rot_number):
        disk_img += rotate(chunk_img, angle=360*i/rot_number)
        disk_img_iwa += rotate(chunk_img_iwa, angle=360*i/rot_number)
        
    # The planet image
    planet_img = rotate(planet_brightness*psfs[planet_psfs_number], angle=planet_angle*180.0/np.pi)
    threshold = 0.5*planet_img.max()
    planet_filter = np.where(planet_img > threshold , 1, 0)

    # Brightness
    planet_b = np.sum(np.multiply(planet_filter, planet_img))*psf_scale**2
    dust_b = np.sum(np.multiply(planet_filter, disk_img))*psf_scale**2
    dust_b_iwa = np.sum(np.multiply(planet_filter, disk_img_iwa))*psf_scale**2
              
    return planet_b, dust_b, dust_b_iwa

def cir_psf(pre_img, planets_pos, planet_brightness, psf_scale, iwa_ignore, add_planet, img_pixel, rot_number, psfs_name):
    """Circular symmetric PSF processing calculation
    
    Calculates the final image of a circular symmetric pre-image through circular symmetric PSF with added planets.

    Args:
        pre_img (np.ndarray): The input image.
        planets_pos (list): List of planet positions.
        planet_brightness (list): List of planet brightness values.
        psf_scale (float): Scale of the PSF.
        add_planet (bool): Whether to add planets to the PSF.
        img_pixel (int): Number of pixels in the image.
        rot_number (int): Number of rotations for generating the circular PSF.
        psfs_name (str): Name of the file containing PSFs.

    Returns:
        np.ndarray: Final image.
    """
    # Load PSFs from the specified file
    psfs = np.load(psfs_name)
    
    # IWA
    iwa = 0
    if iwa_ignore:
        iwa = cir_iwa(psfs)
    
    # Initialize the final image
    final_img = np.zeros([img_pixel, img_pixel])
    
    # Initialize an empty image chunk
    chunk_img = np.zeros([img_pixel, img_pixel])
    
    # Generate the chunk image along the x-axis
    for i in range(iwa, img_pixel//2+1):
        weight = pre_img[255+i][255]
        if weight != 0:
            chunk_img += 2*np.pi*i*weight*psfs[i]/rot_number

    # Rotate and accumulate the chunk image to create the final image
    for i in range(rot_number):
        final_img += rotate(chunk_img, angle=360*i/rot_number)
        
    # Add planets if requested
    if add_planet:
        for i in range(len(planets_pos)):
            psfs_number = int(planets_pos[i][0]/psf_scale)
            if psfs_number<img_pixel/2:
                final_img += rotate(planet_brightness[i]*psfs[psfs_number], angle=planets_pos[i][1]*180.0/np.pi)
            else:
                print("Planet #"+str(i+1)+" is outside the range of the plot")
            
    return final_img
    
def cir_psf_planets(planets_pos, planet_brightness, psf_scale, img_pixel, psfs_name):
    """
    Add one or more planets to an image through a series of circular symmetric PSFs.

    Args:
        planets_pos (list): A list of planet positions, in units of arcsec.
        planet_brightness (list): A list of planet brightnesses, in units of Jy.
        psf_scale (float): The scale of the PSF, in units of arcsec.
        img_pixel (int): The number of pixels in the image.
        psfs_name (str): The file name of the PSFs.

    Returns:
        np.ndarray: The image of planets.
    """
    # Load PSFs from the specified file
    psfs = np.load(psfs_name)
    
    # Initialize the final image
    planet_img = np.zeros([img_pixel, img_pixel])

    for i in range(len(planets_pos)):
        psfs_number = int(planets_pos[i][0]/psf_scale)
        if psfs_number<img_pixel/2:
            planet_img += rotate(planet_brightness[i]*psfs[psfs_number], angle=planets_pos[i][1]*180.0/np.pi)
    return planet_img
    
def core_throughput_vectorized(psfs):
    """
    Calculate the core throughput of a series of PSFs in a vectorized way.

    Args:
        psfs (numpy.ndarray): The normalized PSFs.

    Returns:
        numpy.ndarray: The core throughput of the psfs
    """
    thresholds = 0.5 * np.max(psfs, axis=(1, 2)).reshape(-1, 1, 1)
    #core = (psfs > thresholds).astype(int)
    #core_flux = core * psfs
    core_flux = np.where(psfs > thresholds , psfs, 0)
    core_throughput = np.sum(core_flux, axis=(1, 2))
    return core_throughput
        
def core_throughput_fun(psfs):
    """
    Calculate the core throughput of a series of PSFs.

    Args:
        psfs (numpy.ndarray): The normalized PSFs.

    Returns:
        numpy.ndarray: The core throughput of the psfs
    """
    psf_len = len(psfs)
    core_throughput = np.empty(psf_len)
    for i in range(psf_len):
        psf_i = psfs[i]
        threshold = 0.5*psf_i.max()
        core_flux = np.where(psf_i > threshold , psf_i, 0)
        core_throughput[i] = np.sum(core_flux)
    return core_throughput

def cir_iwa(psfs):
    """
    Calculate the inner working angle (IWA) of a series of circular symmetric PSFs.

    Args:
        psfs (numpy.ndarray): The normalized PSFs.

    Returns:
        int: The IWA of the PSFs, in units of pixel.
    """
    core_throughput = core_throughput_fun(psfs)
    iwa = np.searchsorted(core_throughput, 0.5*core_throughput.max())
    return iwa

def cir_core_throughput_plot(psfs, plot_dpi=300):
    """
    Plot the core throughput of a series of circular symmetric PSFs.

    Args:
        psfs (numpy.ndarray): The normalized PSFs.
        plot_dpi (int): Dots per inch (DPI) for the plot.

    Return:
        Core_throughput.png
    """
    core_throughput = core_throughput_fun(psfs)
    iwa = np.searchsorted(core_throughput, 0.5*core_throughput.max())
    # Create a new figure with a specified DPI
    fig = plt.figure(dpi=plot_dpi)
    axs = plt.gca()

    # Set axis labels
    axs.set_ylabel('Core Throughput')
    axs.set_xlabel('distance [pixel]')
    
    # Plot the core throughput
    plt.plot(range(len(core_throughput)), core_throughput, color='orange')

    # Mark the IWA with a black verticle line
    plt.axvline(iwa, color='black', linestyle='--')
    
    # Save the figure with the specified name
    fig.savefig("Core_throughput.png", format='png', bbox_inches='tight')

    # Display the plot
    plt.show()