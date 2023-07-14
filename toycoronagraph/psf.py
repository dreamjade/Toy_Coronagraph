import numpy as np
import multiprocessing as mp
!set NUMEXPR_MAX_THREADS = 64
from hcipy import *

def Wavefront_pos(x,y,pupil_grid):
    aperture = evaluate_supersampled(make_circular_aperture(1), pupil_grid, 4)
    return aperture * np.exp(2j * np.pi * (pupil_grid.x * x + pupil_grid.y * y))

def psf_chunk(i_psf_range, img_pixel, pupil_grid, prop, lyot_stop, coro):
    x = 2*i_psf_range / img_pixel
    wf = Wavefront(Wavefront_pos(x, 0,pupil_grid))
    img = prop(lyot_stop(coro(wf))).intensity
    return i, img.to_dict()["values"].reshape(img_pixel, img_pixel)

def psf_calculation(charge, img_pixel=512, psf_range=16, num_cores = 4):
    pupil_grid = make_pupil_grid(1024, 1.5)
    focal_grid = make_focal_grid(16, 16)
    prop = FraunhoferPropagator(pupil_grid, focal_grid)
    lyot_mask = evaluate_supersampled(make_circular_aperture(0.5), pupil_grid, 4)
    coro = VortexCoronagraph(pupil_grid, charge=10)
    lyot_stop = Apodizer(lyot_mask)
    
    chunk_size = img_pixel // (2*num_cores)
    psfs = np.empty((img_pixel//2+1, img_pixel, img_pixel))
    pool = mp.Pool(processes=num_cores)
    results = [pool.apply_async(psf_chunk, args=(i*psf_range, img_pixel, pupil_grid, prop, lyot_stop, coro)) for i in range(img_pixel//2+1)]
    pool.close()
    pool.join()
    for result in results:
        i, psf = result.get()
        psfs[i] = psf
    np.save('psfs_c'+str(charge)+'.npy', psfs)
    return psfs

def cir_psf(pre_img, img_pixel=512, psf_range=16, rot_number=360, psfs_name=None):
    chunk_img = np.zeros([img_pixel, img_pixel])
    if psfs_name==None:
        for i in range(img_pixel//2+1):
            x = 2*i*psf_range / img_pixel
            weight = pre_img[255+i][255]
            if weight != 0:
                wf = Wavefront(Wavefront_pos(x, 0))
                img = prop(lyot_stop(coro(wf))).intensity
                chunk_img += 2*np.pi*i*weight*img.to_dict()["values"].reshape(img_pixel, img_pixel)/rot_number
    else:
        psfs = np.load(psfs_name)
        for i in range(img_pixel//2+1):
            weight = pre_img[255+i][255]
            if weight != 0:
                chunk_img += 2*np.pi*i*weight*psfs[i]/rot_number
    final_img = np.zeros([img_pixel, img_pixel])
    for i in range(rot_number):
        final_img += rotate(chunk_img, angle=360*i/rot_number)
    return final_img