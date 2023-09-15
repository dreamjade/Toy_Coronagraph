import os

def create_py_file(file_name, data):
    """
    Creates a .py file with the given name and data.

    Args:
        file_name (str): The name of the .py file.
        data (str): The data to be included in the .py file.

    Returns:
        *.py
    """

    with open(file_name, "w") as f:
        f.write(data)

def example_para():
    """
    Create an example parameter file named 'toycoronagraph_para.py'.

    The example parameter file includes settings related to the toy coronagraph simulation.

    Returns:
        toycoronagraph_para.py
    """
    # Define the name of the parameter file and the data to be written
    file_name = "toycoronagraph_para.py"
    data = \
'''
fits_filename = None # The address of the FITS file containing the target data. If None, the default file `DATADIR/example.fit` will be used.
taregt_pixel = 512 # The pixel number of the target image in the x-direction and y-direction
num_cores = 16 # The number of cores used in parallel calculation.
rot_number = taregt_pixel # The number of rotations for generating the circular symmetry final image.

# Calculate the psf scale in arcsecs per pixel
wavelength = 1.0e-6 # meter
D = 2.4 # aperture size in meter
rad_to_arcsecond = 206264.806247
visual_range = 32 # in lambda/D
psf_scale = wavelength/D*rad_to_arcsecond*visual_range/taregt_pixel # arcsec per pixel

# Example .fits file has unit W/m^2/pixel
# vF_v(W/m^2/pixel) to F_v(Jy/arcsec^2)
light_speed = 299792458 # m/s
jy=10**26 # The conversion factor from W / m^2 / sr / Hz to Jy
F_transfer = jy*wavelength/light_speed/psf_scale**2

# Coronagraph setting
coronagraph_type='vortex'
lyot_mask_size = 0.95
psf_range=16
'''
    # Create the parameter file using the defined name and data
    create_py_file(file_name, data)