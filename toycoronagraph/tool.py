import numpy as np
from skimage.transform import warp, SimilarityTransform
from scipy.spatial import distance

def convert_to_polar(pos_car):
    """
    Converts cartesian coordinates to polar coordinates.

    Args:
        pos_car (numpy.ndarray): A numpy array containing the x-coordinate and y-coordinate of the point.

    Returns:
        numpy.ndarray: A numpy array containing the radius and angle of the point.
    """
    r = np.sqrt(pos_car[0]**2 + pos_car[1]**2)
    theta = np.arctan2(pos_car[1], pos_car[0])
    return np.array([r, theta])

def is_positive_even_integer(number):
    """
    Checks if the input number is a positive even integer.

    Args:
        number (int): The number to check.

    Returns:
        bool: True if the number is a positive even integer, False otherwise.
    """
    if (not isinstance(number, (int, np.int64)) 
        or number <= 0):
        return False
    if number % 2 != 0:
        return False
    return True

def is_planet_pos_allowed(pos, mode):
    """
    Checks if the input planet position is valid based on the specified mode.

    Args:
        pos (list or numpy.ndarray): The planet position.
        mode (str): The mode to determine the validity of the position.

    Returns:
        bool: True if the planet position is valid, False otherwise.
    """
    if mode == "moving":
        if len(pos)!=5:
            print("The motion of a planet is described by five parameters: the semi-major axis, the eccentricity, the position angle, the inclination, and the time in units of the period.")
            return False
        if (not isinstance(pos[0], (float, int, np.float64, np.int64))
            or pos[0]<=0
            or not isinstance(pos[1], (float, int, np.float64, np.int64))
            or pos[1]<0
            or pos[1]>=1
            or not isinstance(pos[2], (float, int, np.float64, np.int64))
            or not isinstance(pos[3], (float, int, np.float64, np.int64))
            or not isinstance(pos[4], (float, int, np.float64, np.int64))):
            return False
        return True
    elif mode == "cartesian":
        if len(pos)!=2:
            return False
        if (not isinstance(pos[0], (float, int, np.float64, np.int64)) 
            or not isinstance(pos[1], (float, int, np.float64, np.int64))):
            return False
        return True
    elif mode == "polar":
        if len(pos)!=2:
            return False
        if (not isinstance(pos[0], (float, int, np.float64, np.int64)) 
            or not isinstance(pos[1], (float, int, np.float64, np.int64)) 
            or pos[0]<=0):
            return False
        return True 
    
    else:
        print("unknown add planet mode")

def frame_interpolation(positions, images, pos):
    """
    Find an intermediate image between images using optical flow methods.
    """
    pos = np.array(pos).reshape(1, -1)
    positions = np.array(positions).reshape(-1, np.size(pos))
    distances = distance.cdist(positions, pos, 'euclidean').flatten()
    rel_pos = positions-pos
    if np.size(pos)==1:
        rel_pos = np.concatenate((rel_pos, np.zeros(rel_pos.shape)), axis=1)
    total_weight = 0
    final_image = np.zeros(np.shape(images[0]))
    for i in range(len(distances)):
        d=distances[i]
        if d==0:
            return images[i]
        else:
            weight = 1/d
            final_image = final_image + weight*warp(images[i], SimilarityTransform(translation=rel_pos[i]))
            total_weight = total_weight + weight
    return final_image/total_weight