import numpy as np

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

    if not isinstance(number, (int, np.int64)) or number <= 0:
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
        if not isinstance(pos[0], (float, int, np.float64, np.int64)) or not isinstance(pos[1], (float, int, np.float64, np.int64)):
            return False
        return True
    elif mode == "polar":
        if len(pos)!=2:
            return False
        if not isinstance(pos[0], (float, int, np.float64, np.int64)) or not isinstance(pos[1], (float, int, np.float64, np.int64)) or pos[0]<=0:
            return False
        return True 
    
    else:
        print("unknown add planet mode")