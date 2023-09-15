import numpy as np
import matplotlib.pyplot as plt

def t_period_convert(t):
    """
    Convert a time value to its fractional part in the interval [0, 1].

    Args:
        t (float): The input time value.

    Returns:
        float: The fractional part of the input time.
    """
    # Check if the input time is a valid numeric type
    if not isinstance(t, (float, int, np.int64, np.float64)):
        print("time/period must be a positive number, auto set to initial position")
        return 1.0

    # Handle the special case when t is 0
    if t == 0:
        return 1.0

    # Handle the case when t is positive
    elif t > 0:
        integer_t = int(t)
        if np.isclose(t, integer_t):
            return 1.0
        else:
            return t-integer_t

    # Handle the case when t is negative
    else: #t<0
        integer_t = int(t)
        if np.isclose(t, integer_t):
            return 1.0
        else:
            return t+(1-integer_t)

def planet_position(a, e, pa, inc, t, mode, res):
    """
    Calculate the position of a planet in an elliptical orbit at a given time.

    Args:
        a (float): Semi-major axis of the orbit.
        e (float): Eccentricity of the orbit.
        pa (float): Position angle of the orbit (degrees).
        inc (float): Inclination of the orbit (degrees).
        t (float): Time value for which to calculate the planet's position.
        mode (str): Coordinate mode ("polar" or "cartesian") for the output.
        res (int): Number of points used for orbit calculations.

    Returns:
        tuple: Planet's position in the specified mode.
    """
    # Convert inputs from degree to rad
    pa = pa*2*np.pi/360.0
    inc = inc*2*np.pi/360.0

    # Keep the fractional_part of t
    t = t_period_convert(t)
    
    # Calculate the orbit
    theta = np.linspace(0, 2 * np.pi, res)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))
    
    # Calculate cumulated r^2
    cr2 = np.cumsum(r**2)

    # Find the nearest points
    target_theta = 0
    target_cr2 = t*cr2[-1]
    near_pos = np.searchsorted(cr2, target_cr2)
    if np.isclose(target_cr2, cr2[near_pos]):
        target_theta = theta[near_pos]
    else:
        theta_w1 = 1/(target_cr2 - cr2[near_pos-1])
        theta_w2 = 1/(target_cr2 - cr2[near_pos])
        target_theta = (theta[near_pos-1]*theta_w1 + theta[near_pos]*theta_w2)/(theta_w1+theta_w2)

    # Find the planet's location
    planet_r = a * (1 - e**2) / (1 + e * np.cos(target_theta))
    project_x, project_y = np.cos(target_theta+pa), np.sin(target_theta+pa)*np.cos(inc)
    if mode == "polar":
        return planet_r*abs(project_x+1j*project_y), np.arctan2(project_y, project_x)
    elif mode == "cartesian":
        return planet_r*project_x, planet_r*project_y
    else:
        print("Function planet_position got unknown mode")

def orbit_position(a, e, pa, inc, res):
    """
    Calculate the Cartesian positions of points along an elliptical orbit.

    Args:
        a (float): Semi-major axis of the orbit.
        e (float): Eccentricity of the orbit.
        pa (float): Position angle of the orbit (degrees).
        inc (float): Inclination of the orbit (degrees).
        res (int): Number of points used for orbit calculations.

    Returns:
        tuple: Arrays of x and y positions along the orbit.
    """
    # Convert angle inputs from degree to rad
    pa = pa*2*np.pi/360.0
    inc = inc*2*np.pi/360.0

    # Calculate the orbit
    theta = np.linspace(0, 2 * np.pi, res)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))

    # Find the projections
    orbit_x = r * np.cos(theta+pa)
    orbit_y = r * np.sin(theta+pa)*np.cos(inc)
    return orbit_x, orbit_y

def orbit_plot(a, e, pa, inc, planet_pos, plot_dpi, res, flip, mode, name=''):
    """
    Plot the orbit of a planet around a star.

    Args:
        a (float): Semi-major axis of the orbit.
        e (float): Eccentricity of the orbit.
        pa (float): Position angle of the orbit (degrees).
        inc (float): Inclination of the orbit (degrees).
        planet_pos (tuple): Planet's position (x, y) in the specified mode.
        mode (str): Coordinate mode ("polar" or "cartesian").
        name (str): Additional name for the saved plot file.
        plot_dpi (int): Dots per inch (DPI) for the plot.
        res (int): Number of points used for orbit calculations.

    Returns:
        orbit_*.png
    """
    # Convert position inputs to Cartesian coordinate system
    if mode == "polar":
        planet_pos = planet_pos[0]*np.cos(planet_pos[1]), planet_pos[0]*np.sin(planet_pos[1])
    elif mode == "cartesian":
        pass
    else:
        print("Function orbit_plot got unknown mode")
        
    # Convert angle inputs from degree to rad
    pa = pa*2*np.pi/360.0
    inc = inc*2*np.pi/360.0

    # Calculate the orbit
    theta = np.linspace(0, 2 * np.pi, res)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))

    # Find the projections
    orbit_x = r * np.cos(theta+pa)
    orbit_y = r * np.sin(theta+pa)*np.cos(inc)

    # Create a new figure with a specified DPI
    fig = plt.figure(dpi=plot_dpi)
    axs = plt.gca()

    # Set aspect ratio and axis labels, and invert the y-axis
    axs.set_aspect('equal', 'box')
    
    if flip:
        axs.invert_yaxis()
        
    axs.set_ylabel('y [arcsec]')
    axs.set_xlabel('x [arcsec]')
    
    # Plot the orbit in orange with a dashed line
    plt.plot(orbit_x, orbit_y, color='orange', linestyle = '--', alpha = 0.5) # Plot the orbit.

    # Mark the central star with a dark orange 'x'
    plt.scatter([0], [0], color='darkorange', marker='x')

    # Mark the planet with a black circle
    plt.scatter([planet_pos[0]], [planet_pos[1]], color='black', marker='o')
    
    # Save the figure with the specified name
    fig.savefig("oribit"+name+".png", format='png', bbox_inches='tight')

    # Display the plot
    plt.show()