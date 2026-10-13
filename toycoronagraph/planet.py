"""
Planet orbits: the projected position of a planet on an elliptical orbit at a given time, and orbit plots.

The orbit is an ellipse with the star at a focus, r(theta) = a(1-e^2)/(1+e cos(theta)), with
theta the true anomaly (theta = 0 at periapsis). It is rotated by the position angle pa in its
plane and then inclined about the x axis, so the sky position is
x = r cos(theta+pa), y = r sin(theta+pa) cos(inc). Lengths are in the unit of a (arcsec in the package).
"""
import numpy as np
import matplotlib.pyplot as plt

def t_period_convert(t):
    """
    Convert a time value to its fractional part in the interval (0, 1].

    Whole numbers of periods (including 0, up to np.isclose) map to 1.0 (the initial position),
    not to 0; negative times count back from the initial position (e.g. -0.3 -> 0.7).

    Args:
        t (float): The input time, in units of the orbital period.

    Returns:
        float: The fractional part of the input time, in (0, 1]; 1.0 (after printing a message)
            if t is not an int/float.
    """
    # Check if the input time is a valid numeric type
    if not isinstance(t, (float, int, np.int64, np.float64)):
        print("time/period must be a positive number, auto set to initial position")
        return 1.0

    if t == 0:
        return 1.0

    # int() truncates towards zero, so the fractional part of a negative t is t - int(t) + 1
    integer_t = int(t)
    if np.isclose(t, integer_t):
        return 1.0
    if t > 0:
        return t-integer_t
    return t+(1-integer_t)

def planet_position(a, e, pa, inc, t, mode, res):
    """
    Calculate the position of a planet in an elliptical orbit at a given time.

    Kepler's second law: the swept area grows as the integral of r^2 dtheta, so the true anomaly
    at time t is where the cumulative sum of r^2 over res samples of theta in [0, 2 pi] reaches
    the fraction t of its total (linear interpolation between the two nearest samples). t = 1
    (or any whole number of periods) is the periapsis.

    Args:
        a (float): Semi-major axis of the orbit (arcsec).
        e (float): Eccentricity of the orbit, 0 <= e < 1.
        pa (float): Position angle of the orbit (degrees): rotation of the periapsis from +x in the orbital plane.
        inc (float): Inclination of the orbit (degrees), about the x axis (0 = face-on).
        t (float): Time in units of the orbital period (see t_period_convert).
        mode (str): Coordinate mode ("polar" or "cartesian") for the output.
        res (int): Number of points used for orbit calculations.

    Returns:
        tuple or None: Planet's projected position: (r in arcsec, angle in radian from +x towards +y)
            for "polar", (x, y) in arcsec for "cartesian"; None (after printing a message) for any other mode.
    """
    # Convert inputs from degree to rad
    pa = pa*2*np.pi/360.0
    inc = inc*2*np.pi/360.0

    # Keep the fractional_part of t
    t = t_period_convert(t)

    # Calculate the orbit
    theta = np.linspace(0, 2 * np.pi, res)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))

    # Calculate cumulated r^2 (proportional to the area swept since periapsis)
    cr2 = np.cumsum(r**2)

    # Find the true anomaly where the swept area is the fraction t of the total:
    # the sample itself if it matches, else inverse-distance (i.e. linear) interpolation between neighbours
    target_cr2 = t*cr2[-1]
    near_pos = np.searchsorted(cr2, target_cr2)
    if np.isclose(target_cr2, cr2[near_pos]):
        target_theta = theta[near_pos]
    else:
        theta_w1 = 1/(target_cr2 - cr2[near_pos-1])
        theta_w2 = 1/(cr2[near_pos] - target_cr2)
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
        a (float): Semi-major axis of the orbit (arcsec).
        e (float): Eccentricity of the orbit, 0 <= e < 1.
        pa (float): Position angle of the orbit (degrees), see planet_position.
        inc (float): Inclination of the orbit (degrees), see planet_position.
        res (int): Number of points used for orbit calculations.

    Returns:
        tuple: orbit_x, orbit_y: (res,) float64 arrays of the projected x and y positions (arcsec),
            at true anomalies evenly spaced over [0, 2 pi] (first and last points coincide).
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
    Plot the orbit of a planet around a star, save it to "oribit<name>.png" and show it.

    Args:
        a (float): Semi-major axis of the orbit (arcsec).
        e (float): Eccentricity of the orbit, 0 <= e < 1.
        pa (float): Position angle of the orbit (degrees), see planet_position.
        inc (float): Inclination of the orbit (degrees), see planet_position.
        planet_pos (tuple): Planet's position in the specified mode: (r arcsec, angle radian) for
            "polar", (x, y) arcsec for "cartesian".
        plot_dpi (int): Dots per inch (DPI) for the plot.
        res (int): Number of points used for orbit calculations.
        flip (bool): True: +y upwards (as Target's images); False: the y axis turned over.
        mode (str): Coordinate mode of planet_pos ("polar" or "cartesian"); any other value prints a
            message and plots planet_pos as cartesian.
        name (str): Suffix of the saved plot file name.

    Returns:
        None. Saves "oribit"+name+".png" in the current directory.
    """
    # Convert position inputs to Cartesian coordinate system
    if mode == "polar":
        planet_pos = planet_pos[0]*np.cos(planet_pos[1]), planet_pos[0]*np.sin(planet_pos[1])
    elif mode != "cartesian":
        print("Function orbit_plot got unknown mode")

    orbit_x, orbit_y = orbit_position(a, e, pa, inc, res)

    # Create a new figure with a specified DPI
    fig = plt.figure(dpi=plot_dpi)
    axs = plt.gca()

    # Set aspect ratio and axis labels, and turn the y-axis over if requested
    axs.set_aspect('equal', 'box')

    if not flip:
        axs.invert_yaxis()

    axs.set_ylabel('y [arcsec]')
    axs.set_xlabel('x [arcsec]')

    # Plot the orbit in orange with a dashed line
    plt.plot(orbit_x, orbit_y, color='orange', linestyle = '--', alpha = 0.5, label="orbit (as seen on the sky)")

    # Mark the central star with a dark orange 'x'
    plt.scatter([0], [0], color='darkorange', marker='x', label="star")

    # Mark the planet with a black circle
    plt.scatter([planet_pos[0]], [planet_pos[1]], color='black', marker='o', label="planet now")
    axs.legend(fontsize=8)
    axs.set_title("Orbit: a = %.3g arcsec, e = %.2g, PA %.3g deg, inclination %.3g deg" % (a, e, pa, inc), fontsize=9)

    # Save the figure with the specified name
    fig.savefig("oribit"+name+".png", format='png', bbox_inches='tight')

    # Display the plot
    plt.show()
