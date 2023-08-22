import numpy as np
import matplotlib.pyplot as plt

def t_period_convert(t):
    if not isinstance(t, (float, int, np.int64, np.float64)) or t == 0:
        print("time/period must be a positive number, auto set to initial position")
        return 1.0
    elif t > 0:
        integer_t = int(t)
        if np.isclose(t, integer_t):
            return 1.0
        else:
            return t-integer_t
    else: #t<0
        integer_t = int(t)
        if np.isclose(t, integer_t):
            return 1.0
        else:
            return t+(1-integer_t)

def planet_position(a, e, pa, inc, t, mode="polar"):
    # Convert inputs from degree to rad
    pa = pa*2*np.pi/360.0
    inc = inc*2*np.pi/360.0

    # Keep the fractional_part of t
    t = t_period_convert(t)
    
    # Calculate the orbit
    theta = np.linspace(0, 2 * np.pi, 1000)
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
    if mode == "polar":
        return planet_r, np.arctan2(np.sin(target_theta+pa)*np.cos(inc), np.cos(target_theta+pa))
    elif mode == "cartesian":
        return planet_r * np.cos(target_theta+pa), planet_r * np.sin(target_theta+pa)*np.cos(inc)
    else:
        print("Function planet_position got unknown mode")    

def orbit_plot(a, e, pa, inc, planet_pos, mode="polar", name='', plot_dpi=300):
    # Convert position inputs to Cartesian coordinate system
    if mode == "polar":
        planet_pos = planet_pos[0]*np.cos(planet_pos[1]), planet_pos[0]*np.sin(planet_pos[1])
    elif mode == "cartesian":
        pass
    else:
        print("Function planet_position got unknown mode")  
    # Convert angle inputs from degree to rad
    pa = pa*2*np.pi/360.0
    inc = inc*2*np.pi/360.0

    # Calculate the orbit
    theta = np.linspace(0, 2 * np.pi, 1000)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))

    # Find the projections
    orbit_x = r * np.cos(theta+pa)
    orbit_y = r * np.sin(theta+pa)*np.cos(inc)

    # Make the unit length of x and y axis show the same pixel in the final picture.
    fig=plt.figure(dpi=plot_dpi)
    axs = plt.gca()
    axs.set_aspect('equal', 'box')
    
    # Plot
    plt.plot(orbit_x, orbit_y, color='orange', linestyle = '--', alpha = 0.5) # Plot the orbit.
    plt.scatter([0], [0], color='darkorange', marker='x') #Mark the focus.
    plt.scatter([planet_pos[0]], [planet_pos[1]], color='black', marker='o') #Mark the planet. 
    fig.savefig("oribit"+name+".png", format='png', bbox_inches='tight')
    plt.show()