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
    else t < 0:
        integer_t = int(t)
        if np.isclose(t, integer_t):
            return 1.0
        else:
            return t-integer_t

def planet_position(a, e, pa = 0, inc = 0, t = 0):
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
    
    # Find the projections
    x = r * np.cos(theta-pa)
    y = r * np.sin(theta-pa)*np.cos(inc)

    # Find the planet's location
    planet_r = a * (1 - e**2) / (1 + e * np.cos(target_theta))
    planet_x = planet_r * np.cos(target_theta-pa)
    planet_y = planet_r * np.sin(target_theta-pa)*np.cos(inc)

    return x,y,planet_x, planet_y

def orbit_plot(positions):
    orbit_x, orbit_y, planet_x, planet_y = positions
    
    # Make the unit length of x and y axis show the same pixel in the final picture.
    fig=plt.figure(dpi=300)
    axs = plt.gca()
    axs.set_aspect('equal', 'box')
    
    # Plot
    plt.plot(orbit_x, orbit_y, color='orange', linestyle = '--', alpha = 0.5) # Plot the orbit.
    plt.scatter([0], [0], color='darkorange', marker='x') #Mark the focus.
    plt.scatter([planet_x], [planet_y], color='black', marker='o') #Mark the planet. 
    fig.savefig("oribit.png", format='png', bbox_inches='tight')
    plt.show()