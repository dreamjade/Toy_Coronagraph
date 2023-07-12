class Target(object):
    """
    A circular symmetric target, which could be a dust ring, a debris disk, or ice remains
    """
    def __init__(self):
        pass


class FinalImage(Target):
    """
    Simulate the image of a target go through the coronagraph

    Args:
        height (float): a height in meters
        dt (float): timestep of the simulation in seconds
    """
    def __init__(self, height, dt=0.1):
        """
        Function that is run to initialize the class.

        The input `self` is required for functions that belong to an object,
        meaning that you want to make the function access and/or depend on the 
        attributes of the object (e.g., self.time, and self.velocity below)
        """
        # let's initalize it's parent class (empty for now because it is a blank class)
        super().__init__()

        # note that we are not using the astropy.units class here as we haven't talked about it yet! But it could be useful!
        self.height = height # current height [meters]
        self.velocity = 0 # current velocity [meters/second]
        self.time = 0 # time elapsed [seconds]
        self.dt = dt # timestep of the simulation [seconds]
        self.g = -9.8 # gravitational acceleration (Don't change) [meters/second^2]


    def get_num_steps_run(self):
        """
        Function that returns the number of timesteps that have run by comparing self.time with self.dt

        Returns:
            num_steps (int): number of time steps already completed in the simulation
        """
        num_steps = int(self.time / self.dt)
        return num_steps

    ##### Activity ######
    """
    Add functionality to advance the particle's height by one time step at a time. (hint: implement the function below).
    Then use this code to calculate how long it takes for the particle to fall down from a height of 10 meters.

    Some useful equations for how to calculate the particle's new state at the next time step.
    Pseudo code below:
    acceleration = g
    new_velocity = current_velocity + acceleration * dt
    new_height = current_height + new_velocity * dt

    Add inputs and outputs. 
    """
    def simulate_timestep(self):
        """
        Advance the simulation time by a single timestep (self.dt). 
        Update the simulation with the new time, height, and velocity

        Returns:
            height (float): the current height in meters
        """
        ## HINT: Modify code here
        return 0. # currently does nothing