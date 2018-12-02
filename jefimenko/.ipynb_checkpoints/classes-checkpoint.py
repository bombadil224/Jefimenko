import numpy as np
import sys
import tqdm as tqdm

class Grid():

    # to initialize a grid just pass it the number of dimensions of the grid
    def __init__(
            self,
            dimension,
            delta=False,  # this is the size of a step in x, y and z
            size=False,  # this is the size of the grid in meters
            time=.1,  # this is the length of a simulation in secounds
            delta_t=.1  # this is the size of a full time step
    ):

        self.dimension = dimension

        if delta is False:
            self.delta = [1.0 for i in range(self.dimension)]
        else:
            self.delta = delta
        self.delta = np.array(self.delta)

        if size is False:
            self.size = [10] * self.dimension
        else:
            self.size = size
        self.size = np.array(self.size)

        # this is the amount of time the simulation will simulate
        self.time = time
        self.delta_t = delta_t  # this is the size of a time step in seconds

        # this is the size of the grids time axis
        self.time_size = np.array(self.time / self.delta_t).astype(int)
        # next the shape of the grids spacial axis
        self.shape = np.array((self.size / self.delta).astype(int))
        # now the total grid
        self.shape_t = np.append(self.time_size, self.shape)

        self.charges = [] * self.time_size  # this will be the charges
        self.currents = []  # this will be the currents

        # the grid layout is grid.grid[field][time axis][x, y, z axis]
        # the grid will hold the E and H fields as 3D vectors for all casses

        self.grid = {}

        self.grid['E'] = np.zeros(
            tuple(np.append(self.shape_t, [3])), dtype='complex')
        self.grid['H'] = np.zeros(
            tuple(np.append(self.shape_t, [3])), dtype='complex')

    def Add_Charge(self, location, Q=1):
        # this will add a charge to the grid notice that a charge is a class

        print('location = ' + str(location) + 'Q = ' + str(Q))
        if len(location) != self.dimension:
            print('charge grid dimension missmatch')
            sys.exit([2])

        new_charge = charge(self, location, Q)
        for i in range(self.time_size):
            self.charges.append(charge(self, location, Q))

    def Add_Current(self, location, direction):
        # this will add a current to the grid notice that a current is a class
        # not yet implemented
        pass

class charge():  # this is used to define a charge on the grid
    def __init__(self, grid, location, Q):

        self.location = np.array(location).astype(float)
        self.Q = Q