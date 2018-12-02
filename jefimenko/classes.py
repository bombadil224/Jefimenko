#!/usr/bin/env python
""" jefimenko - An EM simulator based on the Jefimenko equations

This file is part of Jefimenko.

Jefimenko is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import sys


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

        # this will be the charges
        self.charges = [[] for i in range(self.time_size)]
        # this will be the currents
        self.currents = [[] for i in range(self.time_size)]

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
            self.charges[i].append(charge(self, location, Q))

    def Add_Current(self, location, direction):
        # this will add a current to the grid notice that a current is a class
        # not yet implemented
        pass


class charge():  # this is used to define a charge on the grid
    def __init__(self, grid, location, Q):

        self.location = np.array(location).astype(float)
        self.Q = Q
