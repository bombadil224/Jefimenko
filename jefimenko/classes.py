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
        elif len(size) != self.dimension:
            print('size dimension missmatch grid not genorated')
            sys.exit()
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

    def Add_Charge(self, location, Q=1, velocity=False, acceleration=False):
        # this will add a charge to the grid notice that a charge is a class

        print('location = ' + str(location) + 'Q = ' + str(Q))
        if len(location) != self.dimension:
            print('charge grid dimension missmatch')
            sys.exit([2])

        if velocity is False:
            velocity = self.dimension * [0]
        elif len(veocity) != self.dimension:
            print('velocity dimension missmatch charge not added')
            sys.exit()

        if acceleration is False:
            acceleration = self.dimension * [0]
        elif len(acceleration) != self.dimension:
            print('acceleration dimension missmatch charge not added')
            sys.exit()

        # add the charge to every time step
        for i in range(self.time_size):
            new_charge = Charge(self, location, Q, velocity, acceleration)
            self.charges[i].append(new_charge)

    def modify_charge(self,
                      N,    # what charge you want to modify
                      time,  # at what time you want to modify the charge
                      velocity=False,  # the new velocity at time, time
                      acceleration=False,  # the new acceleration of the charge
                      location=False,  # this changes the charges locaiton
                      print_charge=False):  # set true to print charge info

        # if the location of a charge is changed all informaiton
        # about velocity and acceleration is discarded and overwitten.
        if location is not False:
            if len(location) != self.dimension:
                print('locaton dimension missmatch charge not modified')
                sys.exit()

        if velocity is not False:
            if len(veocity) != self.dimension:
                print('velocity dimension missmatch charge not modified')
                sys.exit()

        if acceleration is not False:
            if len(acceleration) != self.dimension:
                print('acceleration dimension missmatch charge not modified')
                sys.exit()

        if velocity is not False:
            self.charges[time][N].velocity = np.array(velocity)
        if acceleration is not False:
            self.charges[time][N].acceleration = np.array(acceleration)
        if location is not False:
            self.charges[time][N].location = np.array(location)

            temp_location = []
            temp_velocity = []
            temp_acceleration = []

            for t in range(len(self.charges)):
                temp_location.append(self.charges[t][N].location)

            temp_velocity = np.gradient(temp_location, axis=0)
            temp_acceleration = np.gradient(temp_velocity, axis=0)

            if print_charge is not False:
                print('location')
                print(temp_location)
                print('velocity')
                print(temp_velocity)
                print('acceleration')
                print(temp_acceleration)

            for t in range(len(self.charges)):
                self.charges[t][N].velocity = temp_velocity[t]
                self.charges[t][N].acceleration = temp_acceleration[t]

    def Add_Current(self, location,
                    direction=[0, 0, 1],
                    Amps=1):
        # this will add a current to the grid notice that a current is a class
        print('location = ' + str(location) + 'amps = ' + str(Amps))
        if len(location) != self.dimension:
            print('current grid dimension missmatch')
            sys.exit([2])
        if len(direction) != 3:
            print('direction error')
            sys.exit([2])

        for i in range(self.time_size):
            new_current = Current(self, location, direction, Amps)
            self.currents[i].append(new_current)
        pass

    # this will allow you to change the amps and direction of a current
    # at this time you can't change the location
    def Modify_Current(self,
                       N,     # what current you want to modify
                       time=0,  # the time at which to make changes
                       direction=False,  # the new direction
                       amps=False,  # the new amps value
                       Hide=True):

        time_N = int(np.rint(time / self.delta_t))
        if direction is not False:
            self.currents[time_N][N].direction = (direction /
                                                  np.linalg.norm(direction))
        if amps is not False:
            self.currents[time_N][N].amps = float(amps)
        if Hide is False:
            print('Modifying Current')
            print('New Time = ' + str(time_N))
            print('New Amps = ' + str(amps))
            print('New direction = ' + str(direction))
            print()


class Charge():  # this is used to define a charge on the grid
    def __init__(self, grid, location, Q, velocity, acceleration):

        self.location = np.array(location).astype(float)
        self.Q = Q
        self.velocity = np.array(velocity)
        self.acceleration = np.array(acceleration)


class Current():  # this is used to define a current on the grid
    def __init__(self, grid, location, direction, amps):

        self.location = np.array(location).astype(float)
        self.direction = (np.array(direction /
                                   np.linalg.norm(direction)).astype(float))
        self.amps = float(amps)
        self.diff_t = 0
