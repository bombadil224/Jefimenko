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
from .constants import *

import math
import pdb


class Grid():

    # to initialize a grid just pass it the number of dimensions of the grid
    def __init__(
            self,
            delta=False,  # this is the size of a step in x, y and z
            size=False,  # this is the size of the grid in meters
            time=.1,  # this is the length of a simulation in secounds
            delta_t=.1  # this is the size of a full time step
    ):

        # self.dimension = dimension
        self.dimension = 3

        if delta is False:
            self.delta = [1.0 for i in range(self.dimension)]
        else:
            delta = np.pad(delta, (0, 3 - len(delta)),
                           'constant',
                           constant_values=1)
            print('delta paded')
            self.delta = delta

        self.delta = np.array(self.delta).astype(float)

        if size is False:
            self.size = [10] * self.dimension
            print('grid size diffalting to ' + str(self.size))

        elif len(size) != self.dimension:
            for i in range(len(size), self.dimension):
                size.append(self.delta[i])

            self.size = size
            print(' grid size is ' + str(size))
            print('recommended to give all three measurements')
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

        # this is the conductors that are on the grid
        # notice that these don't change over the course of the simulation
        self.conductors = []

        # the grid layout is grid.grid[field][time axis][x, y, z axis]
        # the grid will hold the E and H fields as 3D vectors for all casses

        self.grid = {}

        self.grid['E'] = np.zeros(
            tuple(np.append(self.shape_t, [3])), dtype='complex')
        self.grid['H'] = np.zeros(
            tuple(np.append(self.shape_t, [3])), dtype='complex')

    def Add_Charge(self,
                   location,
                   Q=0,
                   velocity=False,
                   acceleration=False,
                   time=0,
                   print_all=False,
                   charge_count=False):
        # this will add a charge to the grid notice that a charge is a class

        Q = float(Q)

        location = np.pad(location, (0, 3 - len(location)),
                          'constant', constant_values=0)
        location = np.array(location).astype(float)

        if velocity is False:
            velocity = self.dimension * [0]
        elif len(velocity) != self.dimension:
            velocity = np.pad(velocity, (0, 3 - len(velocity)),
                              'constant', constant_values=0)

        if acceleration is False:
            acceleration = self.dimension * [0]
        elif len(acceleration) != self.dimension:
            acceleration = np.pad(acceleration, (0, 3-len(velocity)),
                                  'constant', constant_values=0)

        # add the charge to every time step from time to time_size
        for i in range(0, time):
            new_charge = Charge(self, location, 0, velocity, acceleration)
            self.charges[i].append(new_charge)

        for i in range(time, self.time_size):
            new_charge = Charge(self, location, Q, velocity, acceleration)
            self.charges[i].append(new_charge)

        if print_all is not False:
            print('charge added with')
            print('location')
            print(self.charges[i][-1].location)
            print('velocity')
            print(self.charges[i][-1].velocity)
            print('acceleration')
            print(self.charges[i][-1].acceleration)
            print('charge')
            print(self.charges[i][-1].Q)
            print('time')
            print(time)
            print('')

        if charge_count is True:
            # this must return the acual index of the new charge
            return(len(self.charges[0]) - 1)

    def modify_charge(self,
                      N,    # what charge you want to modify
                      time,  # at what time you want to modify the charge
                      velocity=False,  # the new velocity at time, time
                      acceleration=False,  # the new acceleration of the charge
                      location=False,  # this changes the charges locaiton
                      Q=False,  # this changes the coulombs of the charge
                      print_all=False):  # set true to print charge info

        # if the location of a charge is changed all informaiton
        # about velocity and acceleration is discarded and overwitten.

        if velocity is not False:
            if len(velocity) != self.dimension:
                velocity = np.pad(velocity, (0, 3 - len(velocity)),
                                  'constant', constant_values=0)
            self.charges[time][N].velocity = np.array(velocity)

        if acceleration is not False:
            if len(acceleration) != self.dimension:
                acceleration = np.pad(acceleration, (0, 3 - len(acceleration)),
                                      'constant', constant_values=0)
            self.charges[time][N].acceleration = np.array(acceleration)

        if Q is not False:
            # notice that Q is not effecting velocity or location
            for t in range(time, len(self.charges)):
                # this loop will likely efect velocity and acceleration
                self.charges[t][N].Q = Q

        if location is not False:
            if len(location) != self.dimension:
                location = np.pad(location, (0, 3 - len(location)),
                                  'constant', constant_values=0)

            for t in range(time, len(self.charges)):
                self.charges[t][N].location = np.array(location)

            temp_location = []
            temp_velocity = []
            temp_acceleration = []

            # since we modified the locaton of the charge update the gradient
            for t in range(len(self.charges)):
                temp_location.append(self.charges[t][N].location)

            temp_velocity = np.gradient(temp_location, axis=0)
            temp_acceleration = np.gradient(temp_velocity, axis=0)

            if print_all is not False:
                print('location')
                print(temp_location)
                print('velocity')
                print(temp_velocity)
                print('acceleration')
                print(temp_acceleration)
                print('charge')
                print(Q)
                print('')

            # notice that the gradients that we just found
            # are the new velocity and acceleration
            for t in range(len(self.charges)):
                self.charges[t][N].velocity = temp_velocity[t]
                self.charges[t][N].acceleration = temp_acceleration[t]

    def Add_Current(self, location,  # location in acual distance
                    direction=[0, 0, 1],  # the direction the current points
                    Amps=1,  # the amps of the current
                    count=False,   # this will return the index of the current
                    print_all=False):  # print all information on the current

        #  this will add a current to the grid notice that a current is a class
        location = np.pad(location, (0, 3 - len(location)),
                          'constant', constant_values=0)

        direction = np.pad(direction, (0, 3 - len(direction)),
                           'constant', constant_values=0)

        if np.linalg.norm(direction) == 0:
            print("Current direction norm can't be zero")
            sys.exit()

        direction = (np.array(direction).astype(float) /
                     np.linalg.norm(direction))

        location = np.array(location)
        # insure that the current is in the middle of a grid squair
        location = (location / self.delta).astype(int) * self.delta

        for i in range(self.time_size):
            new_current = Current(self, location, direction, Amps)
            self.currents[i].append(new_current)
        pass

        if print_all is True:
            print('location = ' + str(location) + 'amps = ' + str(Amps) +
                  ' direction = ' + str(direction))

        if count is True:
            return(len(self.currents[0]) - 1)

    # this will allow you to change the amps and direction of a current
    # at this time you can't change the location
    def Modify_Current(self,
                       N,     # what current you want to modify
                       time=0,  # the time at which to make changes
                       direction=False,  # the new direction
                       amps=False,  # the new amps value
                       print_all=False):

        time_N = int(np.rint(time / self.delta_t))
        if direction is not False:
            direction = np.pad(direction, (0, 3 - len(location)),
                               'constant', constant_values=0)
            print(direction)

            # notice that the direction is eithor nomalized or zero
            if np.linalg.norm(direction) != 0:
                direction = (np.array(direction) /
                             np.linalg.norm(direction))

                print('Current modified')

                self.currents[time_N][N].direction = direction
            elif np.linalg.norm(direction) == 0:
                self.currents[time_N][N].direction = np.array(direction)

        if amps is not False:
            self.currents[time_N][N].amps = float(amps)

        if print_all is True:
            print('Modifying Current')
            print('New Time = ' + str(time_N))
            print('New Amps = ' + str(amps))
            print('New direction = ' + str(direction))
            print('poresent direction = ' +
                  str(self.currents[time_N][N].direction))
            print('')

    def Add_Conductor(self,
                      location,   # the location of the conductor
                      materiel='Copper',  # the materiel of the conductor
                      print_all=False):

        location = np.array(location)
        if print_all is True:
            print('cunductor added')
            print('location = ' + str(location) +
                  ' materiel = ' + str(materiel))

        if len(location) != self.dimension:
            print('conductor grid dimension missmatch')
            sys.exit([2])

        location = np.array(location).astype(float)
        self.conductors.append(Conductor(materiel, location, self, print_all))


class Charge():  # this is used to define a charge on the grid
    def __init__(self, grid, location, Q, velocity, acceleration):

        self.location = np.array(location).astype(float)
        self.Q = Q
        self.velocity = np.array(velocity)
        self.acceleration = np.array(acceleration)


class Current():  # this is used to define a current on the grid
    def __init__(self, grid, location, direction, amps):

        self.location = np.array(location).astype(float)
        self.direction = np.array(direction)
        self.amps = float(amps)
        self.diff_t = 0


# test for approximate equality (for floating point types)
def arreqclose_in_list(myarr, list_arrays):
    return next((True for elem in list_arrays
                 if elem.size == myarr.size and
                 allclose(elem, myarr)), False)


class Conductor():
    def __init__(self, materiel, location, grid, print_all=False):

        #  first add a current at the location of the conductor
        #  this will become the induction
        self.current = grid.Add_Current(location,
                                        direction=[0, 0, 0],
                                        Amps=0,
                                        count=True,
                                        print_all=False)

        # this is the term needed to calculate the induction currents
        # it is the electrokinetic field
        # self.EK_field = np.zeros( tuple(grid.time_size), dtype='complex')
        # self.EK_field = np.zeros( tuple(grid.time_size, 3), dtype='complex')

        self.EK_field = np.zeros(
            tuple(np.append(grid.time_size, [3])), dtype='complex')

        self.location = location
        # the location must be of the form [i, j, k]
        # with the same dimention as the grid

        # all ohm values are given in Ohm Meters at 20C
        if materiel == 'aluminium':
            self.Ohms = 2.8 * 10**-8
        elif materiel == 'Antimony':
            self.Ohm = 3.9 * 10**-7
        elif materiel == 'Bismuth':
            self.Ohm = 1.3 * 10**-6
        elif materiel == 'Brass':
            self.Ohm = 0.75 * 10**-7
        elif materiel == 'Cadmium':
            self.Ohm = 6 * 10**-8
        elif materiel == 'Cobalt':
            self.Ohm = 5.6 * 10**-8
        elif materiel == 'Copper':
            self.Ohm = 1.7 * 10**-8
        elif materiel == 'Gold':
            self.Ohm = 2.4 * 10**-8
        elif materiel == 'Carbon':  # (Graphite)
            self.Ohm = 1 * 10**-5
        elif materiel == 'Germanium':
            self.Ohm = 4.6 * 10**-1
        elif materiel == 'Iron':
            self.Ohm = 1.0 * 10**-7
        elif materiel == 'Lead':
            self.Ohm = 1.9 * 10**-7
        elif materiel == 'Manganin':
            self.Ohm = 4.2 * 10**-7
        elif materiel == 'Nichrome':
            self.Ohm = 1.1 * 10**-6
        elif materiel == 'Nickel':
            self.Ohm = 7 * 10**-8
        elif materiel == 'Palladium':
            self.Ohm = 1.0 * 10**-7
        elif materiel == 'Platinum':
            self.Ohm = 0.98 * 10**-7
        elif materiel == 'Quartz':
            self.Ohm = 7 * 1017
        elif materiel == 'Silicon':
            self.Ohm = 6.4 * 102
        elif materiel == 'Silver':
            self.Ohm = 1.6 * 10**-8
        elif materiel == 'Tantalum':
            self.Ohm = 1.3 * 10**-7
        elif materiel == 'Tin':
            self.Ohm = 1.1 * 10**-7
        elif materiel == 'Tungsten':
            self.Ohm = 4.9 * 10**-8
        elif materiel == 'Zinc':
            self.Ohm = 5.5 * 10**-8
        else:
            self.Ohm = 1.7 * 10**-8
