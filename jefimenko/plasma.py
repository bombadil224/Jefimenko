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
from .simulation import *

import re
import numpy as np
import sys
from math import isnan as isnan
from scipy.spatial import distance
import time as Time
import types


import multiprocessing
# from multiprocessing import pool as Pool
# from multiprocessing import pool.ThreadPool as Pool
from functools import partial
import os
import pdb


def plasma_simulation(grid, test=False):

    start = Time.time()
    # do any poreporations

    # calculate the E and B boundary conditions of the simulation

    # calculate the particals x,y,z and v boundary conditions at t_0

    sys.stdout.flush()
    for time in timed_range(0, grid.time_size - 1, print_time=True):
        ''' 1. calculate the magnetic and electric
            field at time t_0 using jefimenko
            notice that do to the nature of jefimenko
            this gives the efective field
            at time t_i + dt for all particals '''
        # location_list is a list of all locations of particals
        location_list = []
        for charge in grid.charges[time]:
            location_list.append(charge.location)

        # find the maximum distance between all locations particals
        '''this should be changed ot include the posibilety of the charges
        contracting to a point that is outside of the starting space'''
        test_max = np.max(distance.cdist(location_list,
                                         location_list,
                                         'euclidean'))
        if grid.max_charge_distance < test_max:
            grid.max_charge_distance = test_max

        if test is True:
            print('charge list is')
            print(location_list)
            print()
            print()

        ''' use field calculator to find the field strength
            at the location of every charge in the grid
            notice that this is not a full grid simulaiton
            but a location by location simulation '''
        E_field, H_field, B_field = simulate_location(location_list,
                                                      time,
                                                      grid,
                                                      test)

        # calculate the velocity and location of charges at time + dt
        move_charges(location_list,
                     E_field,
                     H_field,
                     B_field,
                     time,
                     grid,
                     test)

        ''' calculate the effects on E and B from time_i
            giving E and B at time t_i + dt
            this is the first thing in the loop '''
    # now that everything is done figyer out how long it touck
    end = Time.time()
    print('grid simulated in ' + str(end - start) + ' seconds')
    print()
    pass


def move_charges(location_list,
                 E_field,
                 H_field,
                 B_field,
                 time,
                 grid,
                 test=False):

    for n in range(len(grid.charges[0])):

        location = str(location_list[n])

        if (E_field[location][time] != [0, 0, 0]).any():
            force_E = (np.array(grid.charges[time][n].Q *
                                E_field[location][time]))
        else:
            force_E = np.array([0, 0, 0])

        if (B_field[location][time] != [0, 0, 0]).any():
            force_B = (np.array(grid.charges[time][n].Q *
                                np.cross(grid.charges[time][n].velocity,
                                         B_field[location][time])))
        else:
            force_B = np.array([0, 0, 0])

        force = force_E + force_B

        # calculate the acceloration of the charges
        acceleration = force / grid.charges[time][n].mass
        grid.charges[time][n].acceleration = acceleration.astype(float)

        # update the velocity of the charge
        velocity = grid.charges[time][n].velocity + acceleration * grid.delta_t
        grid.charges[time+1][n].velocity = velocity.astype(float)

        # update the position of the charge
        charge_location = (grid.charges[time][n].location +
                           (grid.charges[time][n].velocity +
                            grid.charges[time+1][n].velocity) / 2 *
                           grid.delta_t)
        grid.charges[time+1][n].location = charge_location.astype(float)

        if test is True:

            print('location is ' +
                  str(tuple((grid.charges[time][n].location))))
            print('electric field is ' +
                  str(E_field[str(location_list[n])][time]))
            print('magnetic field is ' +
                  str(B_field[str(location_list[n])][time]))
            print('acceleration is ' + str(acceleration))
            print('force is ' + str(force))
            print('the time is ' + str(time))
            print('the charge is ' + str(n))
            print('old velocity is ' +
                  str(grid.charges[time][n].velocity))
            print('new velocity is ' +
                  str(grid.charges[time+1][n].velocity))
            print(' location is ' + str(location))
            print()


def simulate_location(location_list,  # the list of locations to include
                      end_time,       # the time to find the field at
                      grid,
                      test):

    E_field = {}
    B_field = {}
    H_field = {}

    # location_list is a list of all locations of particals

    for location in location_list:
        (E_field[str(location)],
         H_field[str(location)],
         B_field[str(location)]) = location_process(end_time,
                                                    grid, location)

    # if __name__ == '__main__':

    # pool = multiprocessing.Pool(processes = 9)
    # pool = Pool(4)
    # pool = Pool.ThreadPool(4)

    # func = partial(location_process, end_time, grid)

        # pdb.set_trace()

#    a = pool.map(func, location_list)
#    pool.close()
#    pool.join()
#
#    for i in range(len(location_list)):
        # E_field[str(location_list[i])] = []
        # B_field[str(location_list[i])] = []
        # H_field[str(location_list[i])] = []

        # E_field[str(location_list[i])] = a[i][0]
        # B_field[str(location_list[i])] = a[i][1]
        # H_field[str(location_list[i])] = a[i][2]
#
#        #pdb.set_trace()
#
#    print(a)

    if test is True:
        print()
        print('location test')
        print('E = ' + str(E_field))
        print('H = ' + str(H_field))
        print('B = ' + str(B_field))
        print()
    return(E_field, H_field, B_field)


def location_process(end_time, grid, location):

    E_field = np.zeros(
            tuple(np.append(grid.time_size, [3])), dtype='complex')
    H_field = np.zeros(
            tuple(np.append(grid.time_size, [3])), dtype='complex')
    B_field = np.zeros(
            tuple(np.append(grid.time_size, [3])), dtype='complex')

    loc_key = str(location)

    if grid.constant_E != [0, 0, 0]:

        if type(grid.constant_E) == types.FunctionType:
            loc_key_F = " ".join(loc_key.split())
            loc_key_F = loc_key_F[1:-1].strip()

            loc_key_F = [float(s) for s in loc_key_F.split(' ')]

            for i in range(len(E_field)):
                E_field[i] = grid.constant_E(loc_key_F,
                                             i * grid.delta)
        else:
            E_field = grid.constant_E

    if grid.constant_H != [0, 0, 0]:

        if isinstance(types.FunctionType, type(grid.constant_H)):
            loc_key_F = " ".join(loc_key.split())
            loc_key_F = loc_key_F[1:-1].strip()

            loc_key_F = [float(s)
                         for s in loc_key_F.split(' ')]
            for i in range(len(E_field)):
                H_field[i] = grid.constant_H(loc_key_F,
                                             i * grid.delta)
        else:
            H_field = grid.constant_H
            B_field = grid.constant_H

    for time in range(end_time - int(np.ceil(grid.max_charge_distance / C_0)),
                      end_time + 1):

        E_field, H_field, B_field = (
            field_calculator(E_field,
                             H_field,
                             B_field,
                             grid.charges[time],
                             grid.currents[time],
                             grid.dipoles[time],
                             grid,
                             location,
                             time))

    return([E_field, H_field, B_field])
