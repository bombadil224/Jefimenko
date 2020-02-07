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
import sys
from math import isnan as isnan
from scipy.spatial import distance
import time as Time
import types

from .simulation import simulate_location_list

#jl = Julia(compiled_modules=False)
#j = julia.Julia()
#j.include(script_dir + "/plasma.jl")

def plasma_simulation(grid, test=False):

    start = Time.time()
    # do any poreporations

    sys.stdout.flush()

    for time in timed_range(0, grid.time_size - 1,
                            start, print_time=not(test)):

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
        '''this should be changed to include the posibilety of the charges
        contracting to a point that is outside of the starting space'''

        ''' use field calculator to find the field strength
            at the location of every charge in the grid
            notice that this is not a full grid simulaiton
            but a location by location simulation '''

        E_field = []
        H_field = []
        B_field = []

        for location in location_list:
            E_field.append(list(grid.constant_E(location,
                                                time * grid.delta_t)))

            H_field.append(list(grid.constant_H(location,
                                                time * grid.delta_t)))

            B_field.append(list(grid.constant_H(location,
                                                time * grid.delta_t)))

        E, H, B, K = simulate_location_list(location_list,
                                         time,
                                         grid,
                                         time_end=time+1)

        for i in range(len(E)):
            E_field[i] = [E_F + E_S for (E_F, E_S) in
                          zip(E_field[i], E[i][time])]
            H_field[i] = [H_F + H_S for (H_F, H_S) in
                          zip(H_field[i], H[i][time])]
            B_field[i] = [B_F + B_S for (B_F, B_S) in
                          zip(B_field[i], B[i][time])]

        if test is True:
            print('charge list is')
            print(location_list)
            print('E field is')
            print(E_field)
            print('H field is')
            print(H_field)
            print()
            print()
        # calculate the velocity and location of charges at time + dt
        if time + 1 < grid.time_size:
            boris_method(location_list,
                         E_field,
                         H_field,
                         B_field,
                         time,
                         grid,
                         test)

        if test is True:
            print('location of partical 0')
            print([grid.charges[t][0].location for t in range(grid.time_size)])

        ''' calculate the effects on E and B from time_i
            giving E and B at time t_i + dt
            this is the first thing in the loop '''
    # now that everything is done figyer out how long it touck
    end = Time.time()
    print('grid simulated in ' + str(end - start) + ' seconds         ')
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

        # update the position of the charge
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


def boris_method(location_list,
                 E_field,
                 H_field,
                 B_field,
                 time,
                 grid,
                 test=False):

    charges = grid.charges

    for n in range(len(grid.charges[time])):

        # pdb.set_trace()

        V_minus = [(V + charges[time][n].Q * E /
                    charges[time][n].mass * grid.delta_t / 2)
                   for (V, E) in zip(charges[time][n].velocity, E_field[n])]

        T_vector = [charges[time][n].Q * B / charges[time][n].mass
                    * grid.delta_t / 2 for B in B_field[n]]

        S_vector = [2 * T / (1 + dot(T_vector, T_vector)) for T in T_vector]

        V_prime = [V + V_M for (V, V_M) in zip(V_minus, cross(V_minus, T_vector))]

        V_plus = [V + S for (V, S) in zip(V_minus, cross(V_prime, S_vector))]

        Velocity = [V + charges[time][n].Q *
                    E / charges[time][n].mass *
                    grid.delta_t / 2 for (V, E) in zip(V_plus, E_field[n])]

        grid.charges[time + 1][n].velocity = Velocity
        grid.charges[time + 1][n].location = [V * grid.delta_t + L
                                              for (V, L) in
                                              zip(charges[time + 1][n].velocity,
                                                  charges[time][n].location)]

def dot(A, B):
    return((A[0] * B[0] + A[1] * B[1] + A[2] * B[2]))

def cross(A, B):
    return(([A[1] * B[2] - A[2] * B[1],
             A[2] * B[0] - A[0] * B[2],
             A[0] * B[1] - A[1] * B[0]]))
