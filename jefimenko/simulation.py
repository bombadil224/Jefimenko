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

#from .charge_current_integrals import *
from .constants import *
from .classes import *
from .derivative import *
from .extras import timed_range

import sys
from math import isnan as isnan
import time as Time
import types

import julia
from julia.api import Julia

import re

import os

import pdb

script_path = os.path.abspath(__file__)
# fileDir = os.path.dirname(os.path.realpath('__file__'))
script_dir = os.path.split(script_path)[0]

jl = Julia(compiled_modules=False)
j = julia.Julia()
j.include(script_dir + "/charge_current_integrals.jl")

# this will run the acual simulation of a grid
def simulate(grid,
             charge_currents=False,
             induction=False,
             magnetization=False,
             Test=False):

    # first genorat the list of every locaiotn being simulated
    location_list = genorate_location_list(grid)

    # the grid layout is grid.grid[field][time axis][x, y, z axis]
    # first genorate the grid to save everything in
    for t in range(grid.time_size):
        grid.grid['E'].append([])
        grid.grid['H'].append([])
        grid.grid['B'].append([])
        grid.grid['K'].append([])
        for i in range(int(grid.size[0] / grid.delta[0])):
            grid.grid['E'][t].append([])
            grid.grid['H'][t].append([])
            grid.grid['B'][t].append([])
            grid.grid['K'][t].append([])
            for j in range(int(grid.size[1] / grid.delta[1])):
                grid.grid['E'][t][i].append([])
                grid.grid['H'][t][i].append([])
                grid.grid['B'][t][i].append([])
                grid.grid['K'][t][i].append([])
                for k in range(int(grid.size[2] / grid.delta[2])):

                    grid.grid['E'][t][i][j].append(
                        list(grid.constant_E([i * grid.delta[0],
                                              j * grid.delta[1],
                                              k * grid.delta[2]],
                                             t * grid.delta_t)))

                    grid.grid['H'][t][i][j].append(
                        list(grid.constant_H([i * grid.delta[0],
                                              j * grid.delta[1],
                                              k * grid.delta[2]],
                                             t * grid.delta_t)))

                    grid.grid['B'][t][i][j].append(
                        list(grid.constant_H([i * grid.delta[0],
                                              j * grid.delta[1],
                                              k * grid.delta[2]],
                                             t * grid.delta_t)))

                    grid.grid['K'][t][i][j].append([0, 0, 0])

    # start timing the simulaiton
    start = Time.time()

    # this tells if the differentials need to be updated
    # update_diff = True

    ########################################################

    # this is used to flush the print statments befor the loop starts
    # this works in combination with timed_range and may not be needed
    sys.stdout.flush()

    ''' start a loop that call simulate_location with the list of locaitons for
                       every time step being simulated '''
    for time in timed_range(0, grid.time_size, start):

        currents_time_diff(grid)    # this updates the differentials

        E, H, B, K = simulate_location_list(location_list,
                                         time,
                                         grid)

        for (loc, j) in zip(location_list, range(len(location_list))):
            for t in range(grid.time_size):
                for i in range(3):
                    grid.grid['E'][t][loc[0]][loc[1]][loc[2]][i] += E[j][t][i]

                    grid.grid['H'][t][loc[0]][loc[1]][loc[2]][i] += H[j][t][i]

                    grid.grid['B'][t][loc[0]][loc[1]][loc[2]][i] += B[j][t][i]

                    grid.grid['K'][t][loc[0]][loc[1]][loc[2]][i] += K[j][t][i]
    print()


def genorate_location_list(grid):
    location_list = []
    for i in range(int(grid.size[0] / grid.delta[0])):
        for j in range(int(grid.size[1] / grid.delta[1])):
            for k in range(int(grid.size[2] / grid.delta[2])):
                location_list.append([i, j, k])
    return(location_list)


def simulate_location_list(location_list,
                           time,
                           grid,
                           time_end=False):
    if time_end is False:
        time_end = grid.time_size

    E_field = []
    H_field = []
    B_field = []
    K_field = []

    for i in range(len(location_list)):

        E, H, B, K = (j.field_calculator(grid.charges[time],
                                      grid.currents[time],
                                      grid.dipoles[time],
                                      grid,
                                      location_list[i],
                                      time,
                                      time_end))

        E_field.append(E)
        H_field.append(H)
        B_field.append(B)
        K_field.append(K)
    # print(K_field)

    return(E_field, H_field, B_field, K_field)
