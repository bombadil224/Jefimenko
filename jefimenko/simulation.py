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

from .derivative import *
import numpy as np
import sys

import pdb

C_0 = 299792458  # this is the speed of light in meters per secound
K_e = 8.9875517873681764 * 10 ** 9  # this is Coulomb's constant
E_0 = (4 * np.pi * K_e) ** -1  # free space permittivity
U_0 = (C_0 ** 2 * E_0) ** -1   # free space permeability


def simulate(grid):
    print('simulating grid')
    # this tells weather the differentials need to be updated
    update_diff = True
    # this is used to flush the print statments befor the loop starts
    sys.stdout.flush()
    # this will run the simulation of a grid and charges

    # start a loop that will simulate the system for every point on the grid
    # location is the location on the grid
    # where the E field is being calculated
    for time in range(grid.time_size):
        # update_diff must be tested at the start of each time loop incase
        # something is changed. this will be needed when moving charges
        # are introduced
        if update_diff is True:
            currents_time_diff(grid)
            update_diff = False

        for location in np.ndindex(tuple(grid.shape)):
            # first find the part of E that is dependent on charges
            grid.grid['E'] = electric_charges(grid.grid['E'],
                                              grid.charges[time],
                                              grid,
                                              location,
                                              time)

            # next find the part of E that is dependent on grid.currents
            grid.grid['E'] = electric_currents(grid.grid['E'],
                                               grid.currents[time],
                                               grid,
                                               location,
                                               time)

            # now find H
            # first find the part of H that is dependent on grid.currents
            # notice that this will cumpleatly
            # solve for H in the case of the jefimenko equations
            # pdb.set_trace()
            grid.grid['H'] = currents(grid.grid['H'],
                                      grid.currents[time],
                                      grid,
                                      location,
                                      time)
    print('grid simulated')


def retardation(r, grid):
    return(int(np.rint((r / C_0) / grid.delta_t)))


def electric_charges(E_field, charges, grid, location, time_0):
    # note that charges are considerd to be confined to a singel cell
    scale = 1 / (4 * np.pi * E_0)

    # this will calculate the field genorated by the stationary charges
    # loop over the charges in the grid
    for charge in charges:
        R = location - charge.location
        R = grid.delta * R
        r = np.linalg.norm(R)
        R = np.pad(R, (0, 3 - len(R)), 'constant', constant_values=0)

        if r != 0:
            # the only case of error should be a time that is to big
            # use a try to find and in this case exit the loop
            try:
                # first the terms dependent on charge.Q
                time = time_0 + retardation(r, grid)
                E_field[time][location] = (E_field[time][location] + scale *
                                           charge.Q *
                                           R / (r**3))
            except IndexError:
                break
        else:
            E_field[time_0][location] = E_field[time_0][location]
    return (E_field)


def electric_currents(E_field, currents, grid, location, time_0):
    scale = np.linalg.norm(grid.delta) / (4 * np.pi * E_0 * C_0)

    # this will calculate the field genorated by the constant currents
    # loop over the currents in the grid
    for current in currents:
        R = location - current.location
        R = grid.delta * R
        r = np.linalg.norm(R)
        R = np.pad(R, (0, 3 - len(R)), 'constant', constant_values=0)

        if r != 0:
            try:
                # first find the terms dependent on current.amp
                time = time_0 + retardation(r, grid)
                E_field[time][location] = (E_field[time][location] - scale *
                                           current.amps * current.direction *
                                           2 / r)

                E_field[time][location] = (E_field[time][location] + scale *
                                           2 * R * np.dot(current.amps *
                                           current.direction, R) * 2 /
                                           (r**3))

                # now find the terms dependent on current.dif
                E_field[time][location] = (E_field[time][location] + scale *
                                           R * np.dot(current.diff_t, R) *
                                           np.log(r) / (r**2 * C_0))

                E_field[time][location] = (E_field[time][location] - scale *
                                           current.diff_t * np.log(r) / (C_0))
            except IndexError:
                break
        else:
            E_field[time_0][location] = E_field[time_0][location]
    return(E_field)


def currents(H_field, currents, grid, location, time_0):
    # notice that this is for the H field to finb B maltiply by U_0
    scale = U_0 * np.linalg.norm(grid.delta) / (4 * np.pi)

    # this will calculate the field genorated by the currents
    # loop over the currents in the grid
    for current in currents:
        R = location - current.location
        R = grid.delta * R
        r = np.linalg.norm(R)
        R = np.pad(R, (0, 3 - len(R)), 'constant', constant_values=0)

        if r != 0:
            try:
                # pdb.set_trace()
                # first find the term dependent on current.amp
                time = time_0 + retardation(r, grid)
                H_field[time][location] = (H_field[time][location] + scale *
                                           current.amps *
                                           np.cross(current.direction, R) *
                                           2 / (r**2))

                # now find the term dependent on current.diff_t
                H_field[time][location] = (H_field[time][location] + scale *
                                           np.cross(current.diff_t, R) *
                                           np.log(r) / (r * C_0))
            except IndexError:
                break
        else:
            H_field[time_0][location] = H_field[time_0][location]
    # pdb.set_trace()
    return(H_field)
