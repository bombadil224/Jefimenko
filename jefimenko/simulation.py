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

C_0 = 299792458  # this is the speed of light in meters per secound
K_e = 8.9875517873681764 * 10 ** 9  # this is Coulomb's constant
E_0 = (4 * np.pi * K_e) ** -1  # free space permittivity


def simulate(grid):
    print('simulating grid')
    sys.stdout.flush()
    # this will run the simulation of a grid and charges

    # start a loop that will simulate the system for every point on the grid
    # V is the location on the grid where the E field is being calculated
    for V in np.ndindex(tuple(grid.shape)):
        grid.grid['E'] = electric_charges(grid.grid['E'],
                                          grid.charges,
                                          grid, V)
    print('grid simulated')


def retardation(r, grid):
    return(int(np.rint((r / C_0) / grid.delta_t)))


def electric_charges(E_field, charges, grid, V):
    scale = 1 / (4 * np.pi * E_0)

    # this will calculate the field genorated by the stationary charges
    # loop over the charges in the grid
    for t in range(grid.time_size):
        for charge in charges[t]:
            R = V - charge.location
            r = np.linalg.norm(grid.delta * R)
            R = np.pad(R, (0, 3 - len(R)), 'constant', constant_values=0)

            # now loop over the time axis
            if r != 0:
                # the only case of error should be a T that is to big
                # use a try to find and in this case exit the loop
                # pdb.set_trace()
                try:
                    time = t + retardation(r, grid)
                    E_field[time][V] += (E_field[time][V] + scale * charge.Q *
                                         R / (r**3))
                except IndexError:
                    break
            else:
                E_field[t][V] = E_field[t][V]
    return (E_field)
