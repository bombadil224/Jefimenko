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

# this file holds the moduals for calculationg the jefimenko equations
from .derivative import *
from .constants import *

import numpy as np
import sys


def retardation(r, grid, c=C_0):
    # this function calculates the retardation used in all other functions
    return(int(np.rint((r / c) / grid.delta_t)))


def field_calculator(E_field, H_field, B_field,
                     charges, currents, dipoles,
                     grid, location, time_0):
    # this function maneges the calculaiton of the E and H fields

    # first we need some grids to hold the added field strenth
    field_E = np.zeros((grid.time_size, 3))
    field_H = np.zeros((grid.time_size, 3))
    field_B = np.zeros((grid.time_size, 3))
    time = time_0

    # find the field resulating form the currents at time "time"
    for current in (currents):
        # find the distance from the current to the locaiton of intrest
        R = location - current.location
        r = np.linalg.norm(R)

        if r != 0:
            # calculate the time at which to add the new field quantity
            time = time_0 + retardation(r, grid)

            # if time is less then the length of the simulaiton
            # calculate the added field
            if time < grid.time_size:
                field_E[time] = (field_E[time] +
                                 Electric_field_currents(current,
                                                         location,
                                                         r,
                                                         R,
                                                         grid))
                field_H[time], field_B[time] = (field_H[time] +
                                 Magnetic_field_currents(current,
                                                         r,
                                                         R,
                                                         grid))
                

    time = time_0
    # now repet the above for charges
    for charge in charges:
        # calcualte the distance to the charge
        R = location - charge.location
        r = np.linalg.norm(R)

        # find the field resulating from the charge at time "time"
        if r != 0:
            time = time_0 + retardation(r, grid)

            # if time is less then the length of the
            # simulaiton calculate the added field
            if time < grid.time_size:
                field_E[time] = (field_E[time] +
                                 Electric_field_charges(charge,
                                                        location,
                                                        r,
                                                        R,
                                                        grid))
                field_H[time], field_B[time] = (field_H[time] +
                                 Magnetic_field_charges(charge,
                                                        location,
                                                        r,
                                                        R,
                                                        grid))

    # now repete the above for dipoles
#    for dipole in dipoles:
#        # calcualte the distance to the dipole
#        R = location - dipole.location
#        r = np.linalg.norm(R)
#        # find the field resulating from the dipole at time "time"
#        if r != 0:
#            time = time_0 + retardation(r, grid)
#
#            # if time is less then the length of the
#            # simulaiton calculate the added field
#            if time < grid.time_size:
#                field_E[time] = (field_E[time] +
#                                 Electric_field_dipole(dipole, r, R, grid))

    # E_field_Permittivity(grid)
    # H_field_Permeability(grid)

    # add everything up
    E_field = E_field + field_E
    H_field = H_field + field_H
    B_field = B_field + field_B

    return(E_field, H_field, B_field)


def Electric_field_currents(current, location, r, R, grid):
    # find the E field do to currents on the grid
    # see equation 2-2.12 of
    # electromagnetic retardation and theory of relativity

    scale = 1 / (4 * np.pi * E_0 *
                 # C_0 * grid.get_Permittivity(location))
                 C_0 * grid.get_Permittivity(current.location))
    # first find the terms dependent on current.amp
    # this is the electrostatic field

    electrostatic = (- scale *
                     current.amps *
                     current.direction *
                     2 / r**2)

    electrostatic = (electrostatic + scale *
                     2 * R * np.dot(current.amps *
                                    current.direction,
                                    R) * 2 / (r**4))

    # now find the terms dependent on current.dif_t
    # this is the electrokinetic field

    electrokinetic = (scale *
                      R * np.dot(current.diff_t, R) /
                      (r**3 * C_0))

    electrokinetic = (electrokinetic - scale *
                      current.diff_t / (r * C_0))

    # find the total E field
    E_field = (electrostatic + electrokinetic)
    return(E_field)


def Magnetic_field_currents(current, r, R, grid):
    # this function finds the H field do to currents on the field
    # see equation 2-2.5 of
    # electromagnetic retardation and theory of relativity

    # calculate the H field resulting from the currents
    scale = 1 / (4 * np.pi)

    # first find the terms do to steady currents
    H_field_amps = (scale * current.amps *
                    np.cross(current.direction, R) * 2 / (r**2))

    # now find the term dependent on current.diff_t
    # this is the electrokinetic field
    H_field_diff_t = (scale * np.cross(current.diff_t, R) / (r**2 * C_0))

    # find the total H field
    H_field = H_field_amps + H_field_diff_t

    B_field = U_0 * H_field
    return(H_field, B_field)


def Electric_field_charges(charge, location, r, R, grid):
    # this function finds the E field resulting from charges
    # see equation 4-4.34 of
    # electromagnetic retardation and theory of relativity

    # scale = 1 / (4 * np.pi * E_0 * grid.get_Permittivity(location))
    scale = 1 / (4 * np.pi * E_0 * grid.get_Permittivity(charge.location))

    velocity = np.pad(charge.velocity, (0, 3 - len(charge.velocity)),
                      'constant', constant_values=0)

    acceleration = (np.pad(charge.acceleration,
                           (0, 3 - len(charge.acceleration)),
                           'constant', constant_values=0))

    common_factor = charge.Q / (r**3 *
                                (1 - np.dot(R, (velocity /
                                            (r * C_0))))**3)

    # find the common terms of the remaining parts
    radius_factor = (R - r * velocity / C_0)

    # next find the factor resalting from velocity only
    velocity_term = (scale * common_factor * radius_factor *
                     (1 - velocity ** 2 / C_0 ** 2))

    # now find the factor that is dependent on acceleration
    acceleration_factor = np.cross(radius_factor,
                                   acceleration / C_0**2)

    # now find the full term that is dependent on acceleration
    acceleration_term = (scale * common_factor
                         * np.cross(R, acceleration_factor))

    # find the total E field
    E_field = (velocity_term + acceleration_term)
    return(E_field)


def Magnetic_field_charges(charge, location, r, R, grid):
    # this function finds the H field resulting from charges
    # see equation 4-5.10 of
    # electromagnetic retardation and theory of relativity

    scale = 1 / (4 * np.pi)

    velocity = np.pad(charge.velocity, (0, 3 - len(charge.velocity)),
                      'constant', constant_values=0)

    acceleration = np.pad(charge.acceleration, (0, 3 -
                                                len(charge.acceleration)),
                          'constant', constant_values=0)

    # first find the most used term in the remaning equations
    most_used_factor = 1/(r * (1 - np.dot(R, (velocity / (r * C_0)))))

    # find the common term
    common_factor = charge.Q * scale * most_used_factor**2

    # next find the velocity factor
    velocity_factor = (1 - np.linalg.norm(velocity)**2 / C_0**2
                       + np.dot(R, acceleration) / C_0**2)

    # now find the full velocity term
    velocity_term = (common_factor * most_used_factor
                     * velocity_factor * np.cross(velocity, R))

    # now find the term that is dependent on acceleration
    # this is the electrokinetic field
    acceleration_term = common_factor * np.cross(acceleration, R) / C_0

    # find the total H field
    H_field = (velocity_term + acceleration_term)

    B_field = U_0 * H_field
    return(H_field, B_field)


def Electric_field_dipole(dipole, location, r, R, grid):
    # this is the equation for a electric dipole from page 449
    # of handbook of phisics
    # also see equation 5-4.13 from electricity and magnetism by jefimenko
    scale = 1 / (4 * np.pi * E_0)

    # first find the term requiering a dot product
    first_term = (3 * np.dot(dipole.dipole_moment, R) * R) / r**5

    # now find the remaining term
    last_term = dipole.dipole_moment / r**3

    E_field = scale * (first_term - last_term)
    print('E_field test ' + str(E_field))
    print(dipole.dipole_moment)
    return(E_field)