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

from .constants import *
from .classes import *
from .charge_current_integrals import *
from .derivative import *

import numpy as np
# this is a attempt to add polarization to the model

''' see examples 8-7.1 and 8-7.6 of electricity and magnetisem '''


def normals(Permittivity, shape, delta):
    # this uses the direction of the greatest directional dirivetiv as the
    # normal to the Permittivity

    pol_normals = []  # this is what needs found
    # next is the map of all permitiv regions
    permittivity_map = np.zeros(tuple(shape), dtype='complex')

    for i in Permittivity:
        permittivity_map[i.index] = i.R_Permittivity

    # find the gradient of the permittivity
    per_G = space_grad(permittivity_map, delta)

    for i in np.ndindex(tuple(shape)):
        if (np.linalg.norm(per_G[i]) != 0):
            # find the polarization normals
            pol_normals.append(permittivity_normals(i,
                                                    per_G[i] /
                                                    np.linalg.norm(per_G[i])))

        # else:
        #    pol_normals.append(permittivity_normals([0,0,0], i))
            # pol_normals.append(normal)

    return(pol_normals)


def polarization_effect(E_field,
                        field_P,
                        grid_P_grad,
                        Permittivity_normals,
                        location,
                        time_0,
                        delta,
                        grid):

    # this uses equation 8-7.12

    # since we have the gradient of the polarization
    # and the normals of the Permittivity we are ready to do the calculaiton

    # define the varibels that will be the output of the polarization field
    field_E = np.zeros((grid.time_size, 3))
    field_E1 = np.zeros((grid.time_size, 3))
    field_E2 = np.zeros((grid.time_size, 3))

    # first the volum integreal

    # find the field resulating form the polarization field at time "time"

    # this is an integral over the polorization regions
    for P_location in np.ndindex(tuple(grid.shape)):

        ''' find the distance from the polorization location
          to the locaiton of intrest '''
        R = location - (P_location * delta)
        r = np.linalg.norm(R)

        if r != 0:
            time = time_0 + retardation(r, grid.delta_t)

            if time < grid.time_size:
                field_E1[time] = (field_E1[time] +
                                  Polarization_F_C_V(grid_P_grad[P_location],
                                                     r, R, delta))
                pass

    # now the surfface integreal
    # integration is over the  normals to the Permittivity regions
    for normal in Permittivity_normals:

        R = location - (normal.index * delta)
        r = np.linalg.norm(R)

        if r != 0:
            time = time_0 + retardation(r, grid.delta_t)

            if time < grid.time_size:

                field_E2[time] = (field_E2[time] +
                                  Polarization_F_C_S(field_P[normal.index],
                                                     normal.vector,
                                                     r, R, delta))

    field_E = field_E1 + field_E2 + field_E
    return(field_E)


# Polarization_field_charge_V
def Polarization_F_C_V(grid_P_grad, r, R, delta):
    # see equation 8-7.11 of electromagnetic theory, jefimenko
    scale = 1 / (4 * np.pi * E_0)

    field_E = - scale * np.dot(grid_P_grad, [1, 1, 1]) * R / r**2
    return(field_E)


# Polarization_field_charge_S
def Polarization_F_C_S(field_P,
                       Permittivity_normal_vector,
                       r,
                       R,
                       delta):
    scale = 1 / (4 * np.pi * E_0)

    field_E = (- scale * np.dot(Permittivity_normal_vector, field_P) *
               np.dot(R, delta) * R / r**2)
    return(field_E)


def polarization_field(E_field, P_field, E_r, shape, grid):

    for i in range(len(E_r)):

        # page 455 handbook of physics
        for location in np.ndindex(tuple(shape)):

            ''' Note that if one wishes to add the ability to use
                Anisotropic dielectrics this is the place to do it '''
            P_field[location] = (E_0 * E_field[location] *
                                 (grid.get_Permittivity(location) - 1))
    return(P_field)
