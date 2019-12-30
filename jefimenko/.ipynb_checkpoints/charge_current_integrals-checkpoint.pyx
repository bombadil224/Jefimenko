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

#from .derivative import *

#from .constants import *
#from .py_constants import *
# from .constants import *

import numpy as np

cdef double pi = 3.141592653589793
cdef double C_0 = 299792458  # this is the speed of light in meters per secound
cdef double K_e = 8.9875517873681764 * 10 ** 9  # this is Coulomb's constant
cdef double E_0 = (4 * pi * K_e) ** -1  # free space permittivity
cdef double U_0 = (C_0 ** 2 * E_0) ** -1   # free space permeability


cdef retardation(double r, double delta_t, double c=C_0):

    cdef int ret_val = 0
    ret_val = int(round((r / c) / delta_t))

    return(ret_val)


cpdef field_calculator(E_field,
                       H_field,
                       B_field,
                       charges,
                       currents,
                       dipoles,
                       grid,
                       location,
                       int time_0,
                       int time_size):
    
    '''# add time_size to call'''
    
    # this function maneges the calculaiton of the E and H fields
    cdef double r
    cdef double delta_t = grid.delta_t
    cdef int time = time_0

    # first we need some grids to hold the added field strenth
    field_E = np.zeros((time_size, 3))
    field_H = np.zeros((time_size, 3))
    field_B = np.zeros((time_size, 3))

    # find the field resulating form the currents at time "time"
    for current in (currents):
        # find the distance from the current to the locaiton of intrest
        R = location - current.location
        # r = np.linalg.norm(R)
        r = dot(R, R)**.5

        if r != 0:
            # calculate the time at which to add the new field quantity
            time = time_0 + retardation(r, delta_t)

            # if time is less then the length of the simulaiton
            # calculate the added field
            Permittivity = grid.get_Permittivity
            if time < time_size:
                field_E[time] = (field_E[time] +
                                 Electric_field_currents(current.diff_t,
                                                         list(current.location),
                                                         current.amps,
                                                         list(current.direction),
                                                         r,
                                                         list(R), [1,1,1]))
                                                         #list(grid.get_Permittivity(location)))

                field_H[time], field_B[time] = (field_H[time] +
                                                Magnetic_field_currents(
                                                    current,
                                                    r,
                                                    R,
                                                    grid))

    time = time_0
    # now repet the above for charges
    for charge in charges:
        # calcualte the distance to the charge
        R = location - charge.location
        # r = np.linalg.norm(R)
        r = dot(R, R)**.5

        # find the field resulating from the charge at time "time"
        if r != 0:
            time = time_0 + retardation(r, delta_t)

            # if time is less then the length of the
            # simulaiton calculate the added field
            if time < time_size:
                #print()
                #print()
                #print('Permittivity(charge.location) test ')
                #print(grid.get_Permittivity(charge.location))
                field_E[time] = (field_E[time] +
                                  Electric_field_charges(charge,
                                                         charge.velocity,
                                                         charge.acceleration,
                                                         charge.Q,
                                                         location,
                                                         r,
                                                         R,
                                                         grid,
                                                         # [1,1,1]))
                                                         grid.get_Permittivity(charge.location)))

#                                 Electric_field_charges(charge,
#                                                        location,
#                                                        r,
#                                                        R,
#                                                        grid))
                field_H[time], field_B[time] = (
                    field_H[time] +
                    Magnetic_field_charges(charge,
                                           location,
                                           r,
                                           R,
                                           grid))

    # add everything up
    E_field = E_field + field_E
    H_field = H_field + field_H
    B_field = B_field + field_B

    return(E_field, H_field, B_field)

#Electric_field_currents(grid.current, grid.location, grid.amps, grid.direction, r, R, get_Permittivity)


cdef Electric_field_currents(current_diff_t,
                             list location,
                             double amps,
                             list direction,
                             double r,
                             list R,
                             list Permittivity):
    # find the E field do to currents on the grid
    # see equation 2-2.12 of
    # electromagnetic retardation and theory of relativity
    cdef list scale
    cdef list electrostatic
    
    cdef list electrokinetic
    cdef list E_field

    cdef list permitivity_array = Permittivity
    cdef double perm

    scale = [1 / (4 * pi * E_0 *
                 C_0 * perm) for perm in permitivity_array]
    # first find the terms dependent on current.amp
    # this is the electrostatic field

    electrostatic = [(- sc * amps * dr * 2 / r**2)
                     for (dr, sc) in zip(direction, scale)]

    electrostatic = [(elcs + sc *
                     2 * R_d * amps * dot(
                                    direction,
                                    R) * 2 / (r**4))
                     for (sc, elcs, R_d) in zip(scale, electrostatic, R) ]

    # now find the terms dependent on current.dif_t
    # this is the electrokinetic field

    electrokinetic = [sc *
                      R_d * dot(current_diff_t, R) /
                      (r**3 * C_0) for (R_d, sc) in zip(R, scale)]

    electrokinetic = [(EK - sc *
                      diff_t / (r * C_0))
                      for (EK, sc, diff_t) in zip(electrokinetic, scale, current_diff_t)]

    # find the total E field
    E_field = [ES + EK for (ES, EK) in zip(electrostatic, electrokinetic)]
    
    
    return(E_field)


cdef Magnetic_field_currents(current, r, R, grid):

# cdef Magnetic_field_currents(current, amps, direction, r, R, grid):
    # this function finds the H field do to currents
    # see equation 2-2.5 of
    # electromagnetic retardation and theory of relativity

    # calculate the H field resulting from the currents
    # scale is a factor resulting from geometry
    # and a good place to concider for permibility
    scale = 1 / (4 * pi)

    # first find the terms do to steady currents
    H_field_amps = (scale * current.amps *
                    cross(current.direction, R) * 2 / (r**2))

    # now find the term dependent on current.diff_t
    # this is the electrokinetic field
    H_field_diff_t = (scale * cross(current.diff_t, R) / (r**2 * C_0))

    # find the total H field
    H_field = H_field_amps + H_field_diff_t

    B_field = U_0 * H_field
    return(H_field, B_field)


cdef Electric_field_charges(charge,
                            velocity,
                            acceleration,
                            Q,
                            location,
                            r,
                            R,
                            grid,
                            Permittivity):

    # this function finds the E field resulting from charges
    # see equation 4-4.34 of
    # electromagnetic retardation and theory of relativity

    # scale = 1 / (4 * np.pi * E_0 * grid.get_Permittivity(location))

    # scale = 1 / (4 * np.pi * E_0 * grid.get_Permittivity(charge.location))
    ''' chang this to a sum of equations 4-1.11 and 4-4.32
    this will only requier 2 dot products'''

    # scale = 1 / (4 * pi * E_0 * grid.get_Permittivity(charge.location))
    scale = [1 / (4 * pi * E_0 * perm) for perm in Permittivity]

    #velocity = charge.velocity

    #acceleration = charge.acceleration

    velocity = [v / C_0 for v in velocity]

    velocity_r = [v / r for v in velocity]
    
    dot_velocity_term = dot(R, (velocity_r ))
#    dot_velocity_term = dot(R, (velocity / (r * C_0)))
#
#    velocity = np.array(charge.velocity)
#    
    common_factor = 1 / (1 - dot_velocity_term)

    dot_acceleration_term = dot(R, acceleration)

    E_A_scale = ([(charge.Q * sc *
           common_factor * common_factor
           * (common_factor * dot_acceleration_term * (1/(r*r)))) for sc in scale])

    E_A = [E_A_sc * (R_e / r + v ) + a / r for E_A_sc, R_e, v, a in zip(E_A_scale, R, velocity, acceleration ) ]

    E_B_scale = [charge.Q * sc * common_factor * common_factor * common_factor
           * (1 - dot(velocity, velocity) ) / (r*r*r) for sc in scale ]

    E_B = [E_B_sc * (R_e - r * v) for E_B_sc, R_e, v in zip(E_B_scale, R, velocity)]

    E_field = [E_a + E_b for E_a, E_b in zip(E_A, E_B)]
    # E_field = E_A + E_B
#    return([1,1,1])
    return(E_field)


cdef Magnetic_field_charges(charge, location, r, R, grid):
    # this function finds the H field resulting from charges
    # see equation 4-5.10 of
    # electromagnetic retardation and theory of relativity

    # scale is a factor resulting from geometry
    # and a good place to conscider for permibility
    scale = 1 / (4 * pi)

    velocity = np.array(charge.velocity)

    acceleration = charge.acceleration

    # first find the most used term in the remaning equations
    most_used_factor = 1/(r * (1 - dot(R, (velocity / (r * C_0)))))

    # find the common term
    common_factor = charge.Q * scale * most_used_factor**2

    # next find the velocity factor
    velocity_factor = (1 - np.dot(velocity, velocity) / C_0**2
                       + dot(R, acceleration) / C_0**2)

    field = (cross(acceleration / C_0
                   + most_used_factor * velocity_factor * velocity, R))

    # find the total H field
    H_field = common_factor * (field)

    B_field = U_0 * H_field
    return(H_field, B_field)


cdef Electric_field_dipole(dipole, location, r, R, grid):
    # this is the equation for an electric dipole from page 449
    # of handbook of phisics
    # also see equation 5-4.13 from electricity and magnetism by jefimenko
    scale = 1 / (4 * pi * E_0)

    # first find the term requiering a dot product
    first_term = (3 * dot(dipole.dipole_moment, R) * R) / r**5

    # now find the remaining term
    last_term = dipole.dipole_moment / r**3

    E_field = scale * (first_term - last_term)
    print('E_field test ' + str(E_field))
    print(dipole.dipole_moment)
    return(E_field)


#cdef cross( list A, list B):
cdef cross(A, B):
    return(np.array([A[1] * B[2] - A[2] * B[1],
                     A[2] * B[0] - A[0] * B[2],
                     A[0] * B[1] - A[1] * B[0]]))


cdef dot(A, B):
    return(np.array(A[0] * B[0] + A[1] * B[1] + A[2] * B[2]))
