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


from .constants import *

# cdef double pi = 3.141592653589793
# cdef double C_0 = 299792458  # this is the speed of light in meters
# cdef double K_e = 8.9875517873681764 * 10 ** 9  # this is Coulomb's constant
# cdef double E_0 = (4 * pi * K_e) ** -1  # free space permittivity
# cdef double U_0 = (C_0 ** 2 * E_0) ** -1   # free space permeability


def retardation(r, delta_t, c=C_0):

    # cdef int ret_val = 0
    # ret_val = int(np.rint(((r / c) / delta_t)))
    ret_val = int(round(((r / c) / delta_t)))

    return(ret_val)


def field_calculator(charges,
                     currents,
                     dipoles,
                     grid,
                     location,
                     time_0,
                     time_size):
    '''# add time_size to call'''

    # this function maneges the calculaiton of the E and H fields
    # cdef double r
    # cdef double delta_t = grid.delta_t
    # cdef int time = time_0

    # first we need some grids to hold the added field strenth
    field_E = [[0, 0, 0] for i in range(time_size)]
    field_H = [[0, 0, 0] for i in range(time_size)]
    field_B = [[0, 0, 0] for i in range(time_size)]

    # find the field resulating form the currents at time "time"
    for current in (currents):
        # print(vars(current))

        # find the distance from the current to the locaiton of intrest
        R = [loc - C_loc for (loc, C_loc) in zip(location, current.location)]
        r = dot(R, R)**.5

        if r != 0:
            # calculate the time at which to add the new field quantity
            time = time_0 + retardation(r, grid.delta_t)

            # if time is less then the length of the simulaiton
            # calculate the added field
            Permittivity = grid.get_Permittivity
            if time < time_size:
                field_E[time] = ([F + E for (F, E) in zip(field_E[time],
                                 E_field_currents(current.diff_t,
                                                  current.location,
                                                  current.amps,
                                                  current.direction,
                                                  r,
                                                  R,
                                                  Permittivity(location)))])

                H_F, B_F = Magnetic_field_currents(current,
                                                   r,
                                                   R,
                                                   grid)

                field_H[time] = [F_H + H for (F_H, H)
                                 in zip(field_H[time], H_F)]
                field_B[time] = [F_B + B for (F_B, B)
                                 in zip(field_B[time], B_F)]

    # time = time_0
    # now repet the above for charges
    for charge in charges:
        # calcualte the distance to the charge
        R = [loc - C_loc for (loc, C_loc) in zip(location, charge.location)]
        r = dot(R, R)**.5

        # find the field resulating from the charge at time "time"
        if r != 0:
            time = time_0 + retardation(r, grid.delta_t)

            # if time is less then the length of the
            # simulaiton calculate the added field
            if time < time_size:
                field_E[time] = [F + E for
                                 (F, E) in
                                 zip(field_E[time],
                                     E_field_charges(charge,
                                                     charge.velocity,
                                                     charge.acceleration,
                                                     charge.Q,
                                                     location,
                                                     r,
                                                     R,
                                                     grid,
                                                     grid.get_Permittivity(
                                                     charge.location)))]

                H_F, B_F = Magnetic_field_charges(charge,
                                                  location,
                                                  r,
                                                  R,
                                                  grid)

                field_H[time] = [F_H + H for (F_H, H) in
                                 zip(field_H[time], H_F)]
                field_B[time] = [F_B + B for (F_B, B) in
                                 zip(field_B[time], B_F)]

    return(field_E, field_H, field_B)


def E_field_currents(current_diff_t,
                     location,
                     amps,
                     direction,
                     r,
                     R,
                     Permittivity):
    # find the E field do to currents on the grid
    # see equation 2-2.12 of
    # electromagnetic retardation and theory of relativity

    # cdef list scale
    # cdef list electrostatic

    # cdef list electrokinetic
    # cdef list E_field

    # cdef list permitivity_array = Permittivity
    # cdef double perm

    scale = [1 / (4 * pi * E_0 *
                  C_0 * perm) for perm in Permittivity]
    # first find the terms dependent on current.amp
    # this is the electrostatic field

    electrostatic = [(- sc * amps * dr * 2 / r**2)
                     for (dr, sc) in zip(direction, scale)]

    electrostatic = [(elcs + sc *
                     2 * R_d * amps * dot(
                                    direction,
                                    R) * 2 / (r**4))
                     for (sc, elcs, R_d) in zip(scale, electrostatic, R)]

    # now find the terms dependent on current.dif_t
    # this is the electrokinetic field

    electrokinetic = [sc *
                      R_d * dot(current_diff_t, R) /
                      (r**3 * C_0) for (R_d, sc) in zip(R, scale)]

    electrokinetic = [(EK - sc *
                      diff_t / (r * C_0))
                      for (EK, sc, diff_t) in zip(electrokinetic,
                                                  scale,
                                                  current_diff_t)]

    # find the total E field
    E_field = [ES + EK for (ES, EK) in zip(electrostatic,
                                           electrokinetic)]

    return(E_field)


# cdef Magnetic_field_currents(current, amps, direction, r, R, grid):
def Magnetic_field_currents(current, r, R, grid):

    # this function finds the H field do to currents
    # see equation 2-2.5 of
    # electromagnetic retardation and theory of relativity

    # calculate the H field resulting from the currents
    # scale is a factor resulting from geometry
    # and a good place to concider for permibility
    scale = 1 / (4 * pi)

    # first find the terms do to steady currents
    H_field_amps = [scale * current.amps * Cross * 2 / (r**2) for Cross in
                    cross(current.direction, R)]

    # now find the term dependent on current.diff_t
    # this is the electrokinetic field
    H_field_diff_t = [scale * Cross / (r**2 * C_0)
                      for Cross in cross(current.diff_t, R)]

    # find the total H field
    H_field = [H_f + H_d for (H_f, H_d) in zip(H_field_amps, H_field_diff_t)]

    B_field = [U_0 * H for H in H_field]
    return(H_field, B_field)


def E_field_charges(charge,
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

    ''' chang this to a sum of equations 4-1.11 and 4-4.32
    this will only requier 2 dot products'''

    scale = [1 / (4 * pi * E_0 * perm) for perm in Permittivity]

    velocity = [v / C_0 for v in velocity]

    velocity_r = [v / r for v in velocity]

    dot_velocity_term = dot(R, (velocity_r))

    common_factor = 1 / (1 - dot_velocity_term)

    dot_acceleration_term = dot(R, acceleration)

    E_A_scale = ([(charge.Q * sc *
                   common_factor * common_factor
                   * (common_factor * dot_acceleration_term * (1/(r*r))))
                  for sc in scale])

    E_A = [E_A_sc * (R_e / r + v) + a / r
           for E_A_sc, R_e, v, a in zip(E_A_scale,
                                        R,
                                        velocity,
                                        acceleration)]

    E_B_scale = [charge.Q * sc * common_factor * common_factor * common_factor
                 * (1 - dot(velocity, velocity)) / (r*r*r) for sc in scale]

    E_B = [E_B_sc * (R_e - r * v) for
           E_B_sc, R_e, v in zip(E_B_scale, R, velocity)]

    E_field = [E_a + E_b for E_a, E_b in zip(E_A, E_B)]
    return(E_field)


def Magnetic_field_charges(charge, location, r, R, grid):
    # this function finds the H field resulting from charges
    # see equation 4-5.10 of
    # electromagnetic retardation and theory of relativity

    # scale is a factor resulting from geometry
    # and a good place to conscider for permibility
    scale = 1 / (4 * pi)

    velocity = charge.velocity

    acceleration = charge.acceleration

    # first find the most used term in the remaning equations
    most_used_factor = 1/(r * (1 - dot(R, velocity) / (r * C_0)))

    # find the common term
    common_factor = charge.Q * scale * most_used_factor**2

    # next find the velocity factor
    velocity_factor = (1 - dot(velocity, velocity) / C_0**2
                       + dot(R, acceleration) / C_0**2)

    field = cross([A / C_0 + V * most_used_factor * velocity_factor
                   for (A, V) in zip(acceleration, velocity)], R)

    # find the total H field
    H_field = [common_factor * F for F in field]

    B_field = [U_0 * H for H in H_field]
    return(H_field, B_field)


def Electric_field_dipole(dipole, location, r, R, grid):
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


# cdef cross( list A, list B):
def cross(A, B):
    return(([A[1] * B[2] - A[2] * B[1],
             A[2] * B[0] - A[0] * B[2],
             A[0] * B[1] - A[1] * B[0]]))


def dot(A, B):
    return((A[0] * B[0] + A[1] * B[1] + A[2] * B[2]))
