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
from .charge_current_integrals import *
from .constants import *
from .classes import *
import numpy as np
import sys
from math import isnan as isnan
import time

import pdb


def timed_range(number_set):  # this will keep trak of how far compleat the simulaiton is
    n = 0
    m = number_set
    for i in range(number_set):
        b = "percent compleat = " + str(100 * (i + 1)/ number_set)
        n = n + 1
        if i != number_set - 1:
            print(b, end="\r")
            yield(i)
        else:
            print(b, end="\r")
            print('compleat, '+ str(b) + '%                    ')
            yield(i)


def simulate(grid, charge_currents=False, induction=False):
    print('simulating grid')
    # this tells if the differentials need to be updated
    update_diff = True
    # this is used to flush the print statments befor the loop starts
    sys.stdout.flush()
    # this will run the simulation of a grid and charges

    # start a loop that will simulate the system for every point on the grid
    # location is the location on the grid
    # where the E field is being calculated
    for time in timed_range(grid.time_size):
        # update_diff must be tested at the start of each time loop incase
        # something is changed. this will be needed when moving charges
        # are introduced
        if update_diff is True:
            currents_time_diff(grid)
            update_diff = False

        for location in np.ndindex(tuple(grid.shape)):

            # find the part of E that depends on charges
            grid.grid['E'] = dynamic_charges_E(grid.grid['E'],
                                               grid.charges[time],
                                               grid,
                                               location,
                                               time)

            # next find the part of E that is dependent on currents
            grid.grid['E'] = electric_currents(grid.grid['E'],
                                               grid.currents[time],
                                               grid,
                                               location,
                                               time)

            # now find H
            # first find the part of H that is dependent on currents
            grid.grid['H'] = currents(grid.grid['H'],
                                      grid.currents[time],
                                      grid,
                                      location,
                                      time)

            # now find the part of H that depends on moving charges
            grid.grid['H'] = dynamic_charges_H(grid.grid['H'],
                                               grid.charges[time],
                                               grid,
                                               location,
                                               time)

        if induction is False:
            pass
        elif induction is True:  # find induction currents in conductors
            induction_currents(grid,
                               grid.currents[time],
                               grid.charges[time],
                               grid.conductors,
                               time)
#######################
#  this modifies charges for each current see part removed from classes cunductors for further details
#        if charge_currents is False:
#            pass
#        elif charge_currents is True:   # this will add charges resulting
#                                        # from the currents in cundoctors
#            current_charges(grid,
#                            grid.currents[time],
#                            grid.conductors, time)
#
#  this part removed untill further notice
#######################
    print('grid simulated')


#def current_charges(grid, currents, conductors, time):
#    # this will make cunductors modify charges and currents
#    for conductor in conductors:
#        for i in range(3):
#
#            # find charges to modify and modifying current
#            charge_1 = grid.charges[time][conductor.charges[2*i]]
#            charge_2 = grid.charges[time][conductor.charges[2*i+1]]
#            current = grid.currents[time][conductor.current]
#
#            charge_modifyer = (np.inner(current.direction, e[i]) *
#                               current.amps * grid.delta_t)
#
#            #  calculate new values for charges on cunductor
#            grid.modify_charge(conductor.charges[2*i],  # charge to modify
#                               time=time + 1,  # time to modify the charge at
#                               velocity=False,  # the new velocity
#                               acceleration=False,  # the new acceleration
#                               location=False,  # the charges locaiton
#                               Q=(charge_1.Q + charge_modifyer),  # coulombs
#                               print_all=False)  # set true to print info
#
#            grid.modify_charge(conductor.charges[2*i+1],  # charge modify
#                               time=time + 1,  # time to make changes
#                               velocity=False,  # new velocity
#                               acceleration=False,  # acceleration
#                               location=False,  # charges locaiton
#                               Q=charge_2.Q - charge_modifyer,  # coulombs
#                               print_all=False)  # set true to print info


def induction_currents(grid,  # calculate induction
                       currents,
                       charges,
                       conductors,
                       time):

    for conductor in conductors:

        if len(charges) != 0:
            conductor.EK_field = EK_charge_field(conductor.EK_field,
                                                 charges,
                                                 grid,
                                                 conductor.location,
                                                 time)
        if len(currents) != 0: 
            conductor.EK_field = EK_current_field(conductor.EK_field,
                                                  currents,
                                                  grid,
                                                  conductor.location,
                                                  time)

        if np.linalg.norm(conductor.EK_field) > 0:

            direction = (conductor.EK_field[time] /
                        np.linalg.norm(conductor.EK_field))
            # amps = (np.linalg.norm(E_field) /
            #        (conductor.Ohm * np.inner(grid.delta, direction)))
            amps = (np.linalg.norm(conductor.EK_field))

            grid.Modify_Current(conductor.current,  # current to modify
                                time=time * grid.delta_t,  # time to modify
                                direction=direction,  # the new direction
                                amps=amps,  # the new amps value
                                print_all=False)

            print('the new amps value at time ' + str(time) + ' is :  ' + str(amps))
            print('the new direction at time ' + str(time) + ' is :  ' + str(direction))
    print('induction_currents functinal')
    print('')

# the electrokinetic field genorated by charges
def EK_charge_field(field,
                    charges,
                    grid,
                    location,
                    time_0):
    
    # note that charges are considerd to be at a porticular location
    scale = 1 / (4 * np.pi * E_0 * C_0**2)

    # this will calculate the field genorated by the charges
    # loop over the charges in the grid
    # this uses equation 4-4.34
    # from "Electromagnetic retardation and theory of relativity"
    # unfortuntly the electrostatic and electrokinetic fields are
    # heavely linked in this situation making seporationg them diffucalt
    for charge in charges:
        R = grid.delta * location
        R = location - charge.location
        # R = grid.delta * R
        r = np.linalg.norm(R)
        R = np.pad(R, (0, 3 - len(R)), 'constant', constant_values=0)

        velocity = np.pad(charge.velocity, (0, 3 - len(charge.velocity)),
                          'constant', constant_values=0)

        # we will need a unit vector in the direction of motion
        velocity_unit = velocity / np.linalg.norm(velocity)

        acceleration = (np.pad(charge.acceleration,
                               (0, 3 - len(charge.acceleration)),
                               'constant', constant_values=0))

        if r != 0:

            # first find the factor common to all terms
            common_factor = 1 / (r * (1 - np.dot(R, (velocity /
                                                    (r * C_0)))))

            # find the common terms of the remaining parts
            radius_factor = charge.Q * commen_factor * scale

            # next find the factor resalting from velocity only
            velocity_term = (radius_factor * velocity**2 * (
                             np.dot(R, velocity.norm) -
                             r * velocity / C_0)) * common_factor**2

            # now find the term that is dependent on acceleration
            acceleration_term = radius_factor * (-1) * (acceleration +
                                  (np.dot(R, acceleration) * R) /
                                  (c_0 * r))

            # the only case of error should be a time that is to big
            # use a try to find and in this case exit the loop
            try:
                time = time_0 + retardation(r, grid)
                field[time] = (field[time] +
                               velocity_term +
                               acceleration_term)

            except IndexError:
                pass
        else:
            field[time_0] = field[time_0]
    return (field)
    pass


# the electrokinetic field genorated by currents
def EK_current_field(EK_field, currents, grid, location, time_0):
    scale = np.linalg.norm(grid.delta) / (4 * np.pi * E_0 * C_0)

    # this will calculate the field genorated by the constant currents
    # loop over the currents in the grid
    for current in currents:
        R = grid.delta * location
        R = location - current.location
        # R = grid.delta * R
        r = np.linalg.norm(R)
        R = np.pad(R, (0, 3 - len(R)), 'constant', constant_values=0)

        if r != 0:
            try:
                # first find the terms dependent on current.amp
                # this is the electrostatic field
                time = time_0 + retardation(r, grid)
                # now find the terms dependent on current.dif_t
                # this is the electrokinetic field

                electrokinetic = (scale *
                                  R * np.dot(current.diff_t, R) /
                                  (r**3 * C_0))

                electrokinetic = (electrokinetic - scale *
                                           current.diff_t / (r * C_0))

                # if we are just finding the field of a cunductor
                # we are only intorested in the electrokinetic field
                # for calculationg currents
                EK_field[time] = EK_field[time] + electrokinetic

            except IndexError:
                pass
        else:
            EK_field[time_0] = EK_field[time_0]

    # this is used to find the electrokinetic field for a conductor

    return(EK_field)
