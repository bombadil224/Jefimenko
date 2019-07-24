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
from .boundary import *
from .polarization import *

import numpy as np
import sys
from math import isnan as isnan
import time as Time


# this file handels the acual simulation


# this will keep trak of how far compleat the simulaiton is
def timed_range(number_set):
    n = 0
    m = number_set
    for i in range(number_set):
        b = "percent compleat = " + str(100 * (i + 1) / number_set)
        n = n + 1
        if i != number_set - 1:
            print(b, end="\r")
            yield(i)
        else:
            print(b, end="\r")
            print('compleat, ' + str(b) + '%                    ')
            yield(i)


# this will run the acual simulation of a grid
def simulate(grid,
             charge_currents=False,
             induction=False,
             magnetization=False,
             Test=True):

    print('simulation started')
    start = Time.time()

    # this tells if the differentials need to be updated
    update_diff = True

    # befor we run the simulation see if we are using Permeability
    if grid.free_space is False:
        # first find the gradient of the polarization
        grid.Permittivity_normals = normals(grid.Permittivity,
                                            grid.shape,
                                            grid.delta)
        pass

    # find every region of the grid
    # regions, boundaries = genorate_regions(grid)

    ########################################################
    # this code uses the idea of using regions for calculationg fields
    '''' for a first use a mothed of polarization will be used
    the porimary use of this will be to insure correct porapugation
    speed when C changes '''
    # this is test code and needs replaced
    regions = []
    region = []

    boundaries = []
    boundary = []

    for i in grid.boundary:
        boundary.append(i)

    for i in (np.ndindex(tuple(grid.shape))):
        region.append(i)
    regions.append(region)
    boundaries.append(boundary)

    ########################################################

    # this is used to flush the print statments befor the loop starts
    # this works in combination with grid.time_size and may not be needed
    sys.stdout.flush()

    # start a loop that will simulate the system for every point on the grid
    # location is the location on the grid
    # where the E and H field is being calculated
    print('stage one charges and currents')
    for time in timed_range(grid.time_size):

        # now calculate the polarization field and all genorated dipoles
        # from it at time, time
        # note that this should be the first thing that is done

        # the polarization field is realy only kept for future analysis

        # update_diff must be tested at the start of each time loop incase
        # something is changed. this will be needed when moving charges
        # are introduced
        if update_diff is True:
            currents_time_diff(grid)
            update_diff = False

        # this will have to be changed to a region by region calculation
        # porabubly best to move it to a new modual

        # for location in np.ndindex(tuple(grid.shape)):
        '''' these next two for loops are used for mothed of regions
        this is needed to insure porapugation speed'''

        for (region, boundary) in zip(regions, boundaries):
            for location in region:

                # make sure that the location that is used
                # is the acual location not just its index
                location_dx = location * grid.delta

                # calculate the fields strength at this location and time
                (grid.grid['E'][:, location[0], location[1], location[2]],
                 grid.grid['H'][:, location[0], location[1], location[2]],
                 grid.grid['B'][:, location[0], location[1], location[2]],) = (
                 field_calculator(grid.grid['E'][:, location[0],
                                                 location[1], location[2]],
                                  grid.grid['H'][:, location[0],
                                                 location[1], location[2]],
                                  grid.grid['B'][:, location[0],
                                                 location[1], location[2]],
                                  grid.charges[time],
                                  grid.currents[time],
                                  grid.dipoles[time],
                                  grid, location_dx,
                                  time))

########################################################
        # this code uses the idea of using regions for calculationg fields
        '''' this code is needed to isure porapugation speed in different media
        note that this needs to be done as the simulation procedes through time
        to insure that things porapugate poraporly as this will efect the
        polarization field '''

        ''' now simulate the boundary of every region this must be done
        inside of the first time loop so that the polarization field is
        poraporly updated '''
        for (region, boundary) in zip(regions, boundaries):
            grid.grid['E'], grid.grid['H'] = boundary_simulation(region,
                                                                 boundary,
                                                                 time,
                                                                 grid.grid['E'],
                                                                 grid.grid['H'],
                                                                 grid)
###############################################################################

    if grid.free_space is False:  # this calculates the effects of permittivity
        print('stage two permittivity')

        # this calculates the effects of permittivity
        ''' we now integrate over the hole thing again
            this time finding the polorization field '''

        for time in timed_range(grid.time_size):

            # first calculate the polarization vector
            grid.grid['P'][time] = polarization_field(grid.grid['E'][time],
                                                      grid.grid['P'][time],
                                                      grid.Permittivity,
                                                      grid.shape, grid)

            grid.grid_P_div = divergence(grid.grid['P'][time], grid.delta)

            for charge in grid.charges[time]:
                charge_loc_index = tuple((charge.location /
                                          grid.delta).astype(int))

                grid.grid_P_div[charge_loc_index] = [0, 0, 0]
                grid.grid['P'][time][charge_loc_index] = [0, 0, 0]

                pass

            for location in region:  # now find the acual polorization effect

                # grid.grid['P_E'] = grid.grid['P_E'] + grid.grid['E']

                location_dx = location * grid.delta

                grid.grid['P_E'][:,
                                 location[0],
                                 location[1],
                                 location[2]] = (
                    grid.grid['P_E'][:,
                                     location[0],
                                     location[1],
                                     location[2]] +
                    polarization_effect(grid.grid['E'][:,
                                        location[0],
                                        location[1],
                                        location[2]],
                                        grid.grid['P'][time],
                                        grid.grid_P_div,
                                        grid.Permittivity_normals,
                                        location_dx,
                                        time,
                                        grid.delta, grid))
                pass

        grid.grid['E'] = grid.grid['E'] + grid.grid['P_E']

    end = Time.time()
    print('grid simulated in ' + str(end - start) + ' seconds')
    print()


def induction_currents(grid,  # calculate induction
                       currents,
                       charges,
                       conductors,
                       time):

    for conductor in conductors:

        # first find the electorkinetic field
        # resulting from the charges and currents
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

        # if this is not zero calculate the induction currents
        # at the next time step

        if np.linalg.norm(conductor.EK_field) > 0:

            direction = (conductor.EK_field[time] /
                         np.linalg.norm(conductor.EK_field[time]))
            # amps = (np.linalg.norm(E_field) /
            #        (conductor.Ohm * np.inner(grid.delta, direction)))
            if time > 0:
                amps = ((np.linalg.norm(conductor.EK_field[time] -
                                        conductor.EK_field[time - 1]))
                        / conductor.Ohm)
            else:
                amps = 0

            if grid.time_size > time + 1:
                grid.Modify_Current(conductor.current,  # current to modify
                                    time=((time + 1) * grid.delta_t),  # time
                                    direction=direction,  # the new direction
                                    amps=amps,  # the new amps value
                                    print_all=False)

            print('the new amps value at time ' + str(time + 1) + ' is :  ' +
                  str(amps))
            print('the new direction at time ' + str(time + 1) + ' is :  ' +
                  str(direction))

    print('induction_currents functinal')
    print('')
