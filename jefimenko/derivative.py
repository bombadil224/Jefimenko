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
import pdb

def currents_time_diff(grid):
    # this will calculate the time derivative of the currents
    # first make a list of all of the currents
    for n in range(len(grid.currents[0])):
        current_time_list = []
        # notice that amps is included in this so it will not be needed later
        for time in range(grid.time_size):
            current_time_list.append(grid.currents[time][n].direction
                                     * grid.currents[time][n].amps)
        # now calculate the gradient with numpy
        gradient = np.gradient(current_time_list, grid.delta_t, axis=0)

        # now put all of the values into grid.currents for storage
        for time in range(grid.time_size):
            grid.currents[time][n].diff_t = gradient[time]
