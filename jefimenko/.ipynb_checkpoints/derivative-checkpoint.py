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


# this finds a time gradient for a field
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


# this finds the gradient of a field
def space_grad(field, delta, axis=tuple([0, 1, 2])):
    if field.shape[0] == 1:
        field = np.repeat(field[:, :, :], 2,
                          axis=0)

    if field.shape[1] == 1:
        field = np.repeat(field[:, :, :], 2,
                          axis=1)

    if field.shape[2] == 1:
        field = np.repeat(field[:, :, :],
                          2,
                          axis=2)

    gradient = np.gradient(field,
                           delta[0],
                           delta[1],
                           delta[2],
                           axis=axis)
    gradient_R = np.zeros(
            tuple(np.append(np.array(gradient[0]).shape,
                            [3])), dtype='complex')

    for i in np.ndindex(np.array(gradient[0]).shape):
        gradient_R[i] = [gradient[0][i], gradient[1][i], gradient[2][i]]

    # pdb.set_trace()
    return(np.array(gradient_R))


def divergence(field, delta):
    "return the divergence of a n-D field"
    gradient = space_grad(field, delta)
    # return np.sum(np.gradient(field),axis=0)
    div = np.sum(gradient, axis=-1)
    # pdb.set_trace()
    return (div)

    
# call delta as grid.delta
def curl_3d(field, delta):
    curl = []

    grad = np.gradient(field,delta[0],delta[1],delta[2],axis = (0,1,2))
    for i in range(len(grad)):
        for j in range(len(grad[0])):
            for k in range(len(grad[0][0])):
                curl.append([grad[1][i ,j ,k ,2] - grad[2][i ,j ,k ,1],
                            grad[2][i ,j ,k ,0] - grad[0][i ,j ,k ,2] ,
                            grad[0][i ,j ,k ,1] - grad[1][i ,j ,k ,0]])
    return(np.array(curl))
