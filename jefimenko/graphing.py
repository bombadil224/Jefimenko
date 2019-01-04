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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


# ths will plot the grid charges and currents
def plot_grid(Grid, time=0):

    # each number of dimen stons is handeld sepratly
    if Grid.dimension == 1:

        x, y, z, Q = [], [], [], []
        # first plot the charges
        for i in range(len(Grid.charges[time])):
            x.append(Grid.charges[time][i].location[0])
            y.append(0)
            z.append(0)
            # Q.append(Grid.charges[i].Q)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c='r', marker='o')

        for i, txt in enumerate(Q):
            ax.annotate(txt, (x[i], y[i]))

        x, y, z, u, v, w = [], [], [], [], [], []
        # secound plot the currents
        for i in range(len(Grid.currents[time])):
            x.append(Grid.currents[time][i].location[0])
            y.append(0)
            z.append(0)

            u.append(Grid.currents[time][i].direction[0])
            v.append(Grid.currents[time][i].direction[1])
            w.append(Grid.currents[time][i].direction[2])

        vlength = .4 * np.linalg.norm(Grid.delta)

        ax.quiver(
                x,
                y,
                z,
                u,
                v,
                w,
                pivot='middle',
                length=vlength,
                arrow_length_ratio=0.3 / vlength)

        # Do any fineshing touches on the grid
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        ax.set_xlim(-1, Grid.size[0])
        ax.set_ylim(-1, 1)
        ax.set_zlim(-1, 1)

        print('shape = 1')
    elif Grid.dimension == 2:

        x, y, z, Q = [], [], [], []
        # first plot the charges
        for i in range(len(Grid.charges[time])):
            x.append(Grid.charges[time][i].location[0])
            y.append(Grid.charges[time][i].location[1])
            z.append(0)
            # Q.append(Grid.charge[i].Q)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c='r', marker='o')

        x, y, z, u, v, w = [], [], [], [], [], []
        # secound plot the currents
        for i in range(len(Grid.currents[time])):
            x.append(Grid.currents[time][i].location[0])
            y.append(Grid.currents[time][i].location[1])
            z.append(0)

            u.append(Grid.currents[time][i].direction[0])
            v.append(Grid.currents[time][i].direction[1])
            w.append(Grid.currents[time][i].direction[2])

        vlength = .4 * np.linalg.norm(Grid.delta)

        ax.quiver(
                x,
                y,
                z,
                u,
                v,
                w,
                pivot='middle',
                length=vlength,
                arrow_length_ratio=0.3 / vlength)

        # Do any fineshing touches on the grid
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        ax.set_xlim(-1, Grid.size[0])
        ax.set_ylim(-1, Grid.size[1])
        ax.set_zlim(-1, 1)

        print('Shape = 2')

    elif Grid.dimension == 3:

        x, y, z = [], [], []
        Q = []
        for i in range(len(Grid.charges[time])):
            x.append(Grid.charges[time][i].location[0])
            y.append(Grid.charges[time][i].location[1])
            z.append(Grid.charges[time][i].location[2])
            # Q.append(Grid.charge[i].Q)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c='r', marker='o')

        x, y, z, u, v, w = [], [], [], [], [], []
        # secound plot the currents
        for i in range(len(Grid.currents[time])):
            x.append(Grid.currents[time][i].location[0])
            y.append(Grid.currents[time][i].location[1])
            z.append(Grid.currents[time][i].location[2])

            u.append(Grid.currents[time][i].direction[0])
            v.append(Grid.currents[time][i].direction[1])
            w.append(Grid.currents[time][i].direction[2])

        vlength = .4 * np.linalg.norm(Grid.delta)

        ax.quiver(
                x,
                y,
                z,
                u,
                v,
                w,
                pivot='middle',
                length=vlength,
                arrow_length_ratio=0.3 / vlength)

        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        ax.set_xlim(-1, Grid.size[0])
        ax.set_ylim(-1, Grid.size[1])
        ax.set_zlim(-1, Grid.size[2])

        print('shape = 3')

    plt.show()
    sys.stdout.flush()


def plot_EM_grid(mode, grid, time=0):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if grid.dimension == 1:

        for i in np.ndindex(tuple(grid.shape)):
            x = i * grid.delta[0]
            y = 0
            z = 0

            u = grid.grid[mode][time][i][0].real
            v = grid.grid[mode][time][i][1].real
            w = grid.grid[mode][time][i][2].real

            norm = np.linalg.norm([u, v, w])
###############
            V = []
            for i in range(len(grid.grid[mode][time])):
                V.append(np.linalg.norm(grid.grid[mode][time][i]))
            V_max = max(V)
###############
            if norm == 0:
                ax.scatter(x, y, z, c='b', marker='o')
            else:
                u, v, w = [u, v, w] / norm
                vlength = np.linalg.norm(grid.delta) * .4 * norm / V_max

                ax.quiver(
                    x,
                    y,
                    z,
                    u,
                    v,
                    w,
                    pivot='middle',
                    length=vlength)

        for i in range(len(grid.charges[time])):
            x = (grid.charges[time][i].location[0])
            y = 0
            z = 0
            ax.scatter(x, y, z, c='r', marker='o')
        ax.set_ylim([-.5, .5])
        ax.set_zlim([-.5, .5])

    elif grid.dimension == 2:

        for i, j in np.ndindex(tuple(grid.shape)):
            x = i * grid.delta[0]
            y = j * grid.delta[1]
            z = 0

            u = grid.grid[mode][time][int(i)][int(j)][0].real
            v = grid.grid[mode][time][int(i)][int(j)][1].real
            w = grid.grid[mode][time][int(i)][int(j)][2].real

            norm = np.linalg.norm([u, v, w])

            if norm == 0:
                ax.scatter(x, y, z, c='b', marker='o')
            else:
                u, v, w = [u, v, w] / norm
                vlength = .4 * np.linalg.norm(grid.delta)

                ax.quiver(
                    x,
                    y,
                    z,
                    u,
                    v,
                    w,
                    pivot='middle',
                    length=vlength)

        for i in range(len(grid.charges[time])):
            x = (grid.charges[time][i].location[0])
            y = (grid.charges[time][i].location[1])
            z = 0
            ax.scatter(x, y, z, c='r', marker='o')

        ax.set_zlim(-1.5, 1.5)

    elif grid.dimension == 3:

        for i, j, k in np.ndindex(tuple(grid.shape)):
            x = i * grid.delta[0]
            y = j * grid.delta[1]
            z = k * grid.delta[2]

            u = grid.grid[mode][time][int(i)][int(j)][int(k)][0].real
            v = grid.grid[mode][time][int(i)][int(j)][int(k)][1].real
            w = grid.grid[mode][time][int(i)][int(j)][int(k)][2].real

            norm = np.linalg.norm([u, v, w])

            if norm == 0:
                ax.scatter(x, y, z, c='b', marker='o')
            else:
                u, v, w = [u, v, w] / norm
                vlength = np.linalg.norm(grid.delta) * .4

                ax.quiver(
                    x,
                    y,
                    z,
                    u,
                    v,
                    w,
                    pivot='middle',
                    length=vlength)

            # u, v, w = [u, v, w] / np.linalg.norm([u, v, w])

            # vlength = np.linalg.norm(np.array([u, v, w]))
            vlength = np.linalg.norm(grid.delta) * .4

            ax.quiver(
                x,
                y,
                z,
                u,
                v,
                w,
                pivot='middle',
                length=vlength,
                arrow_length_ratio=0.3 / vlength)

        for i in range(len(grid.charges[time])):
            x = (grid.charges[time][i].location[0])
            y = (grid.charges[time][i].location[1])
            z = (grid.charges[time][i].location[2])

        ax.set_xlim(0, grid.size[0])
        ax.set_ylim(0, grid.size[1])
        ax.set_zlim(0, grid.size[2])

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

    print('The acutal time is ' + str(time * grid.delta_t))


def field_calculator(input_field, Grid):
    # find the point of minimal field norm and calculate everything form there
    field = {}
    field['vector'] = []
    field['distance'] = []
    field['field'] = []
    field['norm'] = []
    for i in np.ndindex(tuple(Grid.shape)):
        r_vec = i * Grid.delta
        r_norm = np.linalg.norm(r_vec)
        field_r = input_field[tuple(i)]
        field_norm = np.linalg.norm(field_r)

        field['vector'].append(r_vec)
        field['distance'].append(r_norm)
        field['field'].append(field_r)
        field['norm'].append(field_norm)
    return(field)


def print_location(field, location, Grid):
    # locaiton is given in x, y, z distances not in i, j, k array indexes
    # the grid that is pased to this function must not include a time axis
    location = np.array(location)
    index = (location / Grid.delta).astype(int)
    vector = field[tuple(index)] / np.linalg.norm(field[tuple(index)])

    print(' the field is = ' + str(field[tuple(index)]))
    print(' the vector is = ' + str(vector))
    print(' the nomalized field is = ' + str(np.linalg.norm(
                                             field[tuple(index)])))

    return(field[tuple(index)], vector, np.linalg.norm(field[tuple(index)]))
