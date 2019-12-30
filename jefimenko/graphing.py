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

from scipy.interpolate import interp1d

# this file handels the functions needed for graphing


# ths will plot the grid charges currents and conductors
def plot_grid(Grid, time=0):

    # each number of dimensions is handeld sepratly
    if Grid.dimension == 1:

        x, y, z, Q = [], [], [], []
        # first plot the charges
        for i in range(len(Grid.charges[time])):
            x.append(Grid.charges[time][i].location[0])
            y.append(0)
            z.append(0)
            Q.append(Grid.charges[i][time].Q)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c='r', marker='o')

        x, y, z, Q = [], [], [], []
        for i in range(len(Grid.conductors)):
            x.append(Grid.conductors[i].location[0])
            y.append(0)
            z.append(0)
            # Q.append(Grid.charges[i].Q)

        ax.scatter(x, y, z, c='k', marker='s')

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

        x, y, z = [], [], []
        for i in range(len(Grid.conductors)):
            x.append(Grid.conductors[i].location[0])
            y.append(Grid.conductors[i].location[1])
            z.append(0)

        ax.scatter(x, y, z, c='k', marker='s')

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

        print('this is the grid layout')
        x, y, z = [], [], []
        Q = []

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        try:
            for i in range(len(Grid.charges[time])):
                x.append(Grid.charges[time][i].location[0])
                y.append(Grid.charges[time][i].location[1])
                z.append(Grid.charges[time][i].location[2])
                # Q.append((1 ,0 ,Grid.charges[time][i].Q))
                q = Grid.charges[time][i].Q
                if q > 10:
                    q = 10

                Q.append((1, 0, np.interp(q, [0, 1], [0, 10])))

            ax.scatter(x, y, z, c='r', marker='o')
        except:
            pass

        x, y, z = [], [], []
        for i in range(len(Grid.Permittivity)):
            x.append(Grid.Permittivity[i].location[0])
            y.append(Grid.Permittivity[i].location[1])
            z.append(Grid.Permittivity[i].location[2])

        ax.scatter(x, y, z, c='y', marker='d')

        x, y, z = [], [], []
        for i in range(len(Grid.conductors)):
            x.append(Grid.conductors[i].location[0])
            y.append(Grid.conductors[i].location[1])
            z.append(Grid.conductors[i].location[2])

        ax.scatter(x, y, z, c='k', marker='s')

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
                w)  # ,
                # pivot='middle')

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

        V = []
        # for i in range(len(grid.grid[mode][time])):
        for i in np.ndindex(tuple(grid.shape)):
            V.append(np.linalg.norm(grid.grid[mode][time][i]))
        V_max = max(V)
        V_min = min(V)

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

                m = interp1d([V_min, V_max],  # set the length of the EM vector
                             [np.linalg.norm(grid.delta) * .01,
                              np.linalg.norm(grid.delta) * .5])
                vlength = m(norm)

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

        V = []
        # for i in range(len(grid.grid[mode][time])):
        for i, j in np.ndindex(tuple(grid.shape)):
            V.append(np.linalg.norm(grid.grid[mode][time][i][j]))
        V_max = max(V)
        V_min = min(V)

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

                m = interp1d([V_min, V_max],  # set the length of the EM vector
                             [np.linalg.norm(grid.delta) * .01,
                              np.linalg.norm(grid.delta) * .5])
                vlength = m(norm)

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

        V = []
        for i, j, k in np.ndindex(tuple(grid.shape)):
            V.append(np.linalg.norm(grid.grid[mode][time][i][j][k]))
        V_max = max(V)
        V_min = min(V)

        for i, j, k in np.ndindex(tuple(grid.shape)):
            x = i * grid.delta[0]
            y = j * grid.delta[1]
            z = k * grid.delta[2]

            u = grid.grid[mode][time][int(i)][int(j)][int(k)][0]
            v = grid.grid[mode][time][int(i)][int(j)][int(k)][1]
            w = grid.grid[mode][time][int(i)][int(j)][int(k)][2]

            norm = np.linalg.norm([u, v, w])

            if norm == 0:
                ax.scatter(x, y, z, c='b', marker='o')
            else:
                u, v, w = [u, v, w] / norm

                m = interp1d([V_min, V_max],  # set the length of the EM vector
                             [np.linalg.norm(grid.delta) * .01,
                              np.linalg.norm(grid.delta) * .5])
                vlength = m(norm)

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
