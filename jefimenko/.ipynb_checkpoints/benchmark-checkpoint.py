'''
  Copyright (C) 2018

  This file is part of Jefimenko.

  Jefimenko is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software Foundation,
  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
'''

grid = Grid(2,  # creat a 2D grid
            delta=[1, 1],  # this is the size of a step in x,y and z
            size=10,  # this is the size of the grid in meters
            # this is the length of the simulation in secounds
            time=14 * 3.36 * 10**-9,
            delta_t=3.36 * 10**-9  # this is the size of a full time step
            )

grid.Add_Charge([0, 0], Q=1)
plot_grid(grid)

simulate(grid)
for t in range(grid.time_size):
    plot_EM_grid('E', time=t)

for T in range(10):
    print(T)
    field_caculator(grid.grid['E'][T], [1, 0])
