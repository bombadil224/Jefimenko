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

from numpy import pi as pi
from numpy import array as array

# this file holds the constants that are needed

C_0 = 299792458  # this is the speed of light in meters per secound
K_e = 8.9875517873681764 * 10 ** 9  # this is Coulomb's constant
E_0 = (4 * pi * K_e) ** -1  # free space permittivity
U_0 = (C_0 ** 2 * E_0) ** -1   # free space permeability

e = array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # these are the unit vector in R3
