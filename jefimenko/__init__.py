#!/usr/bin/env python
""" Jefimenko - An EM simulator based on the Jefimenko equations

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

# this file initializes the program

__version__ = '0.0.5'
__author__ = 'Nehemiah Null'
__contact__ = '<Bombadil224@gmail.com>'
__license__ = "GPLv3"
__date__ = '2018/12/02'


# from .extras import load_work, save_work
from .simulation import *
# from .plasma import *
from .classes import *
from .graphing import *
from .plasma import *

print('Jefimenko version 0.0.5')
