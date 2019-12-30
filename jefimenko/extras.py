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

import time as Time
import dill as pickle
import os


def timed_range(start,           # the point to start genorationg at
                end,             # the end of needed time steps
                start_time,      # acual time of start of simulaiton
                print_time=True):
    n = 0
    for i in range(start, end):
        b = ("percent compleat = " + str(round(100 * (i + 1) / (end - start)))
             + ' estemated time remaining '
             + str(round((end - start) * (Time.time() - start_time)
                         / (i + 1) - (Time.time() - start_time)))
             + ' seconds              ')
        n = n + 1
        if i != end:
            if print_time is True:
                print(b, end="\r")
            yield(i)
        else:
            if print_time is True:
                print(b, end="\r")
                print('compleat, ' + str(b) + '% estemated time remaining '
                      + str((start_time - Time.time()) * 100 / b + '      '))
            yield(i)


def save_work(mapping, file):
    file_out = os.getcwd() + '/' + file
    pickle_out = open(file_out, "wb")

    pickle.dump(mapping, pickle_out)
    pickle_out.close()
    pass


def load_work(file):
    file_in = os.getcwd() + '/' + file
    pickle_in = open(file_in, "rb")
    mapping = pickle.load(pickle_in)
    return mapping


# test for approximate equality (for floating point types)
def arreqclose_in_list(myarr, list_arrays):
    return next((True for elem in list_arrays
                 if elem.size == myarr.size and
                 allclose(elem, myarr)), False)
