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

# this file contains the classes need to define the grid being simulated
import sys
import types
from .constants import *

# import math

import pdb

class Grid():

    # the grid holds everything being simulated
    def __init__(
            self,
            delta=False,  # this is the size of a step in x, y and z
            size=False,  # this is the size of the grid in meters
            time=.1,  # this is the length of a simulation in secounds
            delta_t=.1,  # this is the size of a full time step
            free_space=True,
            constant_E=[0, 0, 0],
            constant_H=[0, 0, 0],
            plasma_simulaiton=False):

        self.max_charge_distance = 0

        self.constant_E = constant_E   # this defines a constant E field

        self.constant_H = constant_H   # this defines a constant H field

        self.free_space = free_space
        self.dimension = 3

        # define a delta value for the grid in each direction
        if delta is False:
            self.delta = [1.0 for i in range(self.dimension)]
        else:
            delta = np.pad(delta, (0, 3 - len(delta)),
                           'constant',
                           constant_values=1)
            self.delta = delta

        self.delta = np.array(self.delta).astype(float)

        # define the size of the grid in each direction
        if size is False:
            self.size = [10] * self.dimension
            print('grid size diffalting to ' + str(self.size))

        elif len(size) != self.dimension:
            for i in range(len(size), self.dimension):
                size.append(self.delta[i])

            self.size = size
            print(' grid size is ' + str(size))
            print('recommended to give all three measurements')
        else:
            self.size = size

        self.size = np.array(self.size)

        # this is the amount of time the simulation will simulate
        self.time = time
        self.delta_t = delta_t  # this is the size of a time step in seconds

        # this is the size of the grids time axis
        self.time_size = np.array(self.time / self.delta_t).astype(int)

        # next the shape of the grids spacial axis
        self.shape = np.array((self.size / self.delta).astype(int))

        # now the total grid
        self.shape_t = np.append(self.time_size, self.shape)

        # this will be the charges
        self.charges = [[] for i in range(self.time_size)]
        # this will be the currents
        self.currents = [[] for i in range(self.time_size)]
        # this will be the dipoles
        self.dipoles = [[] for i in range(self.time_size)]

        # this will hold the boundary conditions
        self.boundary = [[] for i in range(self.time_size)]

        # this is the conductors that are on the grid
        # notice that these don't change over the course of the simulation
        self.conductors = []

        # now define the relative permittivity grid this is also known as
        # as the dielectric constant

        ''' not used in stead use self.Permittivity_list '''
        # self.R_Permittivity = np.ones(self.shape)
        # self.Permittivity_gradient = np.zeros(self.shape)

        ''' this will be the list of all Permittivity elements
        locations and thier dipoles '''
        self.Permittivity = []

        ''' next the list of all Permittivity unit normals
        this is a class seporate from Permittivity '''
        self.Permittivity_normals = []

        self.Permeability = []

        self.Permeability_normals = []
        # now define the relative permeability grid this is also known as
        # as the dielectric constant
        # notice that Gyromagnetism could be added by making the
        # permeability into a matrix and modifying polarization.py as needed!
        # this will not be done at this time

        # this will be the location of all nonzero Permeability gradients

        # the grid layout is grid.grid[field][time axis][x, y, z axis]
        # the grid will hold the E and H fields as 3D vectors for all casses

        self.grid = {}

        if plasma_simulaiton is False:
            self.grid['E'] = np.zeros(
                tuple(np.append(self.shape_t, [3])), dtype='complex')
            self.grid['H'] = np.zeros(
                tuple(np.append(self.shape_t, [3])), dtype='complex')
            self.grid['B'] = np.zeros(
                tuple(np.append(self.shape_t, [3])), dtype='complex')

            if constant_E is not False:
                for i in np.ndindex(self.grid['E'].shape):
                    if type(self.constant_E) == types.FunctionType:
                        self.grid['E'][i[-3:]] = (constant_E(i[-3:] *
                                                             self.delta,
                                                             i[-4] *
                                                             self.delta_t))
                    else:
                        self.grid['E'][i] = constant_E[i[-1]]
    
            if constant_H is not False:
                for i in np.ndindex(self.grid['H'].shape):
                    if type(self.constant_H) == types.FunctionType:
                        self.grid['H'][i[-3:]] = (constant_H(i[-3:] *
                                                             self.delta, i[-4] *
                                                             self.delta_t))
                    else:
                        self.grid['H'][i] = constant_H[i[-1]]
    
    
            if self.free_space is False:
                # polarization vector grid
                self.grid['P'] = np.zeros(
                    tuple(np.append(self.shape_t, [3])), dtype='complex')
    
                # field genorated by the polarization vector grid
                self.grid['P_E'] = np.zeros(
                    tuple(np.append(self.shape_t, [3])), dtype='complex')
    
                # permitivity gradient grid
                # self.grid_P_grad = np.zeros(
                #    tuple(np.append(self.shape, [3])), dtype='complex')
    
                # magnetization vector grid
                self.grid['M'] = np.zeros(
                    tuple(np.append(self.shape_t, [3])), dtype='complex')
    
                # field genorated by the polarization vector grid
                self.grid['M_H'] = np.zeros(
                    tuple(np.append(self.shape_t, [3])), dtype='complex')
    
                # permebility gradient grid
                self.grid_M_grad = np.zeros(
                    tuple(np.append(self.shape, [3])), dtype='complex')

    # this is from an old attempt to add Permittivity but may be usefull
    def Modify_Permittivity(self,
                            Permittivity,
                            location,  # location is given in acual locaiton
                            print_all=False):

        loc = np.array(location / self.delta).astype(int)
        self.R_Permittivity[tuple(location)] = float(Permittivity)

        if(print_all is True):
            print('Permittivity')
            print(self.R_Permittivity)

    # this is from an old attempt to add Permeability but may be usefull
    def Modify_Permeability(self,
                            Permeability,
                            location,
                            print_all=False):

        loc = np.array(location / self.delta).astype(int)
        self.R_Permeability[tuple(location)] = float(Permeability)

        if(print_all is True):
            print('Permeability')
            print(self.R_Permeability)

    def Add_Charge(self,
                   location,   # the locaiton to place the new charge
                   Q=(-1.6021766208) * 10**(-19),  # charge defalt electron
                   velocity=False,  # the velocity of the charge
                   acceleration=False,  # the acceleration of the charge
                   time=0,     # the time at which to place the charge
                   mass=9.10938356 * 10**(-31),  # kg defalt electron
                   print_all=False,  # print the specs of the charge
                   count=False):  # return the index of the new charge
        # this will add a charge to the grid notice that a charge is a class

        # time when the charg is bing added
        time_N = int(np.rint(time / self.delta_t))

        Q = float(Q)

        location = np.pad(location, (0, 3 - len(location)),
                          'constant', constant_values=0)
        location = np.array(location).astype(float)

        # this defines the charg as a moving with some velocity
        if velocity is False:
            velocity = self.dimension * [0]
        elif len(velocity) != self.dimension:
            velocity = np.pad(velocity, (0, 3 - len(velocity)),
                              'constant', constant_values=0)

        # this defined the charge to have some acceleration
        if acceleration is False:
            acceleration = self.dimension * [0]
        elif len(acceleration) != self.dimension:
            acceleration = np.pad(acceleration, (0, 3-len(velocity)),
                                  'constant', constant_values=0)

        # add the charge to every time step from time to time_size
        for i in range(0, time_N):
            new_charge = Charge(self, location, 0, velocity, acceleration)
            self.charges[i].append(new_charge)

        for i in range(time_N, self.time_size):
            new_charge = Charge(self,
                                location,
                                Q,
                                velocity,
                                acceleration,
                                mass)

            self.charges[i].append(new_charge)

        if print_all is not False:
            print('charge added with')
            print('location')
            print(self.charges[i][-1].location)
            print('velocity')
            print(self.charges[i][-1].velocity)
            print('acceleration')
            print(self.charges[i][-1].acceleration)
            print('charge')
            print(self.charges[i][-1].Q)
            print('time')
            print(time)
            print('')

        if count is True:
            # this must return the acual index of the new charge
            return(len(self.charges[0]) - 1)

    def modify_charge(self,
                      N,    # what charge you want to modify
                      time,  # at what time you want to modify the charge
                      velocity=False,  # the new velocity at time, time
                      acceleration=False,  # the new acceleration of the charge
                      location=False,  # this changes the charges locaiton
                      Q=False,  # this changes the coulombs of the charge
                      print_all=False):  # set true to print charge info

        # if the location of a charge is changed all informaiton
        # about velocity and acceleration is discarded and overwitten.

        time_N = int(np.rint(time / self.delta_t))

        if velocity is not False:
            # if len(velocity) != self.dimension:
            if len(velocity) != 3:
                velocity = np.pad(velocity, (0, 3 - len(velocity)),
                                  'constant', constant_values=0)
            self.charges[time_N][N].velocity = np.array(velocity)

        if acceleration is not False:
            # if len(acceleration) != self.dimension:
            if len(acceleration) != 3:
                acceleration = np.pad(acceleration, (0, 3 - len(acceleration)),
                                      'constant', constant_values=0)
            self.charges[time_N][N].acceleration = np.array(acceleration)

        if Q is not False:
            # notice that Q is not effecting velocity or location
            for t in range(time_N, len(self.charges)):
                # this loop will likely efect velocity and acceleration
                self.charges[t][N].Q = Q

        if location is not False:
            # if len(location) != self.dimension:
            if len(location) != 3:
                location = np.pad(location, (0, 3 - len(location)),
                                  'constant', constant_values=0)
            self.location = location

            for t in range(time_N, len(self.charges)):
                self.charges[t][N].location = np.array(location)

            temp_location = []
            temp_velocity = []
            temp_acceleration = []

            # since we modified the locaton of the charge update the gradient
            for t in range(len(self.charges)):
                temp_location.append(self.charges[t][N].location)

            temp_velocity = np.gradient(temp_location, axis=0)
            temp_acceleration = np.gradient(temp_velocity, axis=0)

            if print_all is not False:
                print('location')
                print(self.location)
                # print(temp_location)
                print('velocity')
                print(temp_velocity)
                print('acceleration')
                print(temp_acceleration)
                print('charge')
                print(Q)
                print('')

            # notice that the gradients that we just found
            # are the new velocity and acceleration
            for t in range(len(self.charges)):
                self.charges[t][N].velocity = temp_velocity[t]
                self.charges[t][N].acceleration = temp_acceleration[t]

    def Add_Current(self, location,  # location in acual distance
                    direction=[0, 0, 1],  # the direction the current points
                    Amps=0,  # the amps of the current
                    count=False,   # this will return the index of the current
                    print_all=False):  # print all information on the current

        #  this will add a current to the grid notice that a current is a class
        location = np.pad(location, (0, 3 - len(location)),
                          'constant', constant_values=0)

        direction = np.pad(direction, (0, 3 - len(direction)),
                           'constant', constant_values=0)

        if np.linalg.norm(direction) == 0:
            print("Current direction norm can't be zero")
            sys.exit()

        direction = (np.array(direction).astype(float) /
                     np.linalg.norm(direction))

        location = np.array(location)
        # insure that the current is in the middle of a grid squair
        location = (location / self.delta).astype(int) * self.delta

        for i in range(self.time_size):
            new_current = Current(self, location, direction, Amps)
            self.currents[i].append(new_current)

        if print_all is True:
            print('location = ' + str(location) + 'amps = ' + str(Amps) +
                  ' direction = ' + str(direction))

        if count is True:
            return(len(self.currents[0]) - 1)

    # this will allow you to change the amps and direction of a current
    # at this time you can't change the location
    def Modify_Current(self,
                       N,     # what current you want to modify
                       time=0,  # the time at which to make changes
                       direction=False,  # the new direction
                       amps=False,  # the new amps value
                       print_all=False):

        time_N = int(np.rint(time / self.delta_t))
        if print_all is True:
            print('the time of modification is = ' + str(time_N))
        if direction is not False:
            direction = np.pad(direction, (0, 3 - len(direction)),
                               'constant', constant_values=0)
            print(direction)

            # notice that the direction is eithor nomalized or zero
            if np.linalg.norm(direction) != 0:
                direction = (np.array(direction) /
                             np.linalg.norm(direction))

                if print_all is True:
                    print('Current modified')

                self.currents[time_N][N].direction = direction
            elif np.linalg.norm(direction) == 0:
                self.currents[time_N][N].direction = np.array(direction)

        if amps is not False:
            self.currents[time_N][N].amps = float(amps)

        if print_all is True:
            print('Modifying Current')
            print('New Time = ' + str(time_N))
            print('New Amps = ' + str(amps))
            print('New direction = ' + str(direction))
            print('poresent direction = ' +
                  str(self.currents[time_N][N].direction))
            print('')

    # rember a dipole is two charges seporated be some distance
    def Add_Dipole(self,
                   location,   # the locaiton to place the dipole
                   Q=1,        # the colembs of the dipole
                   dipole_vector=[1, 0, 0],  # the vector between the charges
                   separation=.01,  # the distance between the charges
                   velocity=False,  # the velocity of the dipole
                   acceleration=False,  # the acceleration of the dipole
                   time=0,     # the time at which to place the dipole
                   print_all=False,  # print the specs of the dipole
                   count=False):  # return the index of the new dipole

        time_N = int(np.rint(time / self.delta_t))

        dipole_vector = dipole_vector / np.linalg.norm(dipole_vector)
        dipole_vector = dipole_vector * separation

        dipole_moment = np.array(dipole_vector) * Q

        location = np.array(location)
        Q = Q
        dipole_vector = np.array(dipole_vector)
        separation = np.array(separation)

        if velocity is not False:
            velocity = np.array(velocity)
        else:
            velocity = np.array([0, 0, 0])

        if acceleration is not False:
            acceleration = np.array(acceleration)
        else:
            acceleration = np.array([0, 0, 0])

        self.dipole_vector = dipole_vector
        # add the charge to every time step from time to time_size
        charge_1 = self.Add_Charge(location +
                                   .5 *
                                   dipole_vector,   # locaiton of charge
                                   Q=Q,        # the colembs of the charge
                                   velocity=False,  # velocity of the charge
                                   acceleration=False,  # charge acceleration
                                   time=time,     # time to place the charge
                                   print_all=print_all,  # print charge specs
                                   count=True)  # return index of the charge

        charge_2 = self.Add_Charge(location -
                                   .5 *
                                   dipole_vector,  # locaiton of charge
                                   Q=-Q,        # the colembs of the charge
                                   velocity=False,  # velocity of the charge
                                   acceleration=False,  # charge acceleration
                                   time=time,     # time to place the charge
                                   print_all=print_all,  # print charge specs
                                   count=True)  # return, index of the charge
        # add the dipole to every time step from time to time_size
        for i in range(0, time_N):
            self.dipoles[i].append(Dipole(location,
                                          0,
                                          dipole_vector,
                                          separation,
                                          velocity,
                                          acceleration,
                                          charge_1,
                                          charge_2))

        for i in range(time_N, self.time_size):
            self.dipoles[i].append(Dipole(location,
                                          Q,
                                          dipole_vector,
                                          separation,
                                          velocity,
                                          acceleration,
                                          charge_1,
                                          charge_2))

        new_dipole = Dipole(location,
                            Q,
                            dipole_vector,
                            separation,
                            velocity,
                            acceleration,
                            charge_1,
                            charge_2)

        if count is True:
            return(len(self.dipoles[0]) - 1)

    def Modify_Dipole(self,
                      N,       # what Dipole you want to modify
                      location=False,   # the locaiton to move the dipole
                      Q=False,        # the new colembs of the dipole
                      dipole_vector=[0, 0, 0],  # new vector between charges
                      separation=False,  # the new distance between the charges
                      velocity=False,  # the new velocity of the dipole
                      acceleration=False,  # the new acceleration of the dipole
                      time=0,     # the time at which to modify the dipole
                      print_all=False):  # print the specs of the dipole

        time_N = int(np.rint(time / self.delta_t))

        if location is not False:  # location of the dipole cinter point
            self.dipoles[time][N].location = np.array(location)

        if Q is not False:  # colembs of the charges
            self.dipoles[time][N].Q = Q

        if separation is not False:  # distance between the charges
            self.dipoles[time][N].separation = separation

        dipole_vector = np.array(dipole_vector)
        if (dipole_vector != [0, 0, 0]).any():
            dipole_vector = dipole_vector / np.linalg.norm(dipole_vector)
            dipole_vector = dipole_vector * self.dipoles[time][N].separation
            self.dipoles[time][N].vector = np.array(dipole_vector)

        if velocity is not False:   # velocity of the dipole
            self.dipoles[time][N].velocity = np.array(velocity)

        if acceleration is not False:  # acceleration of the dipole
            self.dipoles[time][N].acceleration = np.array(accelaration)

        dipole = self.dipoles[time][N]
        self.dipoles[time][N].dipole_moment = (dipole.vector *
                                               dipole.Q *
                                               dipole.separation)

        self.modify_charge(self.dipoles[time][N].charges[0],  # modify charge
                           time*self.delta_t,  # when to modify the charge
                           velocity=False,  # the new velocity
                           acceleration=False,  # the new acceleration
                           location=(location +
                                     .5 *
                                     self.dipoles[time][N].vector),
                           Q=self.dipoles[time][N].Q,  # the new charge value
                           print_all=print_all)  # print info

        if print_all is True:
            print()

        self.modify_charge(self.dipoles[time][N].charges[1],  # modify charge
                           time*self.delta_t,  # when to modify the charge
                           velocity=False,  # the new velocity
                           acceleration=False,  # the new acceleration
                           location=(location -
                                     .5 *
                                     self.dipoles[time][N].vector),
                           Q=-self.dipoles[time][N].Q,  # the new charge value
                           print_all=print_all)  # print info

    def Add_Conductor(self,
                      location,   # the location of the conductor
                      materiel='Copper',  # the materiel of the conductor
                      print_all=False):

        location = np.array(location)
        if print_all is True:
            print('cunductor added')
            print('location = ' + str(location) +
                  ' materiel = ' + str(materiel))

        if len(location) != self.dimension:
            print('conductor grid dimension missmatch')
            sys.exit([2])

        location = np.array(location).astype(float)
        self.conductors.append(Conductor(materiel,
                                         location,
                                         self,
                                         print_all))

    def Add_Boundary(self, time_steps, location, E=[0, 0, 0], H=[0, 0, 0]):
        self.boundary[time_steps].append(Boundary(location, E, H))

    ''' note that permittivity can't be modified
        latter in the simulation it is fixed '''
    def add_Permittivity(self, location, material=False, permittivity=1):
        self.Permittivity.append(Permittivity(self,
                                              location,
                                              material,
                                              permittivity))

    ''' get_Permittivity is only for calculationg the polarization vector
    # this code should be made inoto a genorator '''

    ''' this next code needs tested but first
    i must finish permittivity calculations '''
    def get_Permittivity(self, location):

        for i in range(len(self.Permittivity)):
            return_valu = np.array(location / self.delta).astype(int)
#            if (self.Permittivity[i].index == np.array(location /
#                                                       self.delta).astype(int)).all():
            if (self.Permittivity[i].index == return_valu).all():
                return(self.Permittivity[i].R_Permittivity)
        return(1)


class Charge():  # this is used to define a charge on the grid
    def __init__(self, grid, location, Q, velocity, acceleration, mass):

        self.location = np.array(location).astype(float)
        self.Q = Q
        self.velocity = np.array(velocity)
        self.acceleration = np.array(acceleration)
        self.mass = mass


class Current():  # this is used to define a current on the grid
    def __init__(self, grid, location, direction, amps):

        self.location = np.array(location).astype(float)
        self.direction = np.array(direction)
        self.amps = float(amps)
        self.diff_t = 0


class Dipole():
    def __init__(self,
                 location,
                 Q,
                 vector,
                 separation,
                 velocity,
                 acceleration,
                 charge_1,
                 charge_2):
        self.location = np.array(location).astype(float)
        self.Q = Q
        self.vector = np.array(vector)
        self.separation = separation
        self.velocity = np.array(velocity)
        self.acceleration = np.array(acceleration)
        self.dipole_moment = vector / np.linalg.norm(vector) * Q * separation
        self.charges = [charge_1, charge_2]


class Boundary():
    # this is used to define a boundary condition
    def __init__(self, location, E, H):
        # E and H are boudary values of E and H needed
        self.E = np.array(E)
        self.H = np.array(H)
        # the location of the boundary value
        self.location = np.array(location)


class Permittivity():
    def __init__(self, grid, location, material=False, permittivity=1):
        self.location = np.array(location)
        self.index = tuple((self.location / grid.delta).astype(int))

        ''' later add diffalt values for some materiels'''
        if material is False:
            self.R_Permittivity = permittivity  # Relativ permittivity
        else:
            self.R_Permittivity = permittivity

            self.R_Permeability = permeability


class permittivity_normals():
    def __init__(self,
                 index,  # the index used to define the locaiton of the normal
                 vector):      # the vector that defines the normal
        self.index = tuple(np.array(index).astype(int))
        self.vector = vector


# test for approximate equality (for floating point types)
def arreqclose_in_list(myarr, list_arrays):
    return next((True for elem in list_arrays
                 if elem.size == myarr.size and
                 allclose(elem, myarr)), False)


# this is an old attempt to add Conductors
# and most likely will need heavy modifications
class Conductor():
    def __init__(self, materiel, location, grid, print_all=False):

        #  first add a current at the location of the conductor
        #  this will become the induction
        self.current = grid.Add_Current(location,
                                        direction=[1, 1, 1],
                                        Amps=0,
                                        count=True,
                                        print_all=False)

        # this is the term needed to calculate the induction currents
        # it is the electrokinetic field
        # self.EK_field = np.zeros( tuple(grid.time_size), dtype='complex')
        # self.EK_field = np.zeros( tuple(grid.time_size, 3), dtype='complex')

        self.EK_field = np.zeros(
            tuple(np.append(grid.time_size, [3])), dtype='complex')

        self.location = location
        # the location must be of the form [i, j, k]
        # with the same dimention as the grid

        # all ohm values are given in Ohm Meters at 20C
        if materiel == 'aluminium':
            self.Ohms = 2.8 * 10**-8
        elif materiel == 'Antimony':
            self.Ohm = 3.9 * 10**-7
        elif materiel == 'Bismuth':
            self.Ohm = 1.3 * 10**-6
        elif materiel == 'Brass':
            self.Ohm = 0.75 * 10**-7
        elif materiel == 'Cadmium':
            self.Ohm = 6 * 10**-8
        elif materiel == 'Cobalt':
            self.Ohm = 5.6 * 10**-8
        elif materiel == 'Copper':
            self.Ohm = 1.7 * 10**-8
        elif materiel == 'Gold':
            self.Ohm = 2.4 * 10**-8
        elif materiel == 'Carbon':  # (Graphite)
            self.Ohm = 1 * 10**-5
        elif materiel == 'Germanium':
            self.Ohm = 4.6 * 10**-1
        elif materiel == 'Iron':
            self.Ohm = 1.0 * 10**-7
        elif materiel == 'Lead':
            self.Ohm = 1.9 * 10**-7
        elif materiel == 'Manganin':
            self.Ohm = 4.2 * 10**-7
        elif materiel == 'Nichrome':
            self.Ohm = 1.1 * 10**-6
        elif materiel == 'Nickel':
            self.Ohm = 7 * 10**-8
        elif materiel == 'Palladium':
            self.Ohm = 1.0 * 10**-7
        elif materiel == 'Platinum':
            self.Ohm = 0.98 * 10**-7
        elif materiel == 'Quartz':
            self.Ohm = 7 * 1017
        elif materiel == 'Silicon':
            self.Ohm = 6.4 * 102
        elif materiel == 'Silver':
            self.Ohm = 1.6 * 10**-8
        elif materiel == 'Tantalum':
            self.Ohm = 1.3 * 10**-7
        elif materiel == 'Tin':
            self.Ohm = 1.1 * 10**-7
        elif materiel == 'Tungsten':
            self.Ohm = 4.9 * 10**-8
        elif materiel == 'Zinc':
            self.Ohm = 5.5 * 10**-8
        else:
            self.Ohm = 1.7 * 10**-8
