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
from .constants import *
from .linalg import *
from .extras import arreqclose_in_list

import numpy as np


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

        dimension = 3  # this is the dimenstion of the grid
        self.dimension = 3

        if type(constant_E) == list:
            self.constant_E = lambda location, time: constant_E
        else:
            self.constant_E = constant_E

        if type(constant_H) == list:
            self.constant_H = lambda location, time: constant_H
        else:
            self.constant_H = constant_H

        self.free_space = free_space
        print('the simulaiton is free space ' + str(self.free_space))

        # define a delta value for the grid in each direction
        if delta is False:
            self.delta = [1.0 for i in range(dimension)]
            print('recomended to give delta a valu diffalting to 1 meter')

        elif isinstance(delta, list) is True:
            if len(delta) != dimension:
                self.delta = delta * dimension
                print('Warning delta should be a size 3 array or a number')
                # break()
        else:
            delta = [float(delta)] * dimension

        self.delta = delta
        print('grid delta array is ' + str(self.delta))

        # define the size of the grid in each direction
        if size is False:
            self.size = [10] * 3
            print('grid size diffalting to ' + str(self.size))

        elif isinstance(size, list) is False:
            self.size = [float(size)] * dimensiomn

        elif len(size) == dimension:
            self.size = size
        elif len(size) == 1:
            self.size = size * 3
        else:
            print('error grid size must be a size 3 array or number')
            # break()

            print('the grid has size ' + str(size))

        # this is the amount of time the simulation will simulate
        self.time = time
        self.delta_t = delta_t  # this is the size of a time step in seconds
        print('the simulaiton will sumulate ' + str(self.time) + ' seconds')
        print('each time step will be ' + str(self.delta_t))

        self.time_size = int(self.time / self.delta_t)
        print('the simulation will be composed of ' + str(self.time_size) +
              ' time steps')

        self.shape = [int(size / delta) for
                      (size, delta) in
                      zip(self.size, self.delta)]

        self.shape_t = [self.time_size]
        self.shape_t.append(self.shape)

        # this will be the charges
        self.charges = [[] for i in range(self.time_size)]
        # this will be the currents
        self.currents = [[] for i in range(self.time_size)]
        self.dipoles = [[] for i in range(self.time_size)]

        '''# this will be the dipoles    ?I dont think this will be used?
        self.dipoles = [[] for i in range(self.time_size)]'''

        # this will hold the boundary conditions
        self.boundary = [[] for i in range(self.time_size)]

        # this is the conductors that are on the grid
        # notice that these don't change over the course of the simulation
        self.conductors = []

        # now define the relative permittivity grid this is also known as
        # as the dielectric constant

        ''' this will be the list of all Permittivity elements
        locations and thier dipoles '''
        self.Permittivity = []

        ''' this will be the list of all Permeability elements
        locations and thier dipoles '''
        self.Permeability = []

        # the grid layout is grid.grid[field][time axis][x, y, z axis]
        # the grid will hold the E, H and B fields as 3D vectors for all casses

        self.grid = {}
        self.grid['E'] = []
        self.grid['H'] = []
        self.grid['B'] = []
        self.grid['K'] = []

    def Add_Charge(self,
                   location,   # the locaiton to place the new charge
                   Q=(-1.6021766208) * 10**(-19),  # charge defalt electron
                   velocity=[0, 0, 0],  # the velocity of the charge
                   acceleration=[0, 0, 0],  # the acceleration of the charge
                   time=0,     # the time at which to place the charge
                   mass=9.10938356 * 10**(-31),  # kg defalt electron
                   print_all=False,  # print the specs of the charge
                   count=False,  # return the index of the new charge
                   dipole=False):
        # this will add a charge to the grid notice that a charge is a class
        if isinstance(location, list) is False:
            print('must use a list for location')
        while len(location) < 3:
            location.append(0)

        # time index to add the charg at
        time_N = int(round(time / self.delta_t))

        Q = float(Q)

        # add the charge to every time step from time to time_size
        for i in range(0, time_N):
            new_charge = Charge(self,
                                location,
                                0,
                                velocity,
                                acceleration,
                                dipole)

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

        time_N = int(round(time / self.delta_t))

        if velocity is not False:
            # self.charges[time_N][N].velocity = np.array(velocity)
            self.charges[time_N][N].velocity = velocity

        if acceleration is not False:
            # self.charges[time_N][N].acceleration = np.array(acceleration)
            self.charges[time_N][N].acceleration = acceleration

        if Q is not False:
            # notice that Q is not effecting velocity or location
            for t in range(time_N, len(self.charges)):
                # this loop will likely efect velocity and acceleration
                self.charges[t][N].Q = Q

        if location is not False:
            for t in range(time_N, len(self.charges)):
                self.charges[t][N].location = location

            temp_location = []
            temp_velocity = []
            temp_acceleration = []

            ''' since we modified the locaton of the charge
                             update the gradient '''
            for t in range(len(self.charges)):
                temp_location.append(self.charges[t][N].location)

            temp_velocity = np.gradient(temp_location, axis=0)
            temp_acceleration = np.gradient(temp_velocity, axis=0)

            ''' notice that the gradients that we just found
                are the new velocity and acceleration '''
            for t in range(len(self.charges)):
                self.charges[t][N].velocity = list(temp_velocity[t])
                self.charges[t][N].acceleration = list(temp_acceleration[t])

            if print_all is not False:
                print('location')
                print(location)
                print('time')
                print(time)
                print('time_step')
                print(time_N)
                print('velocity')
                print(temp_velocity)
                print('acceleration')
                print(temp_acceleration)
                print('charge')
                print(self.charges[t][N].Q)
                print('')

    def Add_Current(self, location,  # location in acual distance
                    direction=[0, 0, 1],  # the direction the current points
                    Amps=0,  # the amps of the current
                    count=False,   # this will return the index of the current
                    print_all=False):  # print all information on the current

        #  this will add a current to the grid notice that a current is a class
        for i in range(self.time_size):
            new_current = Current(self, location, direction, Amps)
            self.currents[i].append(new_current)

        if print_all is True:
            print('location = ' + str(location) + 'amps = ' + str(Amps) +
                  ' direction = ' + str(direction))

        if count is True:
            return(len(self.currents[0]) - 1)

    ''' this will allow you to change the amps and direction of a current
        at this time you can't change the location '''
    def Modify_Current(self,
                       N,     # what current you want to modify
                       time=0,  # the time at which to make changes
                       direction=False,  # the new direction
                       amps=False,  # the new amps value
                       print_all=False):

        time_N = int(round(time / self.delta_t))

        # notice that the direction is eithor nomalized or zero
        if direction is not False:
            norm = 0
            for d in direction:
                norm += norm + d**2
            norm = norm**.5
            if norm == 0:
                direction = direction
            else:
                direction = [d / norm for d in direction]

            self.currents[time_N][N].direction = direction

        if amps is not False:
            self.currents[time_N][N].amps = float(amps)

        if print_all is True:
            print('Modifying Current ' + str(N))
            print('the time of modification is = ' + str(time_N))
            print('New Time = ' + str(time_N))
            print('New Amps = ' + str(amps))
            print('New direction = ' + str(direction))
            print('poresent direction = ' +
                  str(self.currents[time_N][N].direction))
            print('')

    def Add_Dipole(self,
                   location,   # the locaiton to place the dipole
                   Q=1,        # the colembs of the dipole
                   dipole_vector=[1, 0, 0],  # the vector between the charges
                   separation=.01,  # the distance between the charges
                   velocity=[0, 0, 0],  # the velocity of the dipole
                   acceleration=[0, 0, 0],  # the acceleration of the dipole
                   time=0,     # the time at which to place the dipole
                   print_all=False,  # print the specs of the dipole
                   count=False):  # return the index of the new dipole
        ''' rember a dipole is two charges seporated be some distance '''

        time_N = int(time / self.delta_t)

        dipole_vector = [D / norm(dipole_vector) for D in dipole_vector]

        for i in range(len(dipole_vector)):
            dipole_vector[i] = dipole_vector[i] * separation

        dipole_moment = [D * Q for D in zip(dipole_vector)]

        self.dipole_vector = dipole_vector
        # add the charge to every time step from time to time_size
        charge_location = [L + .5 * D for (L, D) in zip(location,
                                                        dipole_vector)]

        charge_1 = self.Add_Charge(charge_location,   # locaiton of charge
                                   Q=Q,        # the colembs of the charge
                                   velocity=velocity,  # velocity
                                   acceleration=acceleration,  # acceleration
                                   time=time,     # time to place the charge
                                   print_all=print_all,  # print charge specs
                                   count=True,  # return index of the charge
                                   dipole=True)

        charge_location = [L - .5 * D for (L, D) in zip(location,
                                                        dipole_vector)]

        charge_2 = self.Add_Charge(charge_location,  # locaiton of charge
                                   Q=-Q,        # the colembs of the charge
                                   velocity=velocity,  # velocity
                                   acceleration=acceleration,  # acceleration
                                   time=time,     # time to place the charge
                                   print_all=print_all,  # print charge specs
                                   count=True,  # return, index of the charge
                                   dipole=True)

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
                      dipole_vector=False,  # new vector between charges
                      separation=False,  # the new distance between the charges
                      velocity=False,  # the new velocity of the dipole
                      acceleration=False,  # the new acceleration of the dipole
                      time=0,     # the time at which to modify the dipole
                      print_all=False):  # print the specs of the dipole

        time_N = int(np.rint(time / self.delta_t))

        if location is not False:  # location of the dipole cinter point
            self.dipoles[time][N].location = location

        if Q is not False:  # colembs of the charges
            self.dipoles[time][N].Q = Q

        if separation is not False:  # distance between the charges
            self.dipoles[time][N].separation = separation

        if dipole_vector is not False:
            if (dipole_vector != [0, 0, 0]):
                dipole_vector = [DP / norm(dipole_vector) *
                                 self.dipoles[time][N].separation
                                 for DP in dipole_vector]

                self.dipoles[time][N].vector = dipole_vector

        if velocity is not False:   # velocity of the dipole
            self.dipoles[time][N].velocity = velocity

        if acceleration is not False:  # acceleration of the dipole
            self.dipoles[time][N].acceleration = accelaration

        dipole = self.dipoles[time][N]

        self.dipoles[time][N].dipole_moment = [DP *
                                               self.dipoles[time][N].Q *
                                               self.dipoles[time][N].separation
                                               for DP in dipole.vector]

        charge_location = [L + .5 * V for (L, V) in
                           zip(location,
                               self.dipoles[time][N].vector)]

        self.modify_charge(self.dipoles[time][N].charges[0],  # modify charge
                           time*self.delta_t,  # when to modify the charge
                           velocity=False,  # the new velocity
                           acceleration=False,  # the new acceleration
                           location=charge_location,
                           Q=self.dipoles[time][N].Q,  # the new charge value
                           print_all=print_all)  # print info

        if print_all is True:
            print()

        charge_location = [L - .5 * V for (L, V)
                           in zip(location,
                                  self.dipoles[time][N].vector)]

        self.modify_charge(self.dipoles[time][N].charges[1],  # modify charge
                           time*self.delta_t,  # when to modify the charge
                           velocity=False,  # the new velocity
                           acceleration=False,  # the new acceleration
                           location=charge_location,
                           Q=-self.dipoles[time][N].Q,  # the new charge value
                           print_all=print_all)  # print info

    def Add_Conductor(self,
                      location,   # the location of the conductor
                      materiel='Copper',  # the materiel of the conductor
                      print_all=False):

        location = location
        if print_all is True:
            print('cunductor added')
            print('location = ' + str(location) +
                  ' materiel = ' + str(materiel))

        if len(location) != self.dimension:
            print('conductor grid dimension missmatch')
            sys.exit([2])

        self.conductors.append(Conductor(materiel,
                                         location,
                                         self,
                                         print_all))

#    def Add_Boundary(self, time_steps, location, E=[0, 0, 0], H=[0, 0, 0]):
#        self.boundary[time_steps].append(Boundary(location, E, H))

    ''' note that permittivity can't be modified
        latter in the simulation it is fixed
                 this may need changed '''
    def add_Permittivity(self,
                         location,
                         material=False,
                         permittivity=[1, 1, 1]):
        self.Permittivity.append(Permittivity(self,
                                              location,
                                              material,
                                              permittivity))

    ''' get_Permittivity is only for calculationg the polarization vector
                   this code should be made inoto a genorator? '''

    ''' this next code needs tested but first
    I must finish permittivity calculations '''
    def get_Permittivity(self, location):
        return_valu = [int(L / D) for (L, D)
                       in zip(location, self.delta)]
        for i in range(len(self.Permittivity)):
            if (self.Permittivity[i].index == return_valu):
                return(self.Permittivity[i].R_Permittivity)
        return([1, 1, 1])


class Charge():  # this is used to define a charge on the grid
    def __init__(self,
                 grid,
                 location,
                 Q,
                 velocity,
                 acceleration,
                 mass,
                 dipole=False):

        self.location = [float(location[0]),
                         float(location[1]),
                         float(location[2])]

        self.Q = float(Q)
        self.velocity = [float(velocity[0]),
                         float(velocity[1]),
                         float(velocity[2])]

        self.acceleration = [float(acceleration[0]),
                             float(acceleration[1]),
                             float(acceleration[2])]

        self.mass = float(mass)
        self.dipole = float(dipole)


class Current():  # this is used to define a current on the grid
    def __init__(self, grid, location, direction, amps):

        self.location = location
        self.direction = direction
        self.amps = float(amps)
        self.diff_t = [0, 0, 0]


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
        self.location = location
        self.Q = Q
        self.vector = vector
        self.separation = separation
        self.velocity = velocity
        self.acceleration = acceleration
        self.dipole_moment = [V / norm(vector) * Q * separation
                              for V in vector]
        self.charges = [charge_1, charge_2]


class Boundary():
    # this is used to define a boundary condition
    def __init__(self, location, E, H):
        # E and H are boudary values of E and H needed
        self.E = E
        self.H = H
        # the location of the boundary value
        self.location = location


class Permittivity():
    def __init__(self, grid, location, material=False, permittivity=[1, 1, 1]):
        self.location = location

        self.index = [l / D for (l, D) in zip(location, grid.delta)]
        ''' later add diffalt values for some materiels'''

        # Relativ permittivity
        if material is False:
            if isinstance(permittivity, list) is False:
                self.R_Permittivity = 3 * [permittivity]
            else:
                self.R_Permittivity = permittivity
        else:
            if isinstance(location, list) is False:
                self.R_Permittivity = 3 * [permittivity]
            else:
                self.R_Permittivity = permittivity
            self.R_Permittivity = permittivity


class permittivity_normals():
    def __init__(self,
                 index,  # the index used to define the locaiton of the normal
                 vector):      # the vector that defines the normal
        self.index = index
        self.vector = vector


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

        self.location = location

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
