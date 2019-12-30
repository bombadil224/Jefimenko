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

# this file will simulate the boundary conditions

from .constants import *
from .charge_current_integrals import *

import numpy as np


def genorate_regions(grid):
    ''' this must genorate all regions and hold (pointers?)
    to the values used on the boundary '''
    return([], [])


def boundary_simulation(region,  # the region being simulated
                                 # its acual locations in the grid.
                        boundaries,  # the boundary values of
                                     # the region being simulated.
                        time_0,   # the time of the boundary conditions
                        E_field,
                        H_field,
                        grid):

    # this will simulate the fields do to the boudary of the grid
    E_wave = []
    H_wave = []

    for i in range(len(boundaries[time_0])):
        # pdb.set_trace()
        """ At some point it must be detormined if the E or H field
        points along the pointng vector, if it is then this part will
        travel in the direction of the pointing vector without spreading out.
        """
        # calculate the normalized pointing vector
        pointing = np.cross(boundaries[time_0][i].E, boundaries[time_0][i].H)

        # now break up the E and H fields into its pointing vector
        # and not pointing vector components
        if (pointing != [0, 0, 0]).any():
            pointing = pointing / np.linalg.norm(pointing)

            pointing_E = pointing * np.abs(boundaries[time_0][i].E)
            pointing_H = pointing * np.abs(boundaries[time_0][i].H)

            # now find the part of E and H that will disperce over distance
            for j in range(len(boundaries[time_0])):
                analize_wave_E = []
                analize_wave_H = []
                for time in range(len(boundaries)):

                    E_wave.append(boundaries[time][j].E - pointing *
                                  np.dot(boundaries[time][j].E, pointing))

                    H_wave.append(boundaries[time][j].H - pointing *
                                  np.abs(boundaries[time][j].H, pointing))

                    # now to find the fft all that we will need
                    # is the norm of these vectors

                    analize_wave_E.append(np.linalg.norm(E_wave))
                    analize_wave_H.append(np.linalg.norm(H_wave))

                # now find the wavelengths of the incoming wave that
                # will disperce this must return the wavelengths and
                # the amplitudes of the given wave

                Wave, Amplitude = fft_analysis(analize_wave_E, grid)
                WaveLength = Wave
                NormAmplitude = Amplitude

        else:
            # the part of E and H that will not disperce over distance
            # they are static in nature
            if np.linalg.norm(boundaries[time_0][i].E) != 0:
                pointing_E = (boundaries[time_0][i].E /
                              np.linalg.norm(boundaries[time_0][i].E))
            else:
                pointing_E = np.array([0, 0, 0])

            if np.linalg.norm(boundaries[time_0][i].H) != 0:
                pointing_H = (boundaries[time_0][i].H /
                              np.linalg.norm(boundaries[time_0][i].H))
            else:
                pointing_H = np.array([0, 0, 0])

        if (pointing != [0, 0, 0]).any():
            for location in region:
                # calculate the distance to the boundary condition
                R = location - boundaries[time_0][i].location
                r = np.linalg.linalg.norm(R)

                # calculate the retardation time for the field
                time = time_0 + retardation(r, grid)
                if time < len(E_field):
                    # calculate the field
                    if r != 0:
                        for k in range(len(WaveLength)):
                            diffraction_angle = (np.arccos(np.dot(r, E_wave) /
                                                           (np.linalg.norm(r) *
                                                            np.linalg.norm(E_wave))))

                            E_field[time][tuple(location)] = (
                                E_field[time][tuple(location)] +
                                diffraction(diffraction_angle[k],
                                            E_wave[k] *
                                            NormAmplitude[k],
                                            grid.delta[0],
                                            WaveLength[k]))

                            diffraction_angle = (np.arccos(np.dot(r, H_wave) /
                                                          (np.linalg.norm(r) *
                                                           np.linalg.norm(H_wave))))

                            H_field[time][tuple(location)] = (
                                H_field[time][tuple(location)] +
                                diffraction(diffraction_angle[k],
                                            NormAmplitude[k],
                                            grid.delta[0],
                                            WaveLength[k]))

        location = boundaries[time_0][i].location

        if (pointing_E != [0, 0, 0]).any() or (pointing_H != [0, 0, 0]).any():
            while region_test(location, region) == 1:    # coulemb potential

                R = location - boundaries[time_0][i].location
                # calcualte the distance to the boundary condition
                r = np.linalg.norm(R)

                # calculat the retardation time for the field
                time = time_0 + retardation(r, grid)
                if time < len(E_field):
                    if r != 0:
                        loc = tuple(location.astype(int))

                        E_field[time][loc] = (E_field[time][loc] +
                                              pointing_E / r**2)
                        H_field[time][loc] = (H_field[time][loc] +
                                              pointing_H / r**2)
                location = location + pointing_E
    return(E_field, H_field)


def region_test(location, region):
    for test in region:
        if (test == location.astype(int)).all():
            return(1)
    return(0)


def diffraction(a,      # differaction angle, angle between source angle and r
                I_0,    # intensity at source
                d,      # slit width
                WaveLength):

    ''' this will be the first mothed of calculationg difforaction patterns
    and will just use the equations for diforaction form a singel slit.
    this assumes that the distance between the screen and the aperture is
    veary large compared with the slit width. later this should be exspanded
    to use huygens principle directly. '''

    angle_dependance = np.pi * d * np.sin(a) / WaveLength

    I_a = (np.sin(angle_dependance))**2 / (angle_dependance ** 2)
    return(I_a)


def fft_analysis(wave, grid, c=C_0):
    # first finde the frequencies of the waves
    fft_wave = np.fft.fft(wave)
    # next find the corisponding frequencies

    df = 1 / grid.delta_t      # number of measurments per second

    # now transform the fftfreq function to hertz
    freqs = np.fft.fftfreq(len(fft_wave)) * len(fft_wave) * df
    # now find the wavelengthes
    length = []
    amplitude = []
    for k in range(len(freqs)):
        if freqs[k] != 0:
            length.append(c / freqs[k])
            amplitude.append(np.abs(fft_wave[k]) / len(wave))

    # now find the percintage of the wave with eatch frequencie
    total_amp = sum(amplitude)
    return(length, amplitude)
