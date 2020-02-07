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


using LinearAlgebra
using Debugger


C_0 = 299792458  # this is the speed of light in meters
K_e = 8.9875517873681764 * 10 ^ 9  # this is Coulomb's constant
E_0 = (4 * pi * K_e) ^ -1  # free space permittivity
U_0 = (C_0 ^ 2 * E_0) ^ -1   # free space permeability


function retardation(r, delta_t, c=C_0)::UInt64
    ret_val = ((r / c) / delta_t) + 1
    return(convert(UInt64, round(ret_val)))
end


function field_calculator(charges,
                          currents,
                          dipoles,
                          grid,
                          location,
                          time_0,
                          time_size)

    location = [Float64(location[1]),
                Float64(location[2]),
                Float64(location[3])]

    # this function maneges the calculaiton of the E and H fields

    # first we need some grids to hold the added field strenth
    field_E = []
    field_H = []
    field_B = []
    Electro_Kinetic = []
    #@enter for i in 1:time_size
    for i in 1:time_size
        push!(field_E, [0.0 ,0.0 ,0.0])
        push!(field_H, [0.0 ,0.0 ,0.0])
        push!(field_B, [0.0 ,0.0 ,0.0])
        push!(Electro_Kinetic, [0.0 ,0.0 ,0.0])
    end

    # find the field resulating form the currents at time "time"
    for current in (currents)

        # find the distance from the current to the locaiton of intrest
        R = location - current.location
        r = (dot(R, R)^.5)
        #r = Flote64(dot(R, R)^.5)

        if r != 0
            # calculate the time at which to add the new field quantity
            time = time_0 + retardation(r, grid.delta_t)

            # if time is less then the length of the simulaiton
            # calculate the added field
            Permittivity = grid.get_Permittivity
            if time < time_size + 1
                # field_E[time] = (field_E[time] +
                                 #E_field_currents(current.diff_t,
                F_E = (E_field_currents(current.diff_t,
                                        current.location,
                                        current.amps,
                                        current.direction,
                                        r,
                                        R,
                                        #Permittivity(location)) )])
                                        Permittivity(location)) )
                field_E[time] = F_E[1] + field_E[time]
                #Electro_Kinetic[time] = F_E[2] + Electro_Kinetic[time]
                Electro_Kinetic[time] = F_E[2] + Electro_Kinetic[time]

                #H_F, B_F = Magnetic_field_currents(current,
                #                                   r,
                #                                   R,
                #                                   grid)

                H = Magnetic_field_currents(current,
                                                   r,
                                                   R,
                                                   grid)

                H_F = H[1]
                B_F = H[2]

                field_H[time] = field_H[time] + H_F
                field_B[time] = field_B[time] + B_F
            end
        end
    end


    # time = time_0
    # now repet the above for charges
    for charge in charges
        # calcualte the distance to the charge
        R = location - charge.location
        r = Float64(dot(R, R)^.5)

        # find the field resulating from the charge at time "time"
        if r != 0
            time = time_0 + retardation(r, grid.delta_t)

            # if time is less then the length of the
            # simulaiton calculate the added field
            if time < time_size + 1
                field_E[time] = (field_E[time] +
                                     E_field_charges(charge,
                                     # E_field_charges(charge,
                                                     charge.velocity,
                                                     charge.acceleration,
                                                     charge.Q,
                                                     location,
                                                     r,
                                                     R,
                                                     grid,
                                                     grid.get_Permittivity(
                                                     charge.location)))

##
                M_F_C = Magnetic_field_charges(charge, location, r, R, grid)

                H_F = M_F_C[1]
                B_F = M_F_C[2]

                field_H[time] = field_H[time] + H_F
                field_B[time] = field_B[time] + B_F
            end
        end
    end

    return(field_E, field_H, field_B, Electro_Kinetic)
end


function E_field_currents(current_diff_t,
                         location,
                         amps,
                         direction,
                         r,
                         R,
                         Permittivity) # ::Array{Float64}

    # find the E field do to currents on the grid
    # see equation 2-2.12 of
    # electromagnetic retardation and theory of relativity

    S = 1/(4 * pi * E_0 * C_0)
    scale = [0.0, 0.0, 0.0]
    for i in 1:3
        scale[i] = S / Permittivity[i]
    end

    # first find the term dependent on current.amp
    # this is the electrostatic field

    electrostatic = [0.0, 0.0, 0.0]

    for i in 1:3
        electrostatic[i] = - scale[i] * amps * 2 * direction[i] / r^2
    end

    scale_R = [0.0 0.0 0.0]
    for i in 1:3
        scale_R[i] = scale[i] * R[i]
    end

    electrostatic_2 = (scale_R * 2 * amps * dot(direction, R) * 2 / (r^4))

    for i in 1:3
        electrostatic[i] = electrostatic[i] + electrostatic_2[i]
    end

    scale_diff = [0.0 0.0 0.0]

    # now find the term dependent on current_diff

    # electrokinetic = scale_R * dot(current_diff_t, R) / (r^3  * C_0)
    Current_Kinetic = scale_R * dot(current_diff_t, R) / (r^3  * C_0)

    for i in 1:3
        scale_diff[i] = scale[i] * current_diff_t[i]
    end

    scale_diff = scale_diff / (r * C_0)

    Current_Kinetic = Current_Kinetic - scale_diff

    field = [0.0, 0.0, 0.0]
    electrokinetic = [0.0, 0.0, 0.0]
    for i in 1:3
        electrokinetic[i] = - scale_diff[i]
        field[i] = electrostatic[i] + Current_Kinetic[i]
    end

    return (field, electrokinetic)

end


function Magnetic_field_currents(current, r, R, grid)# ::Array{Float64}i

    # this function finds the H field do to currents
    # see equation 2-2.5 of
    # electromagnetic retardation and theory of relativity

    # calculate the H field resulting from the currents
    # scale is a factor resulting from geometry
    # and a good place to concider for permibility
    scale = 1 / (4 * pi)


    H_field_amps = scale * current.amps * cross(current.direction, R) * 2 / (r^2)

    # first find the terms do to steady currents

    # now find the term dependent on current.diff_t
    # this is the electrokinetic field

     H_field_diff_t = scale * cross(current.diff_t, R) / (r^2 * C_0)


    # find the total H field
    H_field = H_field_amps + H_field_diff_t
    B_field = H_field * U_0

    return(H_field, B_field)
end



function E_field_charges(charge,
                         velocity,
                         acceleration,
                         Q,
                         location,
                         r,
                         R,
                         grid,
                         Permittivity)::Array{Float64}
        # this function finds the E field resulting from charges
    # see equation 4-4.34 of
    # electromagnetic retardation and theory of relativity

    # chang this to a sum of equations 4-1.11 and 4-4.32
    # this will only requier 2 dot products

    #velocity = parse.(Float64, velocity)

    scale = [.0, .0, .0]
    for i in 1:3
        scale[i] = 1 / (4 * pi * E_0 * Permittivity[i])
    end

    velocity = velocity / C_0

    velocity_r = velocity / r

    dot_velocity_term = dot(R, velocity_r)

    common_factor = 1 / (1 - dot_velocity_term)

    dot_acceleration_term = dot(R, acceleration)

    E_A_scale = (charge.Q * scale * common_factor^3 * dot_acceleration_term * (1/(r^2)))

    E_A = [0.0, 0.0, 0.0]

    for i in 1:3
        E_A[i] = E_A_scale[i] * (R[i] / r + velocity[i]) + acceleration[i] / r
    end

    E_B_scale = charge.Q * scale * common_factor^3 * (1 - dot(velocity, velocity)) / (r^3)
    E_B = [0.0, 0.0, 0.0]
   for i in 1:3
       E_B[i] = E_B_scale[i] * (R[i] - r * velocity[i])
   end

    E_field = [0.0, 0.0, 0.0]

    for i in 1:3
        E_field[i] = E_A[i] + E_B[i]
    end
    return(E_field)
end



function Magnetic_field_charges(charge, location, r, R, grid)
    # this function finds the H field resulting from charges
    # see equation 4-5.10 of
    # electromagnetic retardation and theory of relativity

    # scale is a factor resulting from geometry
    # and a good place to conscider for permibility
    scale = 1 / (4 * pi)

    velocity = charge.velocity

    acceleration = charge.acceleration

    # first find the most used term in the remaning equations
    most_used_factor = 1/(r * (1 - dot(R, velocity) / (r * C_0)))

    # find the common term
    common_factor = charge.Q * scale * most_used_factor^2

    # next find the velocity factor
    velocity_factor = (1 - dot(velocity, velocity) / C_0^2
                       + dot(R, acceleration) / C_0^2)

    field = cross(acceleration / C_0 + velocity * most_used_factor * velocity_factor, R)

    # find the total H field
    H_field = common_factor * field

    B_field = U_0 * H_field
    return(H_field, B_field)
end
