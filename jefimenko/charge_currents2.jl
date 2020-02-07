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

# pi = 3.141592653589793
C_0 = 299792458  # this is the speed of light in meters
K_e = 8.9875517873681764 * 10 ^ 9  # this is Coulomb's constant
E_0 = (4 * pi * K_e) ^ -1  # free space permittivity
U_0 = (C_0 ^ 2 * E_0) ^ -1   # free space permeability


function retardation(r, delta_t, c=C_0)::UInt64
    ret_val = ((r / c) / delta_t)
    return(convert(UInt64, round(ret_val)))
end


function E_field_currents(current_diff_t,
                         location,
                         amps,
                         direction,
                         r,
                         R,
                         Permittivity)::Float64

    # find the E field do to currents on the grid
    # see equation 2-2.12 of
    # electromagnetic retardation and theory of relativity

    scale = 1/(4 * pi * E_0 * C_0 * Permittivity)

    # first find the term dependent on current.amp
    # this is the electrostatic field

    electrostatic = - scale .* (amps * 2 / r^2) # .* direction # .* 2 ./ r^2)

    for i in 1:3
        electrostatic[i] = electrostatic[i] * direction[i]
    end

    scale_R = [0.0 0.0 0.0]
    for i in 1:3
        scale_R[i] = scale[i] * R[i]
    end

    electrostatic = electrostatic + (scale_R * 2 * amps * dot(direction, R) * 2 / (r^4))

    scale_diff = [0.0 0.0 0.0]

    # now find the term dependent on current_diff

    for i in 1:3
        scale_diff[i] = scale[i] * current_diff_t[i]
    end

    electrokinetic = scale_R * dot(current_diff_t, R) / (r^3  * C_0)

    scale_diff = scale_diff / (r * C_0)
    
    electrokinetic = electrokinetic  - scale_diff

    field = electrostatic + electrokinetic

    return field 

end







function Magnetic_field_currents(current, r, R, grid)::Float64

    # this function finds the H field do to currents
    # see equation 2-2.5 of
    # electromagnetic retardation and theory of relativity

    # calculate the H field resulting from the currents
    # scale is a factor resulting from geometry
    # and a good place to concider for permibility
    scale = 1 / (4 * pi)

    
    H_field_amps = scale * current.amps * cross(current.direction, R) * 2 / (r^2)

    # first find the terms do to steady currents
#    H_field_amps = [scale * current.amps * Cross * 2 / (r**2) for Cross in
#                    cross(current.direction, R)]

    # now find the term dependent on current.diff_t
    # this is the electrokinetic field

     H_field_diff_t = scale * cross(current.diff_t, R) / (r^2 * C_0)

#    H_field_diff_t = [scale * Cross / (r**2 * C_0)
#                      for Cross in cross(current.diff_t, R)]

    # find the total H field
    H_field = H_field_amps + H_field_diff_t
    B_field = H_field * U_0
    
#    H_field = [H_f + H_d for (H_f, H_d) in zip(H_field_amps, H_field_diff_t)]

#    B_field = [U_0 * H for H in H_field]
    return(H_field, B_field)
end


    
function E_field_charges(charge::Array{Float64},
                         velocit::Array{Float64},
                         acceleration::Float64}} ,
                         Q::Float64,
                         location::Float64,
                         r::Float64,
                         R::Float64,
                         grid::Float64,
                         Permittivity::Float64)::Float64
        # this function finds the E field resulting from charges
    # see equation 4-4.34 of
    # electromagnetic retardation and theory of relativity

    # chang this to a sum of equations 4-1.11 and 4-4.32
    # this will only requier 2 dot products

    scale = 1/(4 * pi * E_0 * C_0 * Permittivity)
    
    # scale = [1 / (4 * pi * E_0 * perm) for perm in Permittivity]

    velocity = velocity / C_0
    # velocity = [v / C_0 for v in velocity]

    velocity_r = velocity / r
    # velocity_r = [v / r for v in velocity]

    dot_velocity_term = dot(R, velocity_r)
    # dot_velocity_term = dot(R, (velocity_r))

    common_factor = 1 / (1 - dot_velocity_term)
    # common_factor = 1 / (1 - dot_velocity_term)

    dot_acceleration_term = dot(R, acceleration)
    # dot_acceleration_term = dot(R, acceleration)

    E_A_scale = ((charge.Q * scale * common_factor^3 * dot_acceleration_term * (1/(r*r))))
    # E_A_scale = ([(charge.Q * sc *
    #                common_factor * common_factor
    #                * (common_factor * dot_acceleration_term * (1/(r*r))))
    #               for sc in scale])

    E_A = [0, 0, 0]

    for i in 1:3
        E_A[i] = E_A_scale[i] * (R[i] / r + velocity[i]) + acceleration[i] / r
    end

#    E_A = [E_A_sc * (R_e / r + v) + a / r
#           for E_A_sc, R_e, v, a in zip(E_A_scale,
#                                        R,
#                                        velocity,
#                                        acceleration)]

    # E_A = [E_A_sc * (R_e / r + v) + a / r
    #        for E_A_sc, R_e, v, a in zip(E_A_scale,
    #                                     R,
    #                                     velocity,
    #                                     acceleration)]

    E_B_scale = charge.Q * scale * common_factor^3 * (1 - dot(velocity, velocity)) / (r^3)
    # E_B_scale = [charge.Q * sc * common_factor * common_factor * common_factor
    #              * (1 - dot(velocity, velocity)) / (r*r*r) for sc in scale]

    E_B = [0, 0, 0]
    for i in 1:3
        E_B[i] = E_B_scale[i] * (R[i] - r * velocity[i])
    end
    # E_B = [E_B_sc * (R_e - r * v) for
    #        E_B_sc, R_e, v in zip(E_B_scale, R, velocity)]

    #E_field = [0, 0, 0]
    #for i in 1:3
    E_field[i] = E_A[i] + E_B[i]
    #end
    # E_field = [E_a + E_b for E_a, E_b in zip(E_A, E_B)]
    return(E_field)
end
