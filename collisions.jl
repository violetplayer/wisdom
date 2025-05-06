##==================================================================================================================
## Ionizing collisions from Janev (1993)
## Some formalism borrowed from Scott Nicks' code, ROTJ: https://git.trialphaenergy.com/theory/rotj
##==================================================================================================================

function relative_energy(m1::Float64, m2::Float64, v1::Array{Float64}, v2::Array{Float64})::Float64
    0.5 * (m1 * m2) / (m1 + m2) * mp * norm(v1 - v2)^2
end

function ST_cold_integral(x1::Array{Float64}, x2::Array{Float64}, K_o_m::Float64, state::PlasmaState)::Float64
    """
    Computes the total (cold, T_target = 0) fast neutral shinethrough for a single pseudo-particle from x1 to x2:
    
    ST(K/m,T=0) = (σ_CX(K/m) + σ_IZ(K/m)) * ∫n_i(t)dt + σ_EZ(K/m) * ∫n_e(t)dt

    Inputs:
        x1    : Initial position
        x2    : Final position
        K     : fast neutral energy [keV]
        m     : fast neutral mass [amu]
        state : PlasmaState struct
    """

    #### get integral distance and normed direction
    D = norm(x2 .- x1)
    n = (x2 .- x1) ./ D

    #### get all cold cross sections
    σE = σEZ_Janev(K_o_m)
    σI = σIZ_Janev(K_o_m) + σCX_Janev(K_o_m)

    #### compute + return value of integral
    quadgk(t -> σE * ne(rproj(x1,n,t), state) + σI * ni(rproj(x1,n,t), state), 
    0.0, D, rtol=1e-2)[1]
end

function σCX_fast_beam(v::Array{Float64}, beam::NeutralBeam, component_index::Integer)::Float64
    σCX_Janev(relative_energy(m_fi, m_fn, v, beam.speeds[component_index] * beam.grid_normal))
end



# function fast_thermal_CX(K::Float64, T::Float64)::Float64
#     throw(ErrorException("not built yet"))
# end

# function τse(r::Float64, K::Float64, state::PlasmaState)::Float64
#     6.27e8 * Te(r, state)^1.5 / (3 * ne(r, state) * log_Λ(r, "ei")) * log(1 + (K / (14.8 * Te(r, state)))^1.5)
# end


# function nσ_path_integral(x1::Array{Float64}, x2::Array{Float64})::Float64
#     D = sqrt(sum((x2 .- x1).^2))
#     n = (x2 .- x1) ./ D

#     function __integrand(t::Float64)::Float64
#         _r = rproj(x1, n, t)
#         return ne(_r) * σ... 
#     end

#     quadgk(f, 0.0, D, rtol=1e-3)[1]
# end

# function collrad(x::Array{Float64}, v::Array{Float64}, plasma::PlasmaState)
#     throw(ErrorException("Not built."))
# end

# function beam_target_reactivity(v::Array{Float64}, T::Float64, σ::AbstractFunction)
# end

# function relative_speed(v1::Array{Float64}, v2::Array{Float64})::Float64
#     sqrt(norm(v1)^2 + norm(v2)^2 - 2 * dot(v1,v2))
# end

# function relative_speed(v1::Float64, v2::Float64, θ::Float64)::Float64
#     sqrt(v1^2 + v2^2 - 2 * v1 * v2 * cos(θ))
# end

# function get_σ0(E::Float64, T::Float64)

#     vb = sqrt(2 * E * 1e3 * q0 / (m * mp))
#     vi = sqrt(2 * T * q0 / (m * mp))
    
#     # Collision energy function
#     function collision_energy_eV(i::Int64,j::Int64,m_target_spec::Float64,vt_spec::Float64)
#         return 0.5 * m_target_spec * (v_arr[i]^2*vt_spec^2 + state.vb^2 - 2 * state.vb * v_arr[i]*vt_spec * cos(theta_arr[j])) / eV_to_J
#     end
       
#     # Initialize cross section dictionary
#     sigma_list = Dict()
    
#     # Loop over beam ion birth processes
#     for code in sigma_codes

#         # Generate 2D integrand array of cross sections
#         if code == "EZ"
#             # Electron-impact ionization
#             m_target = me
#             vt = sqrt(state.T_thermal * eV_to_J / m_target)
#             sigma_array = [σEZ_Janev(collision_energy_eV(i,j,m_target,vt)) for i in 1:v_res for j in 1:theta_res]
#         elseif code == "IZ"
#             # Ion-impact ionization
#             m_target = mp * m_target_amu
#             vt = sqrt(state.T_thermal * eV_to_J / m_target)
#             sigma_array = [σIZ_Janev(collision_energy_eV(i,j,m_target,vt)/1000/m_target_amu) for i in 1:v_res for j in 1:theta_res]
#         elseif code == "CX"
#             # Charge exchange
#             m_target = mp * m_target_amu
#             vt = sqrt(state.T_thermal * eV_to_J / m_target)
#             sigma_array = [σCX_Janev(collision_energy_eV(i,j,m_target,vt)/1000/m_target_amu) for i in 1:v_res for j in 1:theta_res]
#         end

#         # Generate 2D integrand array of Maxwellian distribution function
#         func_array = [exp(-v_arr[i]^2/2) * v_arr[i]^2 * sin(theta_arr[j]) for i in 1:v_res for j in 1:theta_res]

#         # Get total integrand array
#         integrand_array = sigma_array .* func_array

#         # Compute integral and write to dictionary
#         sigma_list[code] = sum(integrand_array) * theta_inc * v_inc / sqrt(2*pi)
#     end
    
#     # Get total cross section and normalize
#     σ0 = sum(values(sigma_list)) * 1.0E+19
    
#     # Return sigma values in prescribed order
#     return σ0, Float32(sigma_list["CX"]), Float32(sigma_list["EZ"]), Float32(sigma_list["IZ"])
# end

