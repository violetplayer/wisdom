##==================================================================================================================
## Method of Effective Potentials - 1D
##==================================================================================================================

const r_mid = get_midpoints(RSgrid.r) # handy constant

#-----------------
# Static ion properties
#-----------------

#### does this orbit-space point represent a physical particle?
function ion_exists(Pθ::Float64, H⊥::Float64, state::FRCState, RSgrid::RealSpaceGrid)::Bool
    if H⊥ >= Uϕ_max
        V⊥_0 = sqrt(2 / (m_fi * mp)) * sqrt(H⊥ - q_fi * q0 * ϕ(rw, state))
        V⊥_w = sqrt(2 / (m_fi * mp)) * sqrt(H⊥ - q_fi * q0 * ϕ(state.r0, state))

        return (Pθ > -m_fi*state.r0*V⊥_0 + q_fi*q0*ψ(state.r0, state)) && 
        (Pθ < m_fi*rw*V⊥_w + q_fi*q0*ψ(rw, state)) && 
        H⊥ >= 1.05 * minimum([U_eff(Pθ, ir, state) for ir in RSgrid.r])
    else
        return false
    end
end

function ion_exists(Pθ::Float64, H⊥::Float64, state::ECM, RSgrid::RealSpaceGrid)::Bool
    return true
    #return (H⊥ >= Uϕ_max) && (H⊥ >= 1.05 * minimum([U_eff(Pθ, ir, state) for ir in RSgrid.r]))
end


#### is the orbit confined?
function orbit_is_confined(Pθ::Float64, H⊥::Float64, state::PlasmaState)::Bool
    Uw = U_eff(Pθ, rw, state)
    return H⊥ < Uw
end

#### find turning points as roots in U_eff given H_⊥
function get_turning_points(Pθ::Float64, H⊥::Float64, state::PlasmaState, RSgrid::RealSpaceGrid)::Array{Float64}

    #### find roots
    rt_init = guess_roots(r_mid, [1 - U_eff(Pθ,x,state)/H⊥ for x in RSgrid.r])
    rts = zeros(length(rt_init))
    for (i,rt) in enumerate(rt_init)
        if rt == r_mid[1]
            rts[i] = rt # no change - fix for near-0 momentum root finding failing
        else
            rts[i] = find_zero(x -> 1 - U_eff(Pθ,x,state)/H⊥, rt)
        end
    end

    #### fixes for missing points
    Uw = U_eff(Pθ, rw, state)
    if length(rts) == 0
       rts = [RSgrid.r[1], rw]
    elseif length(rts) == 1 || length(rts) == 3
        if H⊥ <= Uw # missing interior point
            pushfirst!(rts, RSgrid.r[1])
        elseif H⊥ > Uw # missing exterior point (unconfined orbit, set outer turning point to wall radius)
            push!(rts, rw)
        end
    end

    #### adjust roots to make sure we are 'inside' the turning point on the U_eff curve (keeps nasty errors away)
    out = zeros(length(rts))
    for i=eachindex(rts)
        y_i = 1 - U_eff(Pθ,rts[i],state)/H⊥
        if (y_i <= 0) && (rts[i] < rw)
            exit = false
            j = 1
            while (~exit) & (j < 100)
                #### test...
                δr = j*(RSgrid.Δr / 100)
                δL = (1 - U_eff(Pθ,rts[i]-δr,state)/H⊥)
                δR = (1 - U_eff(Pθ,rts[i]+δr,state)/H⊥)
                if (y_i - δL) < 0
                    if δL > 0
                        out[i] = rts[i] - δr
                        exit = true
                    end
                else
                    if δR > 0
                        out[i] = rts[i] + δr
                        exit = true
                    end
                end
                j += 1
            end
        else
            out[i] = rts[i] # no modification
        end
    end

    return out
end

#### compute the radial bounce time
function get_τr(Pθ::Float64, H⊥::Float64, r1::Float64, r2::Float64, state::PlasmaState)::Float64
    norm = sqrt(2 * m_fi * mp / H⊥)
    I = quadgk(x -> 1 / sqrt(1 - U_eff(Pθ,x,state)/H⊥), r1, r2, rtol=1e-3)[1]
    norm * I
end

#### compute the azimuthal transit time
function get_τθ(Pθ::Float64, H⊥::Float64, τr::Float64, r1::Float64, r2::Float64, state::PlasmaState)::Float64
    _vθ(x) = (Pθ - q_fi * q0 * ψ(x, state)) / (m_fi * mp * x)
    _v⊥(x) = sqrt(2 * (H⊥ - q_fi * q0 * ϕ(x, state)) / (m_fi * mp))
    I = quadgk(x -> (_vθ(x) / (x * sqrt(_v⊥(x)^2 - _vθ(x)^2))), r1, r2, rtol=1e-2)
    2 * π * τr / abs(I)
end

#### compute radial action J_r
function get_ract(Pθ::Float64, H⊥::Float64, r1::Float64, r2::Float64, state::PlasmaState)::Float64
    I = quadgk(x -> sqrt(1 - U_eff(Pθ,x,state)/H⊥), r1, r2, rtol=1e-2)[1]
    sqrt(m_fi * mp * H⊥ / (2 * π^2)) * I
end

function map_J0(OSgrid::OrbitSpaceGrid, RSgrid::RealSpaceGrid, state::PlasmaState)::Array{Float64}
    J0 = zeros(length(OSgrid.Pθ), length(OSgrid.H))
    for (i,p) in enumerate(OSgrid.Pθ)
        for (j,h) in enumerate(OSgrid.H)
            if ion_exists(p*P0, h*eV_to_J*1e3, state, RSgrid)
                rts = get_turning_points(p*P0, h*eV_to_J*1e3, state, RSgrid)
                J0[i,j] = get_ract(p*P0, h*eV_to_J*1e3, rts[1], rts[2], state)
            end
        end
    end
    return J0
end

function get_orbit_bounds(state::PlasmaState, OSgrid::OrbitSpaceGrid)
    Pmin = [-m_fi * mp * state.r0 * sqrt(2 * OSgrid.H[i] * 1e3 * eV_to_J / (m_fi * mp)) + q_fi * q0 * ψ(state.r0, state) for i=eachindex(OSgrid.H)]
    Pmax = [m_fi * mp * rw * sqrt(2 * OSgrid.H[i] * 1e3 * eV_to_J / (m_fi * mp)) + q_fi * q0 * ψ(rw, state) for i=eachindex(OSgrid.H)]
    Pstar = q_fi * q0 * ψ(state.r0, state)
    Psstar = [q_fi * q0 * ψ(state.r0, state) + state.rs * sqrt(2 * m_fi * mp * h) for h in OSgrid.H]
    return Pmin, Pmax, Pstar, Psstar
end

#### get the orbit type
function classify_orbit(Pθ::Float64, H::Float64, ζ::Float64, state::PlasmaState, RSgrid::RealSpaceGrid)::Integer
    """
    -1 : forbidden
    0 : betatron
    1 : figure 8
    2 : type 1
    3 : cyclotron
    """
    Pmin = -m_fi * mp * state.r0 * sqrt(2 * H / (m_fi * mp)) + q_fi * q0 * ψ(state.r0, state)
    Pmax = m_fi * mp * rw * sqrt(2 * H / (m_fi * mp)) + q_fi * q0 * ψ(rw, state)
    Pstar = q_fi * q0 * ψ(state.r0, state)
    Psstar = q_fi * q0 * ψ(state.r0, state) + state.rs * sqrt(2 * m_fi * mp * H)

    if (Pθ < Pmin) || (Pθ > Pmax)
        -1
    else
        if Pθ < Pstar
            0
        elseif Pθ > Psstar
            3
        else
            #### now we need to compute the effective potential well gradient and look for roots
            roots = find_zeros(x -> ∇U_eff(Pθ, x, state) / E0, [RSgrid.r[1], rw])
            if length(roots) == 3
                if H*sin(ζ)^2 / E0 > U_eff(Pθ, roots[2], state) / E0
                    1 # over the central hump - figure 8 orbit
                else
                    3 # low energy cyclotron trapped in half of f8 well
                end
            else
                2 # type-1 case
            end
        end
    end
end

#-----------------
# Local ion values
#-----------------

#### create a monte carlo sample of the ion's phase space
function sample_ion(ion::FastIon, N::Int64, state::PlasmaState, RSgrid::RealSpaceGrid)::Tuple{Vector{Vector{Float64}},Vector{Vector{Float64}},Vector{Float64}}
    #### create a domain sample of the particle
    r_samp = ion.roots[1] .+ (ion.roots[2] - ion.roots[1]) * rand(rng, N) #sqrt.(ion.roots[1]^2 .+ (ion.roots[2]^2 - ion.roots[1]^2) .* rand(rng,N))
    θ_samp = 2 * π * rand(rng, N)
    z_samp = RSgrid.z[1] .+ (RSgrid.z[end] - RSgrid.z[1]) .* rand(rng, N)
    rdir = rand(rng, [-1,1], N)
    zdir = rand(rng, [-1,1], N)

    vθ = [get_vθ(ion.Pθ, ir, state) for ir in r_samp]
    U_samp = (0.5 * m_fi * mp) .* vθ.^2 .+ q_fi * q0 * [ϕ(ir, state) for ir in r_samp]
    vr = sqrt(2/(m_fi*mp)) * (rdir .* sqrt.(ion.H⊥ .- U_samp))

    return [[r_samp[i] * cos(θ_samp[i]), r_samp[i] * sin(θ_samp[i]), z_samp[i]] for i=eachindex(r_samp)], # xyz position
    [v_cyl_to_cart([vr[i], vθ[i], zdir[i] * ion.Vz], r_samp[i] * cos(θ_samp[i]), r_samp[i] * sin(θ_samp[i])) for i=eachindex(r_samp)], # v_xyz velocity
    2.0 ./ (ion.τr * abs.(vr)) # probability

end

#### get the radial probabilitiy
function get_πr(r::Float64, ion::FastIon, state::PlasmaState)::Float64
    if (r < ion.roots[1]) || (r > ion.roots[2])
        0.0
    else
        2 / (ion.τr * sqrt(2/(m_fi*mp)) * sqrt(ion.H⊥ - U_eff(ion.Pθ,r,state)))
    end
end

function get_vθ(Pθ::Float64, r::Float64, state::PlasmaState)::Float64
    if r <= 10^-3
        (Pθ .- q_fi * q0 * ψ(10^-3,state)) / (m_fi * mp * 10^-3)
    else
        (Pθ .- q_fi * q0 * ψ(r,state)) / (m_fi * mp * r)
    end
end

#### compute the effective potential
function U_eff(Pθ::Float64, r::Float64, state::PlasmaState)::Float64
    if r <= 10^-3
        (Pθ - q_fi * q0 * ψ(10^-3, state))^2 / (2 * m_fi * mp * (10^-3)^2) + q_fi * q0 * ϕ(10^-3, state)
    else
        (Pθ - q_fi * q0 * ψ(r, state))^2 / (2 * m_fi * mp * r^2) + q_fi * q0 * ϕ(r, state)
    end
end

#### gradient of the effective potential
function ∇U_eff(Pθ::Float64, r::Float64, state::PlasmaState)::Float64
    if r <= 10^-3
        -1 * (Pθ - q_fi * q0 * ψ(10^-3, state))^2 / (m_fi * mp * (10^-3)^3) + q_fi * q0 * ∇ϕ(10^-3, state)
    else
        -1 * (Pθ - q_fi * q0 * ψ(r, state)) * (Pθ - q_fi * q0 * ψ(r, state) + q_fi * q0 * r * ∇ψ(r, state)) / (m_fi * mp * r^3) + q_fi * q0 * ∇ϕ(r, state) 
    end
end

function Jθfast(r::Float64, ion::FastIon, state::PlasmaState, RSgrid::RealSpaceGrid)::Float64
    vol = (4 * π * r * RSgrid.Δr * (RSgrid.z[end] - RSgrid.z[1]))
    q0 * q_fi * get_πr(r, ion, state) * get_vθ(ion.Pθ, r, state) / vol
end

function get_guiding_center(x::Array{Float64}, v::Array{Float64})::Tuple{Vector{Float64},Float64}

    #### just gets the local guiding center: does not accurately predict ion motion
    v⊥ = sqrt(v[1]^2 + v[2]^2)
    ρL = m_fi * mp * v⊥ / (q_fi * q0 * norm(B(x, state))) # uses |B|
    xc = x[1] - ρL * abs(v[2]) / v⊥
    yc = ρL * (1 - abs(v[2]) / v⊥)
    return [xc, yc, x[3]], ρL
end

function get_action_variables(x::Array{Float64}, v::Array{Float64}, state::PlasmaState)
    r = sqrt(x[1]^2 + x[2]^2)
    vc = v_cart_to_cyl(v, x[1], x[2])
    Pθ = m_fi * mp * r * vc[2] + q_fi * q0 * ψ(r,state)
    H = 0.5 * m_fi * mp * norm(v)^2 + q_fi * q0 * ϕ(r, state)
    ζ = asin(sqrt(1 - 0.5 * m_fi * mp * v[3]^2 / H))
    return Pθ/P0, H*J_to_eV/1e3, rad2deg(ζ)
end

##==================================================================================================================
## Orbit Tracing - also 1D
##==================================================================================================================

#### method for fast calculation of radial action
function boris_action(x_0::Vector{Float64}, 
    v_0::Vector{Float64},
    state::PlasmaState,
    t_max::Float64)

    #### set up
    u_0 = RK04_push(x_0, v_0, -1, state)
    vr_0 = v_cart_to_cyl(v_0, x_0[1], x_0[2])[1]
    Nr = 0
    t = 0.0
    ract = 0.0

    #### run orbit trace
    u_1 = Float64[]
    v_1 = Float64[]
    x_1 = Float64[]
    vr_1 = 0.0
    while true

        if t > 0
            u_0 = u_1
            x_0 = x_1
            v_0 = v_1
            vr_0 = vr_1
        end

        u_1 = boris_push(x_0, u_0, state)
        x_1 = x_0 .+ u_1 * Δt
        v_1 = (u_1 .+ u_0) / 2.0
        vr_1 = v_cart_to_cyl(v_1, x_1[1], x_1[2])[1]
        t += Δt

        if sign(vr_1) != sign(vr_0)
            Nr += 1
        end

        if Nr == 1
            ract += m_fi * mp * vr_1^2 * Δt / (2 * π)
        end

        if !ion_in_domain(x_1) || (Nr > 1) || (t >= t_max)
            break
        end
    end
    return ract
end

function boris_push(x::Vector{Float64}, v::Vector{Float64}, state::PlasmaState)::Vector{Float64}
    vm = v + γ * E(x, state)
    h = γ * B(x, state)
    s = 2 * h /(1 + norm(h)^2)
    return vm + cross(vm + cross(vm, h), s) + γ * E(x, state)
end

function RK04_push(x0::Vector{Float64}, u0::Vector{Float64}, dir::Integer, state::PlasmaState)::Vector{Float64}
    dt = dir * Δt / 2
    k1 = _lorentz_accel(x0, u0, state)
    k2 = _lorentz_accel(x0 + u0 * dt + k1 * dt^2,
                       u0 + k1 * dt, state)
    k3 = _lorentz_accel(x0 + u0 * dt + k2 * dt^2,
                       u0 + k2 * dt, state)
    k4 = _lorentz_accel(x0 + u0 * dt + k3 * dt^2,
                       u0 + k3 * dt, state)
    return u0 + (Δt/6) * (k1 + 2*k2 + 2*k3 + k4)
end

function _lorentz_accel(x::Vector{Float64}, v::Vector{Float64}, state::PlasmaState)::Vector{Float64}
    q_fi*q0 / (m_fi*mp) * (E(x, state) + cross(v, B(x, state)))
end

function ion_in_domain(x::Vector{Float64})::Bool
    if sqrt(x[1]^2 + x[2]^2) >= rw
        false
    else
        true
    end
end