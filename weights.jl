##==================================================================================================================
## Neutral Particle Analyzers
##==================================================================================================================

function geo_factor(x::Array{Float64}, view::NPAChord)::Float64
    D = dot((x .- view.ap_2), view.n)
    R = norm(x .- view.ap_2 .- D * view.n)
    geo_funcs[view.fgeo_index](D,R)
end

function η_factor(K::Float64, view::NPAChord)::Float64
    η_funcs[view.η_index](K) # expects units of keV
end

#### cold charge exchange weights
function nCX_weight(x::Array{Float64}, 
    v::Array{Float64}, 
    view::NPAChord, 
    state::PlasmaState)::Float64
    s = norm(v)
    K_o_m = 0.5 * mp * s^2 * J_to_eV * 1e-3 # keV / amu
    πCX = 1 - exp(-nn(sqrt(x[1]^2+x[2]^2)) * σCX_Janev(K_o_m) * Δt * s)
    πST = exp(-ST_cold_integral(x, view.prc_1, K_o_m, state))
    return πCX * πST 
end 

#### fast-fast charge exchange weights
function nCX_weight(x::Array{Float64}, 
    v::Array{Float64},
    view::NPAChord,
    nfast::Array{Float64},
    beam::NeutralBeam,
    state::PlasmaState)::Float64

    #### fast-beam charge exchange + fast-cold shinethrough
    πCX = 1 - exp(-sum([nfast[i] * σCX_fast_beam(v, beam, i) * norm(v .- beam.speeds[i] * beam.grid_normal) for i=1:3]) * Δt)
    πST = exp(-ST_cold_integral(x, view.prc_1, 0.5 * mp * norm(v)^2 * J_to_eV * 1e-3, state))
    return πST * πCX
end

#### throw particle backwards out of detector
function throw_particle(view::NPAChord, 
    state::PlasmaState, 
    OSgrid::OrbitSpaceGrid,
    RSgrid::RealSpaceGrid)::Tuple{Integer, Integer, Float64}

    good_particle = false
    ch = 0
    Ep = 0.0
    xp = Float64[]
    vp = Float64[]
    actvars = Float64[]
    while !good_particle

        #### get random numbers
        θ1 = rand(rng) * 2 * π
        θ2 = rand(rng) * 2 * π
        r1 = sqrt(rand(rng)) * view.r_a
        r2 = sqrt(rand(rng)) * view.r_a
        Ep = (OSgrid.H[1] + (OSgrid.H[end] - OSgrid.H[1]) * rand(rng)) # keV
        Δd = rand(rng)

        #### get channel
        ch = argmin([abs(Ep - signal.E_ch) for signal in view.signals])
        in_channel = (Ep >= view.signals[ch].E_dn) && (Ep <= view.signals[ch].E_up)

        #### get test particle position & velocity, chuck samples that hit in between channels
        if in_channel
            x1 = view.ap_1 + view.xhat * (r1 * cos(θ1)) + view.yhat * (r1 * sin(θ1))
            x2 = view.ap_2 + view.xhat * (r2 * cos(θ2)) + view.yhat * (r2 * sin(θ2))
            np = normalize!(x1 - x2)
            Li, Lo = get_wall_distances(x1, np)
            xp = x1 + (Li + (Lo - Li) * Δd) * np
            vp = -1 * sqrt(2 * Ep * 1e3 * eV_to_J / (m_fi * mp)) * np
            actvars = get_action_variables(xp, vp, state)
            good_particle = ion_exists(actvars[1] * P0, actvars[2] * sin(deg2rad(actvars[3]))^2 * eV_to_J * 1e3, state, RSgrid)
        end
    end

    #### get W-matrix indices
    OSindex = get_index(actvars[1], actvars[2], actvars[3], OSgrid)
    CHindex = view.signals[ch].index

    #### compute weight
    W = geo_factor(xp, view) * η_factor(Ep, view) * nCX_weight(xp, vp, view, state)

    return CHindex, OSindex, W
end

#### throw particle backwards out of detector [ACTIVE CX]
function throw_particle(view::NPAChord, 
    state::PlasmaState, 
    beams::Array{NeutralBeam}, 
    nf_itp::Any,
    OSgrid::OrbitSpaceGrid,
    RSgrid::RealSpaceGrid)::Tuple{Integer, Integer, Float64}

    good_particle = false
    ch = 0
    Ep = 0.0
    xp = Float64[]
    vp = Float64[]
    nb = Float64[]
    actvars = Float64[]
    while !good_particle

        #### get random numbers
        θ1 = rand(rng) * 2 * π
        θ2 = rand(rng) * 2 * π
        r1 = sqrt(rand(rng)) * view.r_a
        r2 = sqrt(rand(rng)) * view.r_a
        Ep =(OSgrid.H[1] + (OSgrid.H[end] - OSgrid.H[1]) * rand(rng)) # keV
        Δd = rand(rng)

        #### get channel
        ch = argmin([abs(Ep - signal.E_ch) for signal in view.signals])
        in_channel = (Ep >= view.signals[ch].E_dn) && (Ep <= view.signals[ch].E_up)

        if in_channel

            #### get test particle position & velocity
            x1 = view.ap_1 + view.xhat * (r1 * cos(θ1)) + view.yhat * (r1 * sin(θ1))
            x2 = view.ap_2 + view.xhat * (r2 * cos(θ2)) + view.yhat * (r2 * sin(θ2))
            np = normalize!(x1 - x2)
            Li, Lo = get_impact_bounds(x1, np,
            beams[view.beam_index].grid_center, beams[view.beam_index].grid_normal, 0.2)
            xp = x1 .+ (Li + (Lo - Li) * Δd) * np
            vp = -1 * sqrt(2 * Ep * 1e3 * eV_to_J / (m_fi * mp)) * np
            actvars = get_action_variables(xp, vp, state)
            nb = [nf_itp[view.beam_index][i](xp[1],xp[2],xp[3]) for i=1:3]
            good_particle = any(nb .> 0) && ion_exists(actvars[1] * P0, actvars[2] * sin(deg2rad(actvars[3]))^2 * eV_to_J * 1e3, state, RSgrid)
        end
    end

    #### get W-matrix indices
    
    OSindex = get_index(actvars[1], actvars[2], actvars[3], OSgrid)
    CHindex = view.signals[ch].index

    #### compute weight
    W = geo_factor(xp, view) * η_factor(Ep, view) * nCX_weight(xp, vp, view, nb, beams[view.beam_index], state)

    return CHindex, OSindex, W
end

##==================================================================================================================
## Bolometers
##==================================================================================================================

function geo_factor(x::Array{Float64}, v::Array{Float64}, view::PboloSensor)::Float64
    α, β = get_αβ_local(x, v)
    geo_funcs[view.fgeo_index](rad2deg(α), rad2deg(β)) #* (π * view.r^2 / (4 * π * norm(x .- view.x)^2))
end

function η_factor(view::PboloSensor, K::Float64)::Float64
    η_funcs[view.η_index](K) # expects units of keV
end

#### cold charge exchange weights
function nCX_weight(x::Array{Float64}, v::Array{Float64}, view::PboloSensor, state::PlasmaState)::Float64
    s = norm(v)
    K_o_m = 0.5 * mp * s^2 * J_to_eV * 1e-3 # keV / amu
    πCX = 1 - exp(-nn(sqrt(x[1]^2+x[2]^2)) * σCX_Janev(K_o_m) * Δt * s)
    πST = exp(-ST_cold_integral(x, view.x, K_o_m, state))
    return πCX * πST 
end

#### fast-fast charge exchange weights
function nCX_weight(x::Array{Float64}, 
    v::Array{Float64}, 
    view::PboloSensor, 
    state::PlasmaState, 
    beam::NeutralBeam,
    nb::Array{Float64})::Float64


    #### fast-beam charge exchange
    σv = [σCX_fast_beam(v, beam, i) * norm(v .- beam.speeds[i] * beam.grid_normal) for i=1:3]
    πCX = 1 - exp(-sum(nb .* σv) * Δt)

    #### fast-cold shinethrough
    πST = exp(-ST_cold_integral(x, view.x, 0.5 * mp * norm(v)^2 * J_to_eV * 1e-3, state))

    return πST * πCX
end

function throw_neutral(view::PboloSensor, 
    state::PlasmaState,
    OSgrid::OrbitSpaceGrid,
    RSgrid::RealSpaceGrid)::Tuple{Integer, Integer, Float64}

    good_particle = false
    N_tries = 0
    xp = Float64[]
    vp = Float64[]
    actvars = Float64[]
    fg = 0.0
    while !good_particle && N_tries < 100
        θ = rand(rng) * 2 * π
        r = sqrt(rand(rng)) * view.r
        θv = acos(2 * (0.5 + rand(rng) / 2) - 1) # only sample v_z > 0 in local frame
        ϕv = rand(rng) * 2 * π
        Ep = OSgrid.H[1] + (OSgrid.H[end] - (Uϕ_max * J_to_eV / 1e3) - OSgrid.H[1]) * rand(rng) # keV
        Δd = rand(rng)

        x0 = view.x .+ r*cos(θ) * view.xhat .+ r*sin(θ) *view.yhat
        np = sin(θv)*cos(ϕv) * view.xhat .+ sin(θv)*sin(ϕv) * view.yhat .+ cos(θv) * view.n
        vp = -1 * sqrt(2 * Ep * 1e3 * eV_to_J / (m_fi * mp)) * np
        _, Lo = get_wall_distances(x0, np)
        xp = x0 .+ Lo * Δd * np
        fg = geo_factor(x0, vp, view)
        actvars = get_action_variables(xp, vp, state)
        good_particle = (fg > 0) && (sqrt(xp[1]^2+xp[2]^2) < rw) && ion_exists(actvars[1] * P0, actvars[2] * sin(deg2rad(actvars[3]))^2 * eV_to_J * 1e3, state, RSgrid)
        N_tries += 1
    end

    OSindex = get_index(actvars[1], actvars[2], actvars[3], OSgrid)
    CHindex = view.signals[1].index

    if good_particle
        K_keV = actvars[2] - q_fi * q0 * ϕ(sqrt(xp[1]^2+xp[2]^2), state)*J_to_eV/1e3
        η = η_factor(view, K_keV)
        πTOT = nCX_weight(xp, vp, view, state)
        W = fg * η * πTOT
    else
        W = 0.0
    end

    return CHindex, OSindex, W
end

function throw_neutral(view::PboloSensor, 
    state::PlasmaState, 
    beams::Array{NeutralBeam}, 
    nf_itp::Any,
    OSgrid::OrbitSpaceGrid,
    RSgrid::RealSpaceGrid)::Tuple{Integer, Integer, Float64}

    good_particle = false
    Ntries = 0
    xp = Float64[]
    vp = Float64[]
    nb = Float64[]
    actvars = Float64[]
    fg = 0.0
    while !good_particle && Ntries < 100

        #### get random numbers
        θ = rand(rng) * 2 * π
        r = sqrt(rand(rng)) * view.r
        θv = acos(2 * (0.5 + rand(rng) / 2) - 1) # only sample v_z > 0 in local frame
        ϕv = rand(rng) * 2 * π
        Ep = OSgrid.H[1] + (OSgrid.H[end] - (Uϕ_max * J_to_eV / 1e3) - OSgrid.H[1]) * rand(rng) # keV

        #### construct particle
        x0 = view.x .+ r*cos(θ) * view.xhat .+ r*sin(θ) *view.yhat
        np = sin(θv)*cos(ϕv) * view.xhat .+ sin(θv)*sin(ϕv) * view.yhat .+ cos(θv) * view.n
        vp = -1 * sqrt(2 * Ep * 1e3 * eV_to_J / (m_fi * mp)) * np
        fg = geo_factor(x0, vp, view)
        if fg > 0
            Li, Lo = get_impact_bounds(x0, np / norm(np), 
            beams[view.beam_index].grid_center, beams[view.beam_index].grid_normal, 0.2)
            xp = x0 .+ (Li + (Lo - Li) * Δd) * np
            actvars = get_action_variables(xp, vp, state)
            nb = [nf_itp[view.beam_index][i](xp[1],xp[2],xp[3]) for i=1:3]
            good_particle = any(nb .> 0) && (sqrt(xp[1]^2+xp[2]^2) < rw) && ion_exists(actvars[1] * P0, actvars[2] * sin(deg2rad(actvars[3]))^2 * eV_to_J * 1e3, state, RSgrid)
        end
        Ntries += 1
    end

    #### get W-matrix indices
    OSindex = get_index(actvars[1], actvars[2], actvars[3], OSgrid)
    CHindex = view.signals[1].index

    #### compute weight
    if good_particle
        η = η_factor(view, actvars[2] - q_fi * q0 * ϕ(sqrt(xp[1]^2+xp[2]^2), state) * J_to_eV / 1e3)
        πTOT = nCX_weight(xp, vp, view, state, beams[view.beam_index], nb)
        W = η * πTOT * fg
    else
        W = 0.0
    end

    return CHindex, OSindex, W
end

function throw_ion(view::PboloSensor, 
    state::PlasmaState, 
    OSgrid::OrbitSpaceGrid,
    RSgrid::RealSpaceGrid)::Tuple{Integer, Integer, Float64}

    good_particle = false
    Ntries = 0
    xp = Float64[]
    vp = Float64[]
    actvars = Float64[]
    while !good_particle && Ntries < 100

        #### draw a random sample
        θ = rand(rng) * 2 * π
        r = sqrt(rand(rng)) * view.r
        θv = acos(2 * rand(rng) - 1)
        ϕv = rand(rng) * 2 * π
        Ep = (OSgrid.H[1] + (OSgrid.H[end] - OSgrid.H[1]) * rand(rng)) # keV
        
        #### position & velocity
        xp = view.x .+ r*cos(θ) * view.xhat .+ r*sin(θ) *view.yhat
        np = sin(θv)*cos(ϕv) * view.xhat .+ sin(θv)*sin(ϕv) * view.yhat .+ cos(θv) * view.n
        vp = -sqrt(2 * Ep * 1e3 * eV_to_J / (m_fi * mp)) * np
        actvars = get_action_variables(xp, vp, state)
        
        #### create circle to test fin intersection
        xc, ρL = get_guiding_center(xp, vp)
        sc1 = segment_circle_intersect(view.f1_base, view.f1_base .+ view.f1_vec, xc, ρL)
        sc2 = segment_circle_intersect(view.f2_base, view.f2_base .+ view.f2_vec, xc, ρL)

        good_particle = !sc1 & !sc2 & ion_exists(actvars[1] * P0, actvars[2] * sin(deg2rad(actvars[3]))^2 * eV_to_J * 1e3, state, RSgrid)
        Ntries += 1
    end

    OSindex = get_index(actvars[1], actvars[2], actvars[3], OSgrid)
    CHindex = view.signals[1].index

    if good_particle
        fg = geo_factor(xp, vp, view)
        η = η_factor(view, (H - q_fi * q0 * ϕ(rp, state)) * J_to_eV / 1e3)
        W = fg * η
    else
        W = 0.0
    end

    return OSindex, CHindex, W
end

function throw_particle(view::PboloSensor, 
    state::PlasmaState, 
    OSgrid::OrbitSpaceGrid,
    RSgrid::RealSpaceGrid)::Tuple{Integer, Integer, Float64}
    #if rand(rng, [true, false])
    return throw_neutral(view, state, OSgrid, RSgrid)
    #else
        #return throw_ion(view, state, OSgrid)
    #end
end

#### alias for conveninece in main.jl
function throw_particle(view::PboloSensor, 
    state::PlasmaState, 
    beams::Array{NeutralBeam}, 
    nf_itp::Any, 
    OSgrid::OrbitSpaceGrid,
    RSgrid::RealSpaceGrid)::Tuple{Integer,Integer,Float64}
    return throw_neutral(view, state, beams, nf_itp, OSgrid, RSgrid)
end

##==================================================================================================================
## Fusion Product Detectors
##==================================================================================================================

# function geo_factor(x::Array{Float64}, view::NeutronDetector)::Float64 
#     D = dot((x .- view.x), view.n)
#     R = norm(x .- view.x .- D * view.n)
#     geo_funcs[view.fgeo_index](D,R)
# end

function geo_factor(r::Float64, view::NeutronDetector)::Float64
    geo_funcs[view.fgeo_index](r)
    #1 / (4 * π^2 * norm([r,0,0] - (view.x + (view.h/2) * view.n))^2)
end

#### fusion weights
function nFUSION_weight(r::Float64, K::Float64, state::PlasmaState)::Float64
    return 1 - exp(-ni(r, state) * σFUSION_BoschHale(K, "DDn") * sqrt(2 * K * 1e3 * eV_to_J / (m_fi * mp)) * Δt)
end

# function throw_particle_alt(view::NeutronDetector,
#     state::PlasmaState,
#     RSgrid::RealSpaceGrid,
#     OSgrid::OrbitSpaceGrid)::Tuple{Integer,Integer,Float64}

#     good_particle = false
#     xp = Float64[]
#     vp = Float64[]
#     actvars = Float64[]
#     while !good_particle

#         #### get random numbers
#         θ = rand(rng) * 2 * π
#         r = sqrt(rand(rng)) * view.a
#         z = sqrt(rang(rng)) * view.h
#         θv = acos(1 - 2 * rand(rng)) # only sample v_z > 0 in local frame
#         ϕv = rand(rng) * 2 * π
#         Kn = 

#         x_n = view.x .+ (view.n * z) .+ (view.xhat * r * cos(θ)) .+ (view.yhat * r * sin(θ))
#         n_n = sin(θv)*cos(ϕv) * view.xhat .+ sin(θv)*sin(ϕv) * view.yhat .+ cos(θv) * view.n
#         Li, Lo = get_wall_distances(x_n, n_n)
#         if isfinite(Lo) && isfinite(Li) && (Li < Lo) # two ascending intersections required
#             xp = x_n .+ (Li + (Lo - Li) * rand(rng)) * n_n

#         end
#     end

#     #### get W-matrix indices
#     # OSindex = get_index(actvars[1], actvars[2], actvars[3], OSgrid)
#     # CHindex = view.signals[ch].index

#     # #### compute weight
#     # W = geo_factor(xp, view) * nCX_weight(xp, vp, view, state)

# end

function throw_particle(view::NeutronDetector,
    state::PlasmaState, 
    RSgrid::RealSpaceGrid, 
    OSgrid::OrbitSpaceGrid)::Tuple{Integer, Integer, Float64}

    #### get an orbit-space sample
    good_particle = false
    rts = Float64[]
    actvars = [0.0, 0.0, 0.0]
    H⊥ = 0.0
    while !good_particle
        actvars[1] = (OSgrid.Pθ[1] + (OSgrid.Pθ[end] - OSgrid.Pθ[1]) * rand(rng)) * P0
        actvars[2] = (OSgrid.H[1] + (OSgrid.H[end] - OSgrid.H[1]) * rand(rng)) * eV_to_J * 1e3
        actvars[3] = π/18 + (π / 2 - π / 18) * rand(rng)
        H⊥ = actvars[2] * sin(actvars[3])^2
        good_particle = ion_exists(actvars[1], H⊥, state, RSgrid) #&& orbit_is_confined(actvars[1], H⊥, state)
    end

    CHindex = view.signals[1].index
    OSindex = get_index(actvars[1]/P0, actvars[2]*J_to_eV/1e3, rad2deg(actvars[3]), OSgrid)

    # #### compute weight
    rts = get_turning_points(actvars[1], H⊥, state, RSgrid)
    τr = get_τr(actvars[1], H⊥, rts[1], rts[2], state)
    f0 = geo_factor(state.rs, view)
    π0 = nFUSION_weight(state.rs, E0 * J_to_eV / 1e3, state)
    function _integrand(r::Float64)::Float64
        πF = nFUSION_weight(r, (actvars[2] - q0 * q_fi * ϕ(r, state)) * J_to_eV / 1e3, state) / π0
        fg = geo_factor(r, view) / f0
        πR = 1.0 / sqrt(1 - U_eff(actvars[1], r, state) / H⊥)
        return r * πF * πR * fg
    end
    W0 = 2.0 / (τr * sqrt(2 * H⊥ / (m_fi * mp))) * f0 * π0 # prefactor
    W = W0 * quadgk(_integrand, rts[1], rts[2], rtol=0.01)[1]

    return CHindex, OSindex, W
end

##==================================================================================================================
## FIDA
##==================================================================================================================

## oh look, there's nothing here

##==================================================================================================================
## Physics signals
##==================================================================================================================

function throw_particle(view::FastCurrentDensity,
    state::PlasmaState, 
    RSgrid::RealSpaceGrid, 
    OSgrid::OrbitSpaceGrid)

    #### get an orbit-space sample
    good_particle = false
    while !good_particle
        global Pθ = (OSgrid.Pθ[1] + (OSgrid.Pθ[end] - OSgrid.Pθ[1]) * rand(rng)) * P0
        global H = (OSgrid.H[1] + (OSgrid.H[end] - OSgrid.H[1]) * rand(rng)) * eV_to_J * 1e3
        global ζ = deg2rad(OSgrid.ζ[1] + (OSgrid.ζ[end] - OSgrid.ζ[1]) * rand(rng))
        global vθ = get_vθ(Pθ, view.r, state)
        global U⊥ = 0.5 * m_fi * mp * vθ^2 + q_fi * q0 * ϕ(view.r, state)
        good_particle = ion_exists(Pθ, H, ζ, state, RSgrid) && orbit_is_confined(Pθ, H, ζ, state) && (U⊥ < H)
    end

    CHindex = view.signals[1].index
    OSindex = get_index(Pθ / P0, H * J_to_eV / 1e3, rad2deg(ζ), OSgrid)
    roots = get_turning_points(Pθ, H, ζ, state, RSgrid)
    H⊥ = H * sin(ζ)^2
    U⊥ = U_eff(Pθ, view.r, state)
    τr = get_τr(Pθ, H, ζ, roots[1], roots[2], state)
    vr = sqrt(2 * (H⊥ - U⊥) / (m_fi * mp))
    Jθ = (2 / (τr * vr)) * q0 * q_fi * vθ / (4 * π * view.r * view.Δr * view.Δz)

    return CHindex, OSindex, Jθ
end