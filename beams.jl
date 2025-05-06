##==================================================================================================================
## Neutral Beams
## Most of the formalism in MC_fastn is stolen shamelessly from Scott Nicks' code, ROTJ: https://git.trialphaenergy.com/theory/rotj
##==================================================================================================================

#### accumulation for beam neutral target density
function throw_neutral(beam::NeutralBeam, state::PlasmaState, RSgrid::RealSpaceGrid)::Tuple{Integer,Integer,Integer,Vector{Float64}}
    """
    Function that creates a single weighted fast neutral.
    Inputs:
        beam   : NeutralBeam struct, contains beam geometry
        state  : any PlasmaState struct, used for shinethrough & CX probability calculation
        RSgrid : RealSpaceGrid struct
    Returns:
        i : x-axis index
        j : y-axis index
        k : z-axis index
        w : [3,]-vector containing weights for each beam component
    """

    #### Generate neutral that is not duct lost
    good_path = false
    while !good_path
        #### Grid birth position
        rg = sqrt(rand(rng)) * beam.grid_radius
        ϕg = rand(rng) * 2*pi
        rg_vec = rg*sin(ϕg) .* beam.ξ⊥ .+ rg*cos(ϕg) .* beam.ξll

        #### Get divergence angles + vector
        θll = randn(rng) * beam.σll
        θ⊥ = randn(rng) * beam.σ⊥
        d = θ⊥ .* beam.ξ⊥ .+ θll .* beam.ξll

        #### Beam neutral path and birthplace
        global p = normalize!(beam.grid_normal .- (rg_vec ./ beam.focal_length) .+ d)
        global s = rg_vec .+ beam.grid_center

        #### Check for duct losses!
        #### First get scraper projections...
        l_scrap = beam.z_scrapers ./ dot(p,beam.grid_normal)
        r_scrap = [norm(cross(s .+ l_scrap[i] * p .- beam.grid_center, beam.grid_normal)) for i=eachindex(beam.z_scrapers)]
            
        #### Check if neutral makes it through scrapers
        good_path = all(r_scrap .< beam.r_scrapers)
    end

    #### cold charge exchange, sampled speed from cold beam
    v_comp = sqrt.(2 * beam.voltage * [1.0, 1.0/2.0, 1.0/3.0] * 1e3 * eV_to_J / (m_fn * mp))
    E_keV_amu = 0.5 * mp * v_comp.^2 / (q0 * 1e3)  # keV / amu

    #### get random distance r < rw and assign a weight based on ST probability
    #### additionally, weight by neutral beam current + neutralization efficiency
    Li, Lo = get_wall_distances(s, p)
    l_samp = Li + (Lo - Li) * rand(rng)
    x = s .+ l_samp * p
    w_samp = [exp(-ST_cold_integral(s .+ Li * p, x, E_keV_amu[i], state)) for i=eachindex(E_keV_amu)]
    w_samp .*= beam.comp_fractions * beam.current * beam.neutralization_eff / q0

    #### get grid indidces
    i,j,k = get_xyz_grid_indices(x, RSgrid)

    return i,j,k,w_samp # set up for projection into SharedArray
end

#### accumulation for building source distribution
function throw_neutral(beam::NeutralBeam, state::PlasmaState, OSgrid::OrbitSpaceGrid)::Tuple{Integer,Float64}

    #### Generate ion that is confined + not duct lost
    good_ion = false
    #w_samp = 0.0
    E_keV_amu = 0.0
    xi = Float64[]
    xp = Float64[]
    actvars = Float64[]
    while !good_ion
        #### Grid birth position
        rg = sqrt(rand(rng)) * beam.grid_radius
        ϕg = rand(rng) * 2*pi
        rg_vec = rg*sin(ϕg) .* beam.ξ⊥ .+ rg*cos(ϕg) .* beam.ξll

        #### Get divergence angles + vector
        θll = randn(rng) * beam.σll
        θ⊥ = randn(rng) * beam.σ⊥
        d = θ⊥ .* beam.ξ⊥ .+ θll .* beam.ξll

        #### Beam neutral path and birthplace
        p = normalize!(beam.grid_normal .- (rg_vec ./ beam.focal_length) .+ d)
        s = rg_vec .+ beam.grid_center

        #### Check for duct losses!
        #### First get scraper projections...
        l_scrap = beam.z_scrapers ./ dot(p,beam.grid_normal)
        r_scrap = [norm(cross(s .+ l_scrap[i] * p .- beam.grid_center, beam.grid_normal)) for i=eachindex(beam.z_scrapers)]
            
        #### Check if neutral makes it through scrapers
        if all(r_scrap .< beam.r_scrapers)

            #### sampled speed from hot beam, T = 10 eV
            idv = rand(rng)
            if idv <= beam.comp_fractions[1]
                E_comp = rand(Normal(beam.voltage, 10e-3)) # keV
            elseif (idv <= beam.comp_fractions[1] + beam.comp_fractions[2]) && (idv > beam.comp_fractions[1])
                E_comp = rand(Normal(beam.voltage / 2, 10e-3)) # keV
            else
                E_comp = rand(Normal(beam.voltage / 3, 10e-3)) # keV
            end
            v_comp = sqrt(2 * E_comp * 1e3 * eV_to_J / (m_fn * mp))
            E_keV_amu = E_comp / (m_fn / mp) # keV / amu

            #### get distance
            Li, Lo = get_wall_distances(s, p)
            l_samp = Li + (Lo - Li) * rand(rng)
            xi = s .+ Li * p
            xp = s .+ l_samp * p

            # # give velocity a little inelastic scatter
            vp = scatter_nn(v_comp * p, m_fn, m_th)
            #vp = v_comp * normalize!(p .* (0.9 .+ 0.2 * rand(rng,3)))

            # is orbit confined?
            actvars = get_action_variables(xp, vp, state)
            good_ion = orbit_is_confined(actvars[1] * P0, (actvars[2] * sin(deg2rad(actvars[3]))^2) * eV_to_J * 1e3, state)
        end
    end

    #### weight by ST + current + neutralization_eff
    w_samp = exp(-ST_cold_integral(xi, xp, E_keV_amu, state))
    w_samp *= beam.current * beam.neutralization_eff / q0

    #### get OSgrid index
    i = get_index(actvars[1], actvars[2], actvars[3], OSgrid)

    return i, w_samp # set up for projection into SharedArray
end

#### accumulation for building source + scattered particles
function throw_neutral(beam::NeutralBeam, state::PlasmaState, OSgrid::OrbitSpaceGrid, scatter::Bool)::Tuple{Integer,Integer,Float64,Float64}

    #### Generate ion that is confined + not duct lost
    good_ion = false
    σϑ = 0.0
    E_keV_amu = 0.0
    xi = Float64[]
    xp = Float64[]
    actvars = Float64[]
    actvars_s = Float64[] # scattered particle
    while !good_ion
        #### Grid birth position
        rg = sqrt(rand(rng)) * beam.grid_radius
        ϕg = rand(rng) * 2*pi
        rg_vec = rg*sin(ϕg) .* beam.ξ⊥ .+ rg*cos(ϕg) .* beam.ξll

        #### Get divergence angles + vector
        θll = randn(rng) * beam.σll
        θ⊥ = randn(rng) * beam.σ⊥
        d = θ⊥ .* beam.ξ⊥ .+ θll .* beam.ξll

        #### Beam neutral path and birthplace
        p = normalize!(beam.grid_normal .- (rg_vec ./ beam.focal_length) .+ d)
        s = rg_vec .+ beam.grid_center

        #### Check for duct losses!
        #### First get scraper projections...
        l_scrap = beam.z_scrapers ./ dot(p,beam.grid_normal)
        r_scrap = [norm(cross(s .+ l_scrap[i] * p .- beam.grid_center, beam.grid_normal)) for i=eachindex(beam.z_scrapers)]
            
        #### Check if neutral makes it through scrapers
        if all(r_scrap .< beam.r_scrapers)

            #### sampled speed from hot beam, T = 10 eV
            idv = rand(rng)
            if idv <= beam.comp_fractions[1]
                E_comp = sample_beam(beam.voltage, 10) # keV
            elseif (idv <= beam.comp_fractions[1] + beam.comp_fractions[2]) && (idv > beam.comp_fractions[1])
                E_comp = sample_beam(beam.voltage/2, 10) #rand(Normal(beam.voltage / 2, 500e-3)) # keV
            else
                E_comp = sample_beam(beam.voltage/3, 10) #rand(Normal(beam.voltage / 3, 500e-3)) # keV
            end
            v_comp = sqrt(2 * E_comp * 1e3 * eV_to_J / (m_fn * mp))
            E_keV_amu = E_comp / m_fn # keV / amu

            #### get distance
            Li, Lo = get_wall_distances(s, p)
            l_samp = Li + (Lo - Li) * rand(rng)
            xi = s .+ Li * p
            xp = s .+ l_samp * p

            #### scattering check
            check = rand(rng)
            σIZ = σIZ_Janev(E_comp / m_fn)
            σCX = σCX_Janev(E_comp / m_fn)
            if check <= 0.0 #(σIZ / (σCX + σIZ))
                vp, vs, σϑ = scatter_nn(v_comp * p, m_fn, m_th)
            else
                vp = v_comp * p
                vs = v_comp * p # just duplicate it
                σϑ = 1.0
            end

            # is orbit confined?
            actvars = get_action_variables(xp, vp, state)
            actvars_s = get_action_variables(xp, vs, state)
            good_ion = true
            #good_ion = orbit_is_confined(actvars[1] * P0, (actvars[2] * sin(deg2rad(actvars[3]))^2) * eV_to_J * 1e3, state)
        end
    end

    #### weight by ST + current + neutralization_eff
    w_samp = exp(-ST_cold_integral(xi, xp, E_keV_amu, state))
    w_samp *= beam.current * beam.neutralization_eff / q0 / 2.0 # double sample
    
    #### get OSgrid index
    i = get_index(actvars[1], actvars[2], actvars[3], OSgrid)
    is = get_index(actvars_s[1], actvars_s[2], actvars_s[3], OSgrid)

    return i, is, w_samp, w_samp # set up for projection into SharedArray
end

function v_beam_weight(vn::Float64, vth::Float64, beam::NeutralBeam)
    1/sqrt(2*π*vth^2) * sum(beam.comp_fractions .* exp.(-(vn .- beam.speeds).^2 / (2 * vth^2)))
end

function interpolate_fast_neutrals(fastn::SharedArray{Float64}, 
    beams::Array{NeutralBeam}, 
    RSgrid::RealSpaceGrid)
    itp = []
    for i=eachindex(beams)
        push!(itp, [linear_interpolation((RSgrid.x, RSgrid.y, RSgrid.z), 
        fastn[i,j,:,:,:], 
        extrapolation_bc=0.0) for j=1:3])
    end
    return itp
end

function sample_beam(μ, σ)
    rand(Normal(sqrt(2 * μ * 1e3 * eV_to_J / (m_fn * mp)), 
    sqrt(2 * σ * eV_to_J / (m_fn * mp))))^2 * 0.5 * m_fn * mp * J_to_eV / 1e3 # keV
end