##==================================================================================================================
## I/O Functions
##==================================================================================================================

##------------------------
## Loading input files
##------------------------

function load_plasma_state(run_params)
    #### initialize the plasma state
    if run_params["files"]["eqm_type"] == 0
        return _load_rigid_rotor(run_params)
    elseif run_params["files"]["eqm_type"] == 1
        return _load_solovev(run_params)
    elseif run_params["files"]["eqm_type"] == 2
        return _load_lreq(run_params)
    elseif run_params["files"]["eqm_type"] == 3
        return _load_aotc(run_params)
    else
        throw(ErrorException("Enter a valid equilibrium type (0-3)."))
    end
end

function load_orbit_space(run_params)::OrbitSpaceGrid
    bnd = Vector{Float64}[]
    ΔJ = Float64[]
    gridnames = ["ptheta_bounds", "H_bounds", "zeta_bounds"]
    Δgridnames = ["d_ptheta", "d_H", "d_zeta"]
    for i=1:3
        if run_params["orbit_space_grid"][gridnames[i]][1] isa Float64
            push!(bnd, run_params["orbit_space_grid"][gridnames[i]])
        else
            push!(bnd, [Float64(rp) for rp in run_params["orbit_space_grid"][gridnames[i]]])
        end

        if run_params["orbit_space_grid"][Δgridnames[i]] isa Float64
            push!(ΔJ, run_params["orbit_space_grid"][Δgridnames[i]])
        else
            push!(ΔJ, Float64(run_params["orbit_space_grid"][Δgridnames[i]]))
        end
    end
    NJ =  [Int64(round((bnd[i][2] - bnd[i][1]) / ΔJ[i])) + 1 for i=eachindex(ΔJ)]
    NJ[1] += 1 # fix for Pθ
    return OrbitSpaceGrid(bnd, ΔJ,
    ΔJ[1],
    ΔJ[2],
    ΔJ[3],
    ΔJ[1]*ΔJ[2]*ΔJ[3],
    LinRange(bnd[1][1], bnd[1][2], NJ[1]), # Pθ
    LinRange(bnd[2][1], bnd[2][2], NJ[2]), # H
    LinRange(bnd[3][1], bnd[3][2], NJ[3]), # ζ
    NJ[1]*NJ[2]*NJ[3]) # total number

end

function load_diagnostics(file)::Tuple{Vector{AbstractView}, Vector{AbstractSignal}, Vector{AbstractFunction}, Vector{AbstractFunction}}
    fid = h5open(file, "r")
    detectors = keys(fid)
    views = AbstractView[]
    signals = AbstractSignal[]
    geo_funcs = AbstractFunction[]
    η_funcs = AbstractFunction[]
    for group in detectors
        _set_up_diagnostic!(views, signals, geo_funcs, η_funcs, group, read(fid, group))
    end
    close(fid)
    return views, signals, geo_funcs, η_funcs
end

function load_beams(file)::Vector{NeutralBeam}
    fid = h5open(file, "r")
    beams = NeutralBeam[]
    for (i,name) in enumerate(keys(fid))
        _nb = read(fid[name])
        axes = get_beam_axes(_nb["source_grid_center_m"], _nb["source_grid_normal"])
        push!(beams, 
        NeutralBeam(name, i,
        _nb["source_grid_center_m"], 
        _nb["source_grid_normal"], 
        _nb["radius_source_grid_m"], 
        _nb["nominal_focal_length_m"], 
        deg2rad(_nb["nominal_divergence_semi_perp_deg"]), 
        deg2rad(_nb["nominal_divergence_semi_par_deg"]), 
        _nb["z_scrapers"], 
        _nb["r_scrapers"],
        _nb["voltage_kV"], 
        _nb["ion_current_A"], 
        [sqrt(2 * (_nb["voltage_kV"] / i) * 1e3 * eV_to_J / (mp * m_fn)) for i=1:3],
        _nb["component_fractions"], 
        _nb["neutralization_efficiency"], 
        axes[1],
        axes[2]))
    end
    close(fid)
    return beams
end

function load_fastn(file)
    throw(ErrorException("Fast neutral loading not built!"))
end

function load_cold_neutrals(run_params)
    if run_params["files"]["neut_type"] == 0
        nn = _load_constant_neutrals(run_params)
    elseif run_params["files"]["neut_type"] == 1
        nn = _load_lreq_neutrals(run_params)
    elseif run_params["files"]["neut_type"] == 2
        nn = _load_DEGAS2_neutrals(run_params)
    else
        throw(ErrorException("Neutral type must be [0, 1, 2]"))
    end
    return nn
end

function load_particle_data(run_params, OSgrid)
    
    build_source = false

    if run_params["files"]["part_type"] == 0
        build_source = true
        F = spzeros(OSgrid.N)
    elseif run_params["files"]["part_type"] == 1
        F = _load_MC_part_file(run_params["files"]["part_file"], OSgrid)
    elseif run_params["files"]["part_type"] == 2
        F = _load_ROTJ_part_file(run_params["files"]["part_file"], OSgrid)
    elseif run_params["files"]["part_type"] == 3
        F = _load_AOTC_part_file(run_params["files"]["part_file"], OSgrid)
    else
        throw(ErrorException("Enter a valid particle data type (0-3)."))
    end
    return F, build_source
end

##------------------------
## Saving results
##------------------------

function initialize_output_files(dir::String, file_id::String, OSgrid::OrbitSpaceGrid, signals::Vector{AbstractSignal})

    #### initialize weight file
    h5open("$(dir)/weights_" * "$(_get_wisdom_mode_string(wisdom_mode))_" * file_id * ".h5", "w") do file
    
        # Weight data
        weight_group = create_group(file, "weights")
        attributes(weight_group)["description"] = "(i,j,v) values for the W matrix (row-major)"

        # # jacobian determinant
        # jac_group = create_group(file, "jacobian")
        # attributes(jac_group)["description"] = "Jacobain determinant"
        # file["jacobian/values"] = 𝒥

        # Weight axes
        axgrp = create_group(file, "weight_axis")
        attributes(axgrp)["description"] = "Axes are indexed (ptheta, H, pitch)"
        file["weight_axis/values"] = get_full_axis(OSgrid)

        # get the grids, just in case
        create_group(file["weight_axis"], "grid")
        file["weight_axis/grid/ptheta"] = [i for i in OSgrid.Pθ]
        file["weight_axis/grid/H"] = [i for i in OSgrid.H]
        file["weight_axis/grid/pitch"] = [i for i in OSgrid.ζ]

        # signal axes
        siggrp = create_group(file, "signals")
        siggrp["names"] = [sig.name for sig in signals]

    end
end

function store_result!(W::SparseMatrixCSC, result::Tuple{Integer,Integer,Float64})
    result[3] > 0 ? W[result[1],result[2]] += result[3] : nothing
end

function save_fastn(dir::String, file_id::String, beams::Array{NeutralBeam}, fastn::SharedArray)
    h5open("$(dir)/fastn_" * "$(_get_wisdom_mode_string(wisdom_mode))_" * file_id * ".h5", "w") do file
        # Beam data
        for i=eachindex(beams)
            create_group(file, beams[i].name)
            file["$(beams[i].name)/x"] = fastn[i,1:3,:]
            file["$(beams[i].name)/v"] = fastn[i,4:6,:]
            file["$(beams[i].name)/weight"] = 55fastn[i,7,:]
        end
    end
end

function save_weights(dir::String, file_id::String, W::SparseMatrixCSC)
    I,J,V = findnz(W)
    h5open("$(dir)/weights_" * "$(_get_wisdom_mode_string(wisdom_mode))_" * file_id * ".h5", "r+") do file
        file["weights"]["i"] = I
        file["weights"]["j"] = J
        file["weights"]["v"] = V
    end
end

function save_signals(dir::String, file_id::String, S::Vector)
    h5open("$(dir)/weights_" * "$(_get_wisdom_mode_string(wisdom_mode))_" * file_id * ".h5", "r+") do file
        file["signals"]["values"] = S
    end
end

##==================================================================================================================
## Helper functions for loading different file types
##==================================================================================================================

##---------------------
## Equilibrium types
##---------------------

function _load_rigid_rotor(run_params)
    ## paramters
    rs = run_params["model_params"]["r_s"]
    B0 = run_params["model_params"]["B_0"]
    κ = RR_κ(rs)
    nm = run_params["model_params"]["n_max"]
    T0 = B0^2 / (4 * μ0 * nm) # RR temperature
    ψ0 = (rs^2 * B0) / (4 * κ)
    r0 = rs / sqrt(2.0)
    state = RigidRotor(rs, B0, κ, nm, T0, ψ0, r0)
    grid = _load_RSgrid(run_params)
    return state, grid
end

function _load_solovev(run_params)
    ## solovev
    ## really only useful for testing, doesn't prescribe a density profile
    rs = run_params["model_params"]["r_s"]
    B0 = run_params["model_params"]["B_0"]
    ψ0 = (rw^2 * B0 / 2) / (1 - 2 * (rw / rs)^2)
    r0 = rs / sqrt(2.0)
    κ = (1 - 2 * (rw / rs)^2) / 2
    ψ0 = rw^2 * B0 / (4 * κ)
    nm = run_params["model_params"]["n_max"]
    state = Solovev(rs, B0, κ, nm, ψ0, r0)
    grid = _load_RSgrid(run_params)
    return state, grid
end

function _load_lreq(run_params)

    ## LReq dumpfile format (2D)
    fid = h5open(run_params["files"]["eqm_file"], "r")
    ze = read(fid, "mesh/z")
    re = read(fid, "mesh/r")
    idz = argmin(abs.(ze))

    ## grids
    r = LinRange(re[1], re[end], length(re))
    Δr = abs(r[2] - r[1])
    Nz = Int(round(diff(run_params["misc"]["z_bounds"])[1] / Δr))
    z = LinRange(run_params["misc"]["z_bounds"][1], run_params["misc"]["z_bounds"][2], Nz)

    ## set up interpolations
    f_keys = ["vacuum/psi", "state/phi", "state/en", "state/te", "state/tm01"]
    itp = Dict()
    for key in f_keys
        _d = fid[key][:,idz]
        if key in ["vacuum/psi", "state/phi"]
            push!(itp, key=>cubic_spline_interpolation(r, _d, extrapolation_bc=Line()))
        else
            push!(itp, key=>interpolate((r,), _d, Gridded(Constant())))
        end

        if key == "vacuum/psi" # crude gradient calculation
            Δψ = zeros(length(r))
            Δψ[2:end-1] = [(_d[i+1] - _d[i-1]) / (2*Δr) for i=2:length(r)-1]
            Δψ[1] = (_d[2] - _d[1]) / Δr
            Δψ[end] = (_d[end] - _d[end-1]) / Δr
            push!(itp, "psi_gradient"=>cubic_spline_interpolation(r, Δψ))
        elseif key == "state/phi"
            Δϕ = zeros(length(r))
            Δϕ[2:end-1] = [(_d[i+1] - _d[i-1]) / (2*Δr) for i=2:length(r)-1]
            Δϕ[1] = (_d[2] - _d[1]) / Δr
            Δϕ[end] = (_d[end] - _d[end-1]) / Δr
            push!(itp, "phi_gradient"=>cubic_spline_interpolation(r, Δϕ))
        end
    end

    B0 = read(fid, "state/bbz")[argmin(abs.(re .- rw)),idz]
    close(fid)
    
    ## parameters
    if run_params["misc"]["is_frc"] == 1
        name = "LREQ_FRC"
        rs = find_zero(x -> itp["vacuum/psi"](x), 0.4)
        r0 = rs / sqrt(2)
        κ = (1 / 2) * (1 - 2 * (rw / rs)^2)
        ψ0 = rw^2 * B0 / (4 * κ) # not very useful
    else
        throw(ErrorException("is_frc = 0 not yet supported."))
    end

    state = LREQ(rs, B0, ψ0, r0)
    grid = _load_RSgrid(r, z)
    
    return state, grid, itp
end

function _load_aotc(run_params)
    throw(ErrorException("Not built yet!"))
end

##---------------------
## Background particle data
##---------------------

function _load_constant_neutrals(run_params)
    #### constant neutral density for testing
    nn(r::Float64)::Float64 = run_params["model_params"]["n_neut"]
    return nn
end

function _load_lreq_neutrals(run_params)
    #### dumpfile style neutral density
    fid = h5open(run_params["files"]["neut_file"], "r")
    rn = read(fid, "eqm/dens_n/X")
    zn = read(fid, "eqm/dens_n/Y")
    idr = argmin(abs.(rn .- rw))
    idz = argmin(abs.(zn))
    _r = LinRange(rn[1], rn[idr], length(rn[1:idr]))
    f_nn = read(fid, "eqm/dens_n/data")[1:idr,idz]
    close(fid)
    nn_itp = cubic_spline_interpolation(_r, f_nn)
    nn(r::Float64)::Float64 = nn_itp(r)
    return nn
end

function _load_DEGAS2_neutrals(run_params)
    fid = h5open(run_params["files"]["neut_file"], "r")
    rn = read(fid, "r")
    f_nn = read(fid, run_params["files"]["DEGAS2_HDF_PATH"])
    close(fid)
    _r = LinRange(rn[1], rn[end], length(rn))
    Σf_nn = [sum(f_nn[i,:]) for i=eachindex(rn)]
    nn_itp = cubic_spline_interpolation(_r, Σf_nn)
    nn(r::Float64)::Float64 = nn_itp(r)
    return nn
end

##---------------------
## Fast particle distributions
##---------------------

function _load_MC_part_file(file, OSgrid)

    #### load data
    fid = h5open(file, "r")
    z = read(fid, "part/ions/z")
    idz = abs.(z) .<= maximum(abs.(run_params["misc"]["z_bounds"]))
    Pθ = read(fid, "part/ions/ptheta")[idz] / P0 # normalized units
    H = read(fid, "part/ions/TE")[idz] * J_to_eV / 1e3 # keV
    Uz = read(fid, "part/ions/Uz")[idz]
    w = read(fid, "part/ions/weight")[idz] # particle weight
    close(fid)

    Kll = 0.5 * m_fi * mp * Uz.^2 * J_to_eV / 1e3 # keV
    ζ = rad2deg.(asin.(1 .- (Kll ./ H)))

    #### build F
    F = spzeros(OSgrid.N)
    for i=eachindex(Pθ)
        OSindex = get_index(Pθ[i], H[i], ζ[i], OSgrid)
        F[OSindex] += w[i]
    end

    return F
end

function _load_ROTJ_part_file(file, OSgrid)
    F = spzeros(OSgrid.N)
    return F
end 

function _load_AOTC_part_file(file, OSgrid)
    F = spzeros(OSgrid.N)
    return F
end 

##---------------------
## Grids!
##---------------------

function _load_RSgrid(run_params)
    Δx = run_params["misc"]["dr"]
    zbnd = [Float64(i) for i in run_params["misc"]["z_bounds"]]
    return RealSpaceGrid(LinRange(Δx, rw, Int((rw - Δx) / Δx) + 1), 
    LinRange(-rw, rw, Int(2*rw / Δx) + 1), 
    LinRange(-rw, rw, Int(2*rw / Δx) + 1), 
    LinRange(zbnd[1], zbnd[2], Int(diff(zbnd)[1] / Δx) + 1), 
    Δx, Δx)
end

function _load_RSgrid(r, z)
    Δr = abs(r[2] - r[1])
    Δz = abs(z[2] - z[1])
    Nr = Int(round(2*rw / Δr))+1
    return RealSpaceGrid(r, LinRange(-rw, rw, Nr), LinRange(-rw, rw, Nr),
    z, Δr, Δz)
end

##---------------------
## Diagnostics!
##---------------------

function _integrate_neutron_geo(D::AbstractAxis, R::AbstractAxis, f_geo::Array{Float64}, M::Integer, RSgrid::RealSpaceGrid, x0::Array{Float64}, n::Array{Float64}, h::Float64)::Array{Float64}
    xc = x0 .+ (h / 2) * n
    ΔZ = RSgrid.z[end] - RSgrid.z[1]
    V = 2 * π * ΔZ / M
    f_itp = linear_interpolation((D,R),f_geo)
    f_igt = zeros(length(RSgrid.r))
    for _=1:M
        r = rw * sqrt(rand(rng))
        θ = 2 * π * rand(rng)
        x = [r * cos(θ), r * sin(θ),  RSgrid.z[1] + ΔZ * rand(rng)]
        _d = (h / 2) + abs(dot((x .- xc), n))
        _r = norm((x .- x0) .- _d * n)
        f_igt[argmin(abs.(RSgrid.r .- r))] += f_itp(_d, _r) * V
    end
    return f_igt
end

function _set_up_diagnostic!(views::Array{AbstractView}, signals::Array{AbstractSignal}, 
    geo_funcs::Array{AbstractFunction}, η_funcs::Array{AbstractFunction}, name, obj)

    if name[1:3] == "NPA"
        
        #### collect geo function
        f_geo = obj["geo"]["f_geo"]
        R = obj["geo"]["grid"]["R"]
        D = obj["geo"]["grid"]["D"]
        push!(geo_funcs, linear_interpolation((D,R),f_geo))

        #### collect efficiency function
        f_η = permutedims(obj["efficiency"]["func"], reverse(1:ndims(obj["efficiency"]["func"]))) 
        E_η = obj["efficiency"]["energy"]
        push!(η_funcs, linear_interpolation(E_η, f_η))

        #### set up channels + views
        beam_numbers = [nb.name == "DNB-01" ? 0 : parse(Int, nb.name[end]) for nb in beams]
        Nb = length(beams)
        for i=1:Nb+1
            sig_index = length(signals)
            if i == 1
                nb_index = -1 # no beam
                view_name = name*"_passive"
            else
                nb_index = i-1
                nb_number = beam_numbers[nb_index]
                if nb_number in obj["active_beams"]
                    view_name = name*"_active_$(beams[nb_index].name)"
                else
                    continue
                end
            end
            for j=eachindex(obj["E_ch"])
                push!(signals, 
                NPAChannel(view_name*"_"*lpad(string(j),2,"0"), 
                view_name, 
                sig_index+j, 
                obj["E_ch"][j],
                obj["dE_ch"][j], 
                1.0, 
                obj["gain"],
                obj["E_ch"][j] - obj["dE_ch"][j] / 2,
                obj["E_ch"][j] + obj["dE_ch"][j] / 2))
            end
            #### derived view parameters
            h = norm(obj["ap_2"] .- obj["ap_1"]) # length of aperture cell
            n = (obj["ap_1"] .- obj["ap_2"]) ./ h # normal vector
            p1 = line_wall_intersect(obj["ap_2"], n)[1,:] # vessel pierce point (near-side)
            p2 = line_wall_intersect(obj["ap_2"], n)[2,:] # vessel pierce point (far side)
            xc = obj["ap_2"] .+ h * n
            xhat, yhat = get_⊥_axes(n)
            push!(views, 
            NPAChord(view_name, 
            obj["ap_1"], 
            obj["ap_2"], 
            obj["ap_radius"],
            nb_index,
            length(geo_funcs),
            length(η_funcs),
            signals[(sig_index+1):end],
            h, n, p1, p2, xc, xhat, yhat))
        end

    elseif name[1:4] == "FIDA"
        throw(ErrorException("FIDA not built."))

    elseif name[1:3] == "DND"

        #### collect geo function
        f_geo = obj["geo"]["f_geo"] #permutedims(obj["geo"]["f_geo"], reverse(1:ndims(obj["geo"]["f_geo"]))) 
        R = obj["geo"]["grid"]["R"]
        D = obj["geo"]["grid"]["D"]
        f_itg = _integrate_neutron_geo(D, R, f_geo, 10^3, RSgrid, obj["center"], obj["normal"], obj["length"])
        push!(geo_funcs, linear_interpolation(RSgrid.r, f_itg))

        #### push signal and view
        view_name = name
        sig_index = length(signals)
        push!(signals,
        NeutronChannel(view_name, sig_index+1))
        push!(views,
        NeutronDetector(view_name,
        obj["center"],
        obj["normal"],
        obj["length"],
        obj["radius"],
        length(geo_funcs),
        signals[sig_index+1:end]))

    elseif name[1:3] == "DBL"

        #### collect geo function
        f_geo = permutedims(obj["geo"]["f_geo"], reverse(1:ndims(obj["geo"]["f_geo"]))) 
        α = obj["geo"]["alpha"]
        β = obj["geo"]["beta"]
        push!(geo_funcs, linear_interpolation((α,β),f_geo))

        #### collect η
        f_η = permutedims(obj["efficiency"]["func"], reverse(1:ndims(obj["efficiency"]["func"]))) 
        E_η = obj["efficiency"]["energy"]
        push!(η_funcs, linear_interpolation(E_η, f_η))

        #### set up channels + views
        beam_numbers = [nb.name == "DNB-01" ? 0 : parse(Int, nb.name[end]) for nb in beams]
        Nb = length(beams)
        for i=1:Nb+1
            sig_index = length(signals)
            if i == 1
                nb_index = -1 # no beam
                view_name = name*"_passive"
            else
                nb_index = i-1
                nb_number = beam_numbers[nb_index]
                if nb_number in obj["active_beams"]
                    view_name = name*"_active_$(beams[nb_index].name)"
                else
                    continue
                end
            end
            push!(signals,
            PboloChannel(view_name, sig_index+1))
        
            #### derived view parameters
            n = obj["normal"]
            xhat, yhat = get_⊥_axes(n)
            push!(views, 
            PboloSensor(view_name,
            obj["center"],
            obj["normal"],
            obj["radius"],
            length(η_funcs),
            nb_index,
            length(geo_funcs),
            signals[sig_index+1:end],
            obj["alpha_bounds"],
            obj["f1_base"],
            obj["f1_vec"],
            obj["f2_base"],
            obj["f2_vec"],
            xhat, yhat))
        end
    end
end