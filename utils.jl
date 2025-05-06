##==================================================================================================================
## Types and data structures
##==================================================================================================================

abstract type AbstractState end
abstract type AbstractView end
abstract type AbstractSignal end
abstract type AbstractParticle end
AbstractFunction = Union{Function, AbstractInterpolation} # generic function type
AbstractAxis = Union{Vector, AbstractRange} # generic axis type

#-----------------
# Particles
#-----------------

#### MEP fast ion
struct FastIon{X<:Float64, Xc<:ComplexF64, Y<:Integer, Z<:Bool, Xv<:Array{Float64}}<:AbstractParticle

    #### known ion properties
    Pθ::X # canonical angular momentum
    H::X # Total energy H = K + U
    ξ::Xc # sqrt(1 - Jr / J0); can be complex
    Ξ::X # ξ.re - ξ.im
    orbit_index::Y # index in J-space for W matrix

    #### constructed ion properties
    orbit_type::Y
    roots::Xv
    τr::X
    τθ::X
    Vz::X
    H⊥::X
    Jr::X
    J0::X
    confined::Z
    τse::X
    weight::X # typically 1, but is τr / ⟨τse⟩ for lost ions
end

#### Traced fast ion
struct FastIonTrace{Y<:Integer, Z<:Bool, Xv<:Vector{Float64}, Xvv<:Vector{Vector{Float64}}} <: AbstractParticle

    #### orbit parameters: known at start, although tracked
    Pθ::Xv # canonical angular momentum
    H::Xv  # Total energy H = K + U
    H⊥::Xv # Total perp. energy, H_⊥ = H - P_z^2/2m

    #### meta parameters
    solved::Z # has the orbit been solved?
    confined::Z # is the orbit confined?
    orbit_type::Y # what type of orbit is it?
    bounce_class::Y # what is the bounce classification

    #### solved orbit {x, v}
    x::Xvv # position
    v::Xvv # full-step velocity
    u::Xvv # half-step velocity
    x_cyl::Xvv # cylindrical position
    v_cyl::Xvv # full-step cylindrical velocity
    time::Xv

    #### derived parameters from orbit solve
    t_r::Xv         # radial bounce times
    t_z::Xv         # axial bounce times
    t_ract::Xv      # time points for radial action
    ract::Xv        # radial action
    t_zact::Xv      # time points for axial action
    zact::Xv        # axial action
end

struct FastNeutral{X<:Float64,Xv<:Array{Float64}}<:AbstractParticle
    weight::X
    x::Xv
    v::Xv
end

struct WarmNeutral{X<:Float64,Xv<:Array{Float64}}<:AbstractParticle
    weight::X
    x::Xv
    v::Xv
end

#-----------------
# States
#-----------------

abstract type PlasmaState <: AbstractState end
abstract type AbstractBeam <: AbstractState end

struct RigidRotor{X<:Float64}<:PlasmaState
    ## plasma equilibrium variables
    rs::X # separatrix radius
    B0::X # wall field, z = 0
    κ::X # field reversal variable. Different meanings depending on equilibrium
    nm::X # max density
    T0::X # nominal temperature 
    ψ0::X # minimum flux
    r0::X # null point radius
end

struct Solovev{X<:Float64}<:PlasmaState
    ## plasma equilibrium variables
    rs::X # separatrix radius
    B0::X # wall field, z = 0
    κ::X # field reversal variable. Different meanings depending on equilibrium
    nm::X # max density
    ψ0::X # minimum flux. Only applies to FRCs
    r0::X # null point radius
end

struct LREQ{X<:Float64}<:PlasmaState
    ## plasma equilibrium variables
    rs::X # separatrix radius
    B0::X # wall field, z = 0
    ψ0::X # minimum flux. Only applies to FRCs
    r0::X # null point radius
end

## beam geometry
struct NeutralBeam{X<:Float64, Y<:Integer, Z<:String, Xv<:Array{Float64}}<:AbstractBeam

    name::Z
    index::Y

    #### grid information
    grid_center::Xv
    grid_normal::Xv
    grid_radius::X
    focal_length::X
    σ⊥::X # perpendicular divergence
    σll::X # parallel divergence

    #### scrapers
    z_scrapers::Xv
    r_scrapers::Xv

    #### time-point information
    voltage::X
    current::X
    speeds::Xv # speed of 1,1/2,1/3 fractions

    #### other stuff
    comp_fractions::Xv # 1, 1/2, 1/3 component fractions
    neutralization_eff::X
    ξ⊥::Xv # semi-perpendicular axis
    ξll::Xv # semi-parallel axis
end

FRCState = Union{RigidRotor, Solovev, LREQ}

#-----------------
# Grids
#-----------------

struct RealSpaceGrid{X<:Float64, Y<:AbstractAxis}<:AbstractState
    r::Y
    x::Y
    y::Y
    z::Y
    Δr::X
    Δz::X
end

struct OrbitSpaceGrid{X<:Float64, Xv1<:Vector{Vector{Float64}}, Xv2<:Vector{Float64}, Xv3<:AbstractRange, Y<:Integer}<:AbstractState
    bnd::Xv1
    ΔJ::Xv2
    ΔPθ::X
    ΔH::X
    Δζ::X
    V::X
    Pθ::Xv3
    H::Xv3
    ζ::Xv3
    N::Y
end

#-----------------
# Diagnostic Channels
#-----------------

struct NPAChannel{X<:Float64, Y<:Integer, Z<:String}<:AbstractSignal
    name::Z
    view_name::Z
    index::Y
    E_ch::X
    ΔE_ch::X
    η_ch::X
    gain::X

    #### convenience
    E_dn::X
    E_up::X
end

struct PboloChannel{Y<:Integer, Z<:String}<:AbstractSignal
    name::Z
    index::Y
end

struct FIDAChannel{X<:Float64, Y<:Integer, Z<:String}<:AbstractSignal
    name::Z
    index::Y
end

struct NeutronChannel{Y<:Integer, Z<:String}<:AbstractSignal
    name::Z
    index::Y
end

struct ProtonChannel{X<:Float64, Y<:Integer, Z<:String}<:AbstractSignal
    name::Z
    index::Y
end

#-----------------
# Diagnostic Views
#-----------------

struct NPAChord{X<:Float64, Xv<:Array{Float64}, Xs<:Array{AbstractSignal}, Y<:Integer, Z<:String}<:AbstractView

    #### npa chord label
    name::Z

    #### view properties
    ap_1::Xv
    ap_2::Xv
    r_a::X

    #### indices
    beam_index::Y
    fgeo_index::Y
    η_index::Y

    #### associated signals
    signals::Xs

    #### derived parameters
    h::X # length of aperture cell
    n::Xv # normal vector
    prc_1::Xv # vessel pierce point (near-side)
    prc_2::Xv # vessel pierce point (far side)
    xc::Xv # aperture cone apex
    xhat::Xv # arbitrary x-axis of detector plane
    yhat::Xv # arbitrary y-axis of detector plane
end

struct PboloSensor{X<:Float64, Xv<:Array{Float64}, Xs<:Array{AbstractSignal}, Y<:Integer, Z<:String} <: AbstractView
    name::Z # name of Detector
    x::Xv # center point of sensor
    n::Xv # normal vector of sensor
    r::X # radius of sensor
    η_index::Y # efficiency function index η(E)
    beam_index::Y
    fgeo_index::Y
    signals::Xs

    αbnd::Xv # bounds of alpha [min, max]
    f1_base::Xv # base of first fin, centered
    f1_vec::Xv # vector of first fin, norm(f1_vec) = nominal length of fin
    f2_base::Xv
    f2_vec::Xv

    #### derived parameters
    xhat::Xv
    yhat::Xv
end

struct FIDAChord{X<:Float64, Xv<:Array{Float64}, Xs<:Array{AbstractSignal}, Y<:Integer, Z<:String} <: AbstractView
    name::Z # name of detector
    fgeo_index::Z
    signals::Xs

    #### derived parameters
    xhat::Xv
    yhat::Xv
end

#### assumes a cylindrical detector
struct NeutronDetector{X<:Float64, Xv<:Array{Float64}, Xs<:Array{AbstractSignal}, Y<:Integer, Z<:String} <: AbstractView
    name::Z
    x::Xv # base position of cylinder
    n::Xv # normal of cylinder
    h::X # length of cylinder
    a::X # radius of cylinder
    fgeo_index::Y 
    signals::Xs
end

struct ProtonDetector <: AbstractView
    name::String
    fgeo_index::Integer
    signals::Array{ProtonChannel}
end

CXView = Union{NPAChord, PboloSensor, FIDAChord} # views which can see charge exchange
FILDView = Union{PboloSensor} # views that see direct fast ion losses
FusionProductView = Union{NeutronDetector, ProtonDetector} # views that see fusion products

#-----------------
# "Physics" views
#-----------------

abstract type PhysicsView <: AbstractView end

struct FastCurrentDensity{X<:Float64} <: PhysicsView
    r::X # position
    Jθ::X # value [A/m^3]
end

##==================================================================================================================
## PlasmaState methods
##==================================================================================================================

#-----------------
# Rigid Rotor
#-----------------

function ψ(r::Float64, state::RigidRotor)::Float64
    (state.rs^2 * state.B0 / (4 * state.κ)) * log(cosh(state.κ * (2 * r^2 / state.rs^2 - 1)) / cosh(state.κ))
end

function ∇ψ(r::Float64, state::RigidRotor)::Float64
    state.B0 * r * tanh(state.κ * (2 * r^2 / state.rs^2 - 1))
end

function ϕ(r::Float64, state::RigidRotor)::Float64
    0.0
end

function ∇ϕ(r::Float64, state::RigidRotor)::Float64
    0.0
end

function ne(r::Float64, state::RigidRotor)::Float64
    state.nm / cosh(state.κ * (2 * r^2 / state.rs^2 - 1))^2
end

function ni(r::Float64, state::RigidRotor)::Float64
    state.nm / cosh(state.κ * (2 * r^2 / state.rs^2 - 1))^2
end

function Ti(r::Float64, state::Solovev)::Float64
    state.T0
end

function Te(r::Float64, state::Solovev)::Float64
    state.T0
end

#-----------------
# Solo'vev
#-----------------

function ψ(r::Float64, state::Solovev)::Float64
    ρ = r / state.rs
    state.ψ0 * ρ^2 * (1 - ρ^2)
end

function ∇ψ(r::Float64, state::Solovev)::Float64
    ρ = r / state.rs
    2 * state.ψ0 * ρ / state.rs * (1 - 2 * ρ^2)
end

function ϕ(r::Float64, state::Solovev)::Float64
    0.0
end

function ∇ϕ(r::Float64, state::Solovev)::Float64
    0.0
end

function ne(r::Float64, state::Solovev)::Float64
    state.nm
end

function ni(r::Float64, state::Solovev)::Float64
    state.nm
end

function Ti(r::Float64, state::Solovev)::Float64
    0.0
end

function Te(r::Float64, state::Solovev)::Float64
    0.0
end

#-----------------
# Lamy Ridge
#-----------------

function ψ(r::Float64, state::LREQ)::Float64
    itp["vacuum/psi"](r)
end

function ∇ψ(r::Float64, state::LREQ)::Float64
    itp["psi_gradient"](r)
end

function ϕ(r::Float64, state::LREQ)::Float64
    itp["state/phi"](r)
end

function ∇ϕ(r::Float64, state::LREQ)::Float64
    itp["phi_gradient"](r)
end

function ne(r::Float64, state::LREQ)::Float64
    itp["state/en"](r)
end

function ni(r::Float64, state::LREQ)::Float64
    itp["state/en"](r)
end

function Ti(r::Float64, state::LREQ)::Float64
    itp["state/tm01"](r)
end

function Te(r::Float64, state::LREQ)::Float64
    itp["state/te"](r)
end

##==================================================================================================================
## Field quantities
##==================================================================================================================

function B(x::Array{Float64}, state::PlasmaState)::Array{Float64}
    r = sqrt(x[1]^2 + x[2]^2)
    [0.0, 0.0, ∇ψ(r, state)/r]
end

function E(x::Array{Float64}, state::PlasmaState)::Array{Float64}
    [-∇ϕ(sqrt(x[1]^2 + x[2]^2),state), 0.0, 0.0]
end

##==================================================================================================================
## Geometry
##==================================================================================================================

function x_cart_to_cyl(xyz::Array{Float64})::Array{Float64}
    return [sqrt(xyz[1]^2 + xyz[2]^2), atan(xyz[2],xyz[1]), xyz[3]]
end

function x_cyl_to_cart(rθz::Array{Float64})::Array{Float64}
    return [rθz[1] * cos(rθz[2]), rθz[1] * sin(rθz[2]), rθz[3]]
end

function v_cart_to_cyl(v_xyz::Array{Float64}, x::Float64, y::Float64)::Array{Float64}
    r = sqrt(x^2 + y^2)
    vr = (x/r) * v_xyz[1] + (y/r) * v_xyz[2]
    vθ = (-y/r) * v_xyz[1] + (x/r) * v_xyz[2]
    return [vr, vθ, v_xyz[3]]
end

function v_cyl_to_cart(v_rθz::Array{Float64}, x::Float64, y::Float64)::Array{Float64}
    r = sqrt(x^2 + y^2)
    vx = (x/r) * v_rθz[1] - (y/r) * v_rθz[2]
    vy = (y/r) * v_rθz[1] + (x/r) * v_rθz[2]
    return [vx, vy, v_rθz[3]]
end

function line_circle_intersect(x::Array{Float64}, v::Array{Float64}, r_c::Float64)::Matrix{Float64}
    A = (v[1]^2 + v[2]^2)
    B = 2 * (x[1] * v[1] + x[2] * v[2])
    C = (x[1]^2 + x[2]^2) - r_c^2
    sqrt_term = B^2 - 4 * A * C
    res = NaN .* zeros(2,3)
    if sqrt_term > 1e-6
        res[1,:] = x .+ ((-B - sqrt(sqrt_term)) / (2 * A)) .* v
        res[2,:] = x .+ ((-B + sqrt(sqrt_term)) / (2 * A)) .* v
    elseif (sqrt_term < 1e-6) & (sqrt_term >= 0)
        res[1,:] .=  x .- (B / A) .* v
    end
    return res
end

function segment_circle_intersect(x1::Array{Float64}, x2::Array{Float64}, xc::Array{Float64}, R::Float64)::Bool

    A = [x1[1] - xc[1], x1[2] - xc[2]]
    B = [x2[1] - xc[1], x2[2] - xc[2]]
    d = norm(x1 - x2)

    α = (B[1] - A[1])^2 + (B[2] - A[2])^2
    β = 2 * (A[1] * (B[1] - A[1]) + A[2] * (B[2] - A[2]))
    γ = A[1]^2 + A[2]^2 - R^2
    Δ = β^2 - 4 * α * γ

    if Δ < 0
        return false
    elseif Δ == 0
        return -β/(2*α) <= d
    else
        return ((-β + sqrt(Δ))/(2*α) <= d) || ((-β - sqrt(Δ))/(2*α) <= d)
    end
end

function line_wall_intersect(x::Vector{Float64}, v::Vector{Float64})::Array{Float64}
    #### convenience function
    line_circle_intersect(x, v, rw)
end

function line_plane_intersect(x::Array{Float64}, v::Array{Float64}, p::Array{Float64}, n::Array{Float64})::Array{Float64}
    (x - p) - dot((x - p), n)/dot(v, n) * v
end

function point_line_distance(x::Array{Float64}, p::Array{Float64}, n::Array{Float64})::Float64
    norm(x .- p .- (dot(x .- p, n) * n))
end

function get_wall_distances(x::Vector{Float64}, v::Vector{Float64})::Tuple{Float64,Float64}
    res = line_wall_intersect(x, v)
    d1 = sqrt(sum((res[1,:] - x).^2))
    d2 = sqrt(sum((res[2,:] - x).^2))
    return d1, d2
end

function get_impact_bounds(x1::Vector{Float64}, v1::Vector{Float64}, x2::Vector{Float64}, v2::Vector{Float64}, Δ::Float64)
    n2 = cross(v2, cross(v1, v2))
    d = dot(x2 - x1, n2) / dot(v1, n2)
    return d - Δ, d + Δ
end

function rproj(x::Vector{Float64}, n::Vector{Float64}, t::Float64)
    sqrt(x[1]^2 + x[2]^2 + 2 * (x[1] * n[1] + x[2] * n[2]) * t + (n[1]^2 + n[2]^2) * t^2)
end

function get_αβ(vr::Float64, vθ::Float64, vz::Float64)::Tuple{Float64,Float64}
    return atan(vθ / vr), atan(vz / vr)
end

function get_αβ_local(x::Array{Float64}, v::Array{Float64})::Tuple{Float64,Float64}
    vcyl = v_cart_to_cyl(v, x[1], x[2])
    return get_αβ(vcyl[1], vcyl[2], vcyl[3])
end

function get_αβ_wall(x::Array{Float64}, v::Array{Float64})::Tuple{Float64,Float64}
    return get_αβ_local(line_wall_intersect(x, v)[1,:], v)
end

function get_⊥_axes(n::Array{Float64})::Tuple{Vector{Float64},Vector{Float64}}
    v = [0.0,0,0]
    v[1] = n[2] / sqrt(n[1]^2 + n[2]^2)
    v[2] = sqrt(1 - v[1]^2)
    
    xhat = v / norm(v)
    u = cross(n, xhat)
    yhat = u / norm(u)

    return xhat, yhat
end

function rodrigues_rotation(v::Array{Float64}, k::Array{Float64}, θ::Float64)::Array{Float64}
    v .* cos(θ) .+ cross(k, v) .* sin(θ) .+ k .* dot(k, v) .* (1 - cos(θ)) # mind the handedness of rotation!
end

##==================================================================================================================
## Array operations
##==================================================================================================================

function compose_index(indices::Tuple{Integer,Integer,Integer}, OSgrid::OrbitSpaceGrid)::Integer
    n2 = length(OSgrid.H)
    n3 = length(OSgrid.ζ)
    (indices[1]-1) * n2 * n3 + (indices[2]-1) * n3 + (indices[3]-1) + 1
end

function decompose_index(index::Integer, OSgrid::OrbitSpaceGrid)::Tuple{Integer, Integer, Integer}
    n2 = length(OSgrid.H)
    n3 = length(OSgrid.ζ)
    i = 1 + (index - 1) ÷ (n2 * n3)
    j = 1 + ((index - 1) - (i-1)*n2*n3) ÷ n3
    k = 1 + ((index - 1) - (i-1)*n2*n3 - (j-1)*n3)
    return i, j, k
end

function get_orbit_values(index::Integer, OSgrid::OrbitSpaceGrid)::Tuple{Float64, Float64, Float64}
    i,j,k = decompose_index(index, OSgrid)
    return OSgrid.Pθ[i], OSgrid.H[j], OSgrid.ζ[k]
end

function get_index(Pθ::Float64, H::Float64, ζ::Float64, OSgrid::OrbitSpaceGrid)::Integer
    i = argmin(abs.(OSgrid.Pθ .- Pθ))
    j = argmin(abs.(OSgrid.H .- H))
    k = argmin(abs.(OSgrid.ζ .- ζ))
    return compose_index((i,j,k),OSgrid)
end

function get_full_axis(OSgrid::OrbitSpaceGrid)::Array{Float64}
    full = zeros(Float64,(OSgrid.N,3))
    for i=1:OSgrid.N
        ii,ij,ik = decompose_index(i, OSgrid)
        full[i,1] = OSgrid.Pθ[ii]
        full[i,2] = OSgrid.H[ij]
        full[i,3] = OSgrid.ζ[ik]
    end
    return full
end

function OSvector_to_OSarray(v::Any, OSgrid::OrbitSpaceGrid)::Array{Float64}
    A = zeros(length(OSgrid.Pθ), length(OSgrid.H), length(OSgrid.ζ))
    for index=1:OSgrid.N
        i,j,k = decompose_index(index, OSgrid)
        A[i,j,k] = v[index]
    end
    return A
end

function OSarray_to_OSvector(A::Array{Float64}, OSgrid::OrbitSpaceGrid, sp::Bool)
    if sp
        v = spzeros(OSgrid.N)
    else
        v = zeros(OSgrid.N)
    end

    for index=1:OSgrid.N
        i,j,k = decompose_index(index, OSgrid)
        v[index] = A[i,j,k]
    end

    return v
end

function guess_roots(xm::Vector{Float64}, y::Vector{Float64})::Vector{Float64}
    s_ = sign.(y)
    s_[s_ .== 0] .= -1 # fix for where y = 0 exactly in the array
    idx = abs.(diff(s_)) .> 0
    if sum(idx) > 0
        xm[idx]
    else
        Float64[]
    end
end

function get_xyz_grid_indices(x::Array{Float64}, RSgrid::RealSpaceGrid)::Tuple{Integer,Integer,Integer}
    i = argmin(abs.(RSgrid.x .- x[1]))
    j = argmin(abs.(RSgrid.y .- x[2]))
    k = argmin(abs.(RSgrid.z .- x[3]))
    return i, j, k
end

function get_beam_axes(x0, n)
    d = -1 * sum(n .* x0)
    x_int = [0.0, 0.0, -d / n[end]]
    ξll = (x_int .- x0) / sqrt(sum((x_int .- x0).^2))
    return -1 * cross(n, ξll), ξll
end

function get_midpoints(x::AbstractAxis)::Vector{Float64}
    0.5 .* (x[2:end] + x[1:end-1])
end

function get_bin_edges(x::AbstractAxis)::Vector{Float64}
    # assumes uniform bin spacing
    y = zeros(length(x) + 1)
    Δx = x[2] - x[1]
    y[1:end-1] = x .- (Δx / 2)
    y[end] = x[end] + (Δx / 2)
    return y 
end

function create_1D_histogram(x::Vector{Float64}, bins::AbstractRange, w::Vector{Float64})::Vector{Float64}
    res = zeros(length(bins) - 1)
    for i=eachindex(res)
        idx = (x .> bins[i]) .& (x .<= bins[i+1])
        res[i] = sum(w[idx])
    end
    res /= sum(w)
    return res
end

##==================================================================================================================
## Misc mathematical operations
##==================================================================================================================

function RR_κ(rs::Float64)::Float64
    f_κ(x) = tanh(x) - x * (1 - 0.5 * (rs / rw)^2)
    find_zeros(f_κ, 0.1, 1.0)[1]
end

function normal(x, μ, σ)
    (1 / (sqrt(2 * π) * σ)) * exp.(-1 .* (x .- μ).^2 ./ (2 * π * σ^2))
end

function KL_metric(p::Array{Float64}, q::Array{Float64})::Float64
    id = (q .> 0) .&& (p .> 0)
    sum(p[id] .* log.(p[id] ./ q[id]) .+ q[id] .* log.(q[id] ./ p[id]))
end

##==================================================================================================================
## Other stuff
##==================================================================================================================

function _get_orbit_string(i)
    Dict(-1=>"Forbidden", 0=>"Betatron", 1=>"Figure 8",2=>"Type 1", 3=>"Cyclotron")[i]
end

function _get_wisdom_mode_string(i)
    Dict(0=>"wis1D_MEP", 1=>"wis1D_BP", 2=>"wis2D_MEP", 3=>"wis2D_BP")[i]
end