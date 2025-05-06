##==================================================================================================================
## Setup file, loads in parameters and constants
##==================================================================================================================
using JSON
using HDF5
using Interpolations
using LinearAlgebra
using Roots
using Random
using QuadGK
using Distributions
using SharedArrays
using SparseArrays

#### bring in utility and I/O functions for initial setup
include("utils.jl") 
include("IO.jl")

#### load run parameters
const run_params = JSON.parsefile(input_deck)

#### status flags
const have_fastn = !isnothing(run_params["files"]["fastn_file"])
const wisdom_mode = run_params["misc"]["wisdom_mode"] # not presently used
const make_plots = Bool(run_params["misc"]["make_plots"])

#### SI constants
const ε0 = 8.854e-12 # electric constant
const μ0 = 1.257e-6 # magnetic constant
const c = 1/sqrt(ε0*μ0) # speed of light
const q0 = 1.602e-19 # elementary charge
const mp = 1.673e-27 # proton mass
const me = 9.109e-31 # electron mass
const eV_to_J = 1.602e-19 # joule to eV conversion
const J_to_eV = 1 / eV_to_J # vice versa

#### Initialize some run constants
const rw = run_params["misc"]["r_w"] # nominal wall radius
const q_th = Float64(run_params["misc"]["q_th"])   # thermal charge, elementary charge units
const q_fi = Float64(run_params["misc"]["q_fi"])   # fast ion charge
const m_th = Float64(run_params["misc"]["m_th"])   # thermal ion mass, atomic mass units
const m_fi = Float64(run_params["misc"]["m_fi"])   # fast ion mass
const m_fn = Float64(run_params["misc"]["m_fn"])   # fast neutral mass (only used for beam target)
const Nfastn = Int64(run_params["misc"]["N_beam"]) # number of samples to draw *per beam*
const Nion = Int64(run_params["misc"]["N_ion"])    # number of samples to draw *per ion*. Total number of ions determined by signal axis

#### set RNG
const rng = Xoshiro(run_params["misc"]["seed"])

#### plasma state
if run_params["files"]["eqm_type"] in [2, 3, 4]
    state, RSgrid, _itp = load_plasma_state(run_params)
    const itp = _itp # make interpolation functions constants, since they are referenced as globals
else
    state, RSgrid = load_plasma_state(run_params)
end

#### cold neutrals - function
nn = load_cold_neutrals(run_params)

#### derived constants
const Ω0 = (q_fi * q0) * state.B0 / mp
const τ0 = 2 * π / Ω0
const Δt = 0.01 * τ0
const v0 = rw / τ0
const P0 = mp * rw * v0
const E0 = P0 / τ0
const γ = q_fi * q0 * Δt / (2 * m_fi * mp)
const Uϕ_max = maximum([q_fi * q0 * ϕ(ir, state) for ir in RSgrid.r])

#### orbit-space axes
OSgrid = load_orbit_space(run_params)

#### dependent modules
include("atomic.jl")
include("collisions.jl")
include("ions.jl")
include("beams.jl")
include("weights.jl")

#### action normalization
# J0 = map_J0(OSgrid, RSgrid, state)
# const J0_itp = linear_interpolation((OSgrid.Pθ, OSgrid.H), J0)

#### jacobian determinant - needs fixing
# 𝒥 = get_jacobian(OSgrid)

#### beams and diagnostics
beams = load_beams(run_params["files"]["beam_file"])
views, signals, _geo_funcs, _η_funcs = load_diagnostics(run_params["files"]["diagnostic_file"])
const geo_funcs = _geo_funcs
const η_funcs = _η_funcs

#### some numbers
const Nstatus = Int(round(run_params["misc"]["status_frac"] * Nion * length(signals)))# how frequenctly to output information and save results
const N = Nion * length(views) # total number of weights to store

F, build_source = load_particle_data(run_params, OSgrid)

println("Setup completed.")