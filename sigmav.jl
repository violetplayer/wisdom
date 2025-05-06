using HCubature
include("atomic.jl")

#### SI constants
const ε0 = 8.854e-12 # electric constant
const μ0 = 1.257e-6 # magnetic constant
const c = 1/sqrt(ε0*μ0) # speed of light
const q0 = 1.602e-19 # elementary charge
const mp = 1.673e-27 # proton mass
const me = 9.109e-31 # electron mass
const eV_to_J = 1.602e-19 # joule to eV conversion
const J_to_eV = 1 / eV_to_J # vice versa

function get_σv_beam_thermal(Eb::Float64, T::Float64, mb::Float64, mt::Float64)

    vb = sqrt(2 * Eb * 1e3 * eV_to_J / mb)
    vt = sqrt(2 * T * 1e3 * eV_to_J / mt)

    function integrand(v, θ)
        u = sqrt(vb^2 + v^2 - 2 * vb * v * cos(θ))
        σ = σCX_Janev(0.5 * mp * u^2 * J_to_eV * 1e3 / mp) # keV / amu
        return v^2 * exp(-v^2 / (2 * vt^2)) * sin(θ) * σ * u
    end

    a = (0.1, 0)
    b = (10 * vt, π)
    res = hcubature(integrand, a, b, rtol=1e-3)

    return res

end