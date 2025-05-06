##==================================================================================================================
## Ionizing cross sections from Janev (1993)
##==================================================================================================================

function σEZ_Janev(E::Float64)::Float64
    """
    Electron-impact ionization
    Input : Electron energy (keV / amu)
    Output: Cross section (m^2)
    """
    # E represnts energy of electron in collision in units of eV
    E *= 1e3

    # Define fit parameters
    I = 13.6
    A  = +1.8450E-1
    B = zeros(5)
    B[1] = -3.2226E-2
    B[2] = -3.4539E-2
    B[3] = +1.4003E+0
    B[4] = -2.8115E+0
    B[5] = +2.2986E+0

    if E < I
        return 0.0
    end

    # Fit function
    prefac = 10^(-4) * 10^(-13) /(I * E)
    term1 = A * log(E/I)
    term2 = sum([B[j]*(1.0 - I/E)^j for j in 1:5])

    # Return result
    return prefac * (term1 + term2)
end

function σIZ_Janev(E::Float64)::Float64
    """
    Ion-impact ionization
    Input : Ion energy (keV / amu)
    Output: Cross section (m^2)
    """
    # E represnts energy of ion in collision in units of keV/amu

    # Define fit parameters
    A1 = +1.2899E+1
    A2 = +6.1897E+1
    A3 = +9.2731E+3
    A4 = +4.9749E-4
    A5 = +3.9890E-2
    A6 = -1.5900E+0
    A7 = +3.1834E+0
    A8 = -3.7154E+0

    # Fit function
    prefac = 10^(-4) * 10^(-16) * A1
    term1 = exp(-A2/E) * log(1.0 + A3 * E) / E
    term2 = A4 * exp(-A5*E) / (E^A6 + A7 * E^A8)

    # Return result
    return prefac * (term1 + term2)
end

function σCX_Janev(E::Float64)::Float64
    """
    Charge-exchange
    Input : Ion energy (keV / amu)
    Output: Cross section (m^2)
    """
    # E represnts energy of ion in collision in units of keV/amu

    # Define fit parameters
    A1 = +3.2345E+00
    A2 = +2.3588E+02
    A3 = +3.8371E-02
    A4 = +3.8068E-06
    A5 = +1.1832E-10
    A6 = +2.3713E+00

    # Fit function
    prefac = 10^(-4)
    term1 = 10^(-16) * A1 * log(A2/E + A6)
    term2 = 1.0 + A3*E + A4*E^3.5 + A5 * E^5.4

    # Return result
    return prefac * (term1 / term2)
end

function σ0_Janev(E::Float64)::Float64
    """
    Bundles charge exchange and ionization cross sections from Janev
    Input: E [keV / amu]
    Output: Cross section [m^-2]
    """
    σCX_Janev(E) + σIZ_Janev(E) + σEZ_Janev(E) 
end

function shuji_He_strip(K::Float64, m::Float64)::Float64
    function _lowE(_E)
        3.2067e-10*_E^7 - 4.7796e-8*_E^6 + 2.9349e-6*_E^5 - 9.5491e-5*_E^4 + 0.0017619*_E^3 - 0.018335*_E^2 + 0.10174*_E-0.074916
    end
    function _highE(_E)
        5.6482e-10*_E^4 - 2.2556e-7*_E^3 + 2.2012e-5*_E^2 - 0.00079101*_E + 0.19404
    end
    E = K / m # keV / amu
    if (E <= 35) && (E > 1.0)
        return _lowE(E)
    elseif (E < 1.0)
        return _lowE(1.0)
    else
        return _highE(E) 
    end
end

##==================================================================================================================
## Fusion cross sections from Bosch-Hale 
##==================================================================================================================

const BoschHaleParameters = Dict("DDn"=>Dict("A"=>(5.3701E+04, 3.3027E+02, -1.2706E-01, 2.9327E-05, -2.5151E-09),
                                        "B"=>(0.0, 0.0, 0.0, 0.0),
                                        "BG"=>31.3970),
                                "DDp"=>Dict("A"=>(5.5576E+04, 2.1054E+02, -3.2638E-02, 1.4987E-06, 1.8181E-10),
                                        "B"=>(0.0, 0.0, 0.0, 0.0),
                                        "BG"=>31.3970),
                                "DT"=>Dict("A"=>(6.9270E+04, 7.4540E+08, 2.0500E+06, 5.2002E+04, 0.0),
                                    "B"=>(6.3800E+01, -9.9500E-01, 6.9810E-05, 1.7280E-04),
                                    "BG"=>34.3827),
                                "DHe3"=>Dict("A"=>(5.7501E+06, 2.5226E+03, 4.5566E+01, 0.0, 0.0),
                                        "B"=>(-3.1995E-03, -8.5530E-06, 5.9014E-08, 0.0),
                                        "BG"=>68.7508))

function σFUSION_BoschHale(E::Float64, branch::String)::Float64
    # Bosch-Hale S function
    p = BoschHaleParameters[branch]
    numer = p["A"][1] + p["A"][2] * E + p["A"][3] * E^2 + p["A"][4] * E^3 + p["A"][5] * E^4
    denom = 1.0 + p["B"][1] * E + p["B"][2] * E^2 + p["B"][3] * E^3 + p["B"][4] * E^4
    S = numer / denom
    return S / (E * exp(p["BG"] / sqrt(E))) * 1e-28 # m^2
end

##==================================================================================================================
## Scattering processes
##==================================================================================================================

function rutherford(K, ν, ρ)
    0.25 * (q0^2 / (2 * K))^2 / ν^4 * (1 + 2*ρ*ν + ρ^2)^1.5 / (1 + ρ * ν)
end

function scatter_nn(vi::Array{Float64}, m1::Float64, m2::Float64)::Vector{Float64}

    ρ = m1 / m2

    #### assumes uniform σ(θ) in [0,2 degrees]
    θ = rand(rng) * deg2rad(5)
    ν = 2 * rand(rng) - 1
    θ = acos(ν) # scattering angle, COM frame
    ϑ = acos((ν + ρ) / sqrt(1 + 2*ρ*ν + ρ^2)) # scattering angle, lab frame
    γ = rand(rng) * 2 * π # free angle, either frame (no shift)

    #### scattering cross section (Rutherford)
    K0 = 0.5 * m1 * mp * norm(vi)^2
    σϑ = rutherford(K0, ν, ρ)
    
    #### scattered velocity normal - assumes vi_⊥ > 0 
    n = vi / norm(vi)
    ay = sqrt(n[1]^2 / (n[1]^2 + n[2]^2))
    xh = [sqrt(1 - ay^2), ay, 0]
    yh = cross(n, xh)
    nf = sin(ϑ)*cos(γ) * xh .+ sin(ϑ)*sin(γ) * yh .+ cos(ϑ) * n

    #### elastic collision for scattered energy
    Kf = K0 * (1 + 2*ρ*ν+ρ^2) / (1 + ρ)^2
    vf = sqrt(2 * Kf / (m1 * mp)) * nf

    #### conservation of momentum
    vs = ρ * (vi .- vf)

    return vf
end