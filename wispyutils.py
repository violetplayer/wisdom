import numpy as np
from numba import njit

####==================================  NPA functions ======================================####

@njit
def npa_geo_integrand(x, y, xp):
    return xp[2] / ((xp[0] - x)**2 + (xp[1] - y)**2 + xp[2]**2)**(3/2)

@njit
def geo_npa_single(d, z, a, h, N=int(1e3)):
    #### create sample
    r_samp = a * np.sqrt(np.random.rand(N))
    t_samp = np.random.rand(N) * 2*np.pi
    x_samp = r_samp * np.cos(t_samp)
    y_samp = r_samp * np.sin(t_samp)
    #### projected circle from particle position x
    xc = d*h/(h-z)
    Rc = a*z/(z-h)
    #### only use sample points inside projected circle
    dr_c = np.sqrt((xc - x_samp)**2 + y_samp**2)
    idx = dr_c <= Rc
    N_c = np.sum(idx)
    if N_c > 0:
        I = np.sum(npa_geo_integrand(x_samp[idx], y_samp[idx], np.array([d, 0, z])))
        return a**2 * I / (4 * N_c)
    else:
        return 0.0

@njit
def geo_npa_all(d, z, a, h, N=int(1e4)):
    geo = np.zeros((len(d), len(z)))
    for i, _d in enumerate(d):
        for j, _z in enumerate(z):
            geo[i,j] = geo_npa_single(_d, _z, a, h, N)
    return geo

####================================== DNP05 ======================================####
dnp05 = {'D': 0.01,   # channel distance [m]
         'DD': 0.004, # channel width [m]
         'CH_OFFSET': {'c1': -5.797, 'c2': -6.432, 'c3': -7.410}, # [number of channels]
         'MAG_SEP_PARAMS': {'a': 0.152, 'I0': 0, 'V0': 33.4}}  # {[V/A], [A], [V]}
dnp05["CH_PARAMS"] = {'a': 0.00306, **dnp05["CH_OFFSET"]} # [keV^0.5 / A]

def ch_position(ch, c1, c2, c3):
    _ch = np.array(ch)
    idx1 = _ch <= 13
    idx2 = (13 < _ch) & (_ch <= 26)
    idx3 = (26 < _ch) & (_ch <= 39)
    x = np.zeros_like(_ch, dtype=float)
    x[idx1] = dnp05['D'] * (_ch[idx1] - c1)
    x[idx2] = dnp05['D'] * (_ch[idx2] - c2)
    x[idx3] = dnp05['D'] * (_ch[idx3] - c3)
    if x.ndim == 0:
        x = x[()]
    return x

def ch_energy(ch, I, moq):
    E = (ch_position(ch, **dnp05['CH_OFFSET']) * I * dnp05['CH_PARAMS']['a'] / dnp05['D'])**2 / moq * 1.0e3
    return E

def ch_delta_energy(ch, I, moq):
    DE = 2 / moq * (I * dnp05['CH_PARAMS']['a'] / dnp05['D'])**2 * ch_position(ch, **dnp05["CH_OFFSET"]) * dnp05['DD'] * 1.0e3
    return DE

####==================================  Neutron functions ======================================####


@njit
def get_ndet_angle_bounds(R, Z, a, h):
    phi1 = np.pi - np.arctan2(a, R)
    phi2 = np.pi + np.arctan2(a, R)
    if Z <= h / 2:
        th1 = np.pi / 2 - np.arctan2(h/2-Z, R-a)
        th2 = np.pi / 2 + np.arctan2(h/2+Z, R-a)
    else:
        th1 = np.pi / 2 + np.arctan2(Z-h/2, R+a)
        th2 = np.pi / 2 + np.arctan2(Z, R-a)
    return phi1, phi2, th1, th2

#%%===================
@njit
def get_ndet_single(R, Z, L_mfp, rc, hc, N=100):

    if Z < hc / 2:
        raise Exception("z > hc/2 only - assumes mirroring")

    if (Z < hc) & (R <= rc):
        return 0.0 # inside cylinder

    x = np.array([R, 0, Z])

    ### get angle sample intersecting detector cylinder
    phi1, phi2, th1, th2 = get_ndet_angle_bounds(R, Z, rc, hc)
    u_val = np.random.rand(N)
    phi_samp = phi1 + (phi2 - phi1) * u_val
    v_bnd = (1 - np.cos(np.array([th1,th2]))) / 2
    v_val = v_bnd[0] + np.diff(v_bnd) * np.random.rand(N)
    theta_samp = np.arccos(1 - 2 * v_val)

    #### build 3D velocity
    v = np.zeros((N,3))
    v[:,0] = np.sin(theta_samp) * np.cos(phi_samp)
    v[:,1] = np.sin(theta_samp) * np.sin(phi_samp)
    v[:,2] = np.cos(theta_samp)
    v /= np.sqrt(np.sum(v**2, axis=-1)).reshape((N,1))

    I = 0.0
    for i in range(N):
        L = get_cylinder_path_length(x, v[i,:], rc, hc)
        I += (1 - np.exp(-L / L_mfp))

    A_fac = (phi2 - phi1) * (np.cos(th1) - np.cos(th2)) / (4 * np.pi)
    return I * A_fac / (4 * np.pi**2 * N)

@njit
def get_cylinder_path_length(x, v, rc, hc):

    #### assumes cylinder has a base at the origin and is oriented in z-direction
    p = np.array([0, 0, 0])
    n = np.array([0, 0, 1])

    #### get line-cylinder intersections
    A1 = np.sum(np.cross(v, n) * np.cross(p - x, n))
    A2 = np.sum(np.cross(v, n) * np.cross(v, n))
    A3 = np.sum((p - x) * np.cross(v, n))**2
    S = A2 * rc**2 - A3

    if S <= 0: # doesn't hit the infinite cylinder
        return 0.0

    D1 = (A1 - np.sqrt(S)) / A2
    D2 = (A1 + np.sqrt(S)) / A2
    T1 = np.sum(n * (D1 * v - (p - x)))
    T2 = np.sum(n * (D2 * v - (p - x)))

    if T1 > hc: # hits plane first, maybe?
        P = np.sum(((p + hc * n) - x) * n) / np.sum(v * n)
        R = np.linalg.norm(x + P * v - (p + hc * n))
        if R < rc:
            path_length = np.linalg.norm((x + D1 * v) - (x + P * v))
        else:
            path_length = 0

    elif T2 < 0: # hits plane second, maybe?
        P = np.sum((p - x) * n) / np.sum(v * n)
        R = np.linalg.norm(x + P * v - p)
        if R < rc:
            path_length = np.linalg.norm((x + D2 * v) - (x + P * v))
        else:
            path_length = 0

    else: # hits finite cylinder both times
        path_length = np.linalg.norm((x + D1 * v) - (x + D2 * v))

    return path_length

@njit
def geo_ndet_all(R, Z, rc, hc, nH, sH, nC, sC, N=int(1e2)):
    L_mfp = 1 / (nH * sH + nC * sC)
    I = np.zeros((len(R),len(Z)))
    for i in range(len(R)):
        for j in range(len(Z)):
            I[i,j] = (nH * sH * L_mfp) * get_ndet_single(R[i],Z[j],L_mfp,rc,hc,N)
    return I

####==================================  Utility functions ======================================####

def v_cyl_to_cart(vr, vt, vz, x, y):
    r = np.sqrt(x**2 + y**2)
    vx = (x / r) * vr - (y / r) * vt
    vy = (y / r) * vr + (x / r) * vt
    return np.array([vx, vy, vz])