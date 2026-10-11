from cmath import tau
import numpy as np
from numba import njit
from time import time, strftime, gmtime
from scipy.constants import e, m_e, k as kB, hbar, mu_0
pi = np.pi
import matplotlib.pyplot as plt

## This first set of routine computes the energy dispersion and defines the Fermi Surface
###notes to njit:
## the njit function is used to accelerate the computation of the programm (pre-C compilation)
## it can slow down matplotlib, however do not remove if computing a complete T-sweep or B-sweep
## every function called by an njit-ed function needs to have njit prefix itself
@njit()

@njit()
def get_E(kx, ky, s_coeff):
    """Compute and return the in-plane dispersion using a TB representation."""

    a, b, c = s_coeff[:3] # lattice constants
    mu, t1, t2, t3, t4, t5 = s_coeff[3:9] # TB parameters

    E1 = t1/2*(np.cos(kx*a) + np.cos(ky*b)) # first nearest neighbour hopping energy
    E2 = t2*np.cos(kx*a)*np.cos(ky*b) # second n. n.
    E3 = t3/2*(np.cos(2*kx*a) + np.cos(2*ky*b)) # third n. n.
    E4 = t4/2*(np.cos(2*kx*a)*np.cos(ky*b) + np.cos(kx*a)*np.cos(2*ky*b)) # fourth n. n.
    E5 = t5*np.cos(2*kx*a)*np.cos(2*ky*b) # fifth n. n.

    return mu+E1+E2+E3+E4+E5



@njit()
def get_kyF(kxF, s_coeff):
    """At any given kx, compute the corresponding ky on the FS, if it exists."""

    a, b, c = s_coeff[:3] # lattice constants
    mu, t1, t2, t3, t4, t5 = s_coeff[3:9] # TB parameters

    ## Solve the equation E=0 for a given kx on the FS
    x = np.cos(kxF*a)

    u = t3 + t4*x + t5*(4*x**2-2)
    v = t1/2 + t2*x + t4*(x**2-1/2)
    w = mu + t1/2*x + t3*(x**2-1) - t4/2*x + t5*(1-2*x**2)

    y1 = -v/(2*u) + np.sqrt((v/(2*u))**2-w/u)
    y2 = -v/(2*u) - np.sqrt((v/(2*u))**2-w/u)

    ky1 = np.arccos(y1)/b
    ky2 = np.arccos(y2)/b

    if np.isnan(ky1): return ky2
    else: return ky1 # return np.nan if it does not exist

@njit()
def get_kF_lim_1(s_coeff):
    """Compute and return the extreme of the Fermi surface, that is the kx corresponding to ky=π/b."""

    a, b, c = s_coeff[:3] # lattice constants
    mu, t1, t2, t3, t4, t5 = s_coeff[3:9] # TB parameters

    ## Solve the equation E=0 for ky=π/b
    u = t3 - t4 + 2*t5
    v = t1/2 - t2 + t4/2
    w = mu - t1/2 + t4/2 - t5

    x1 = -v/(2*u) + np.sqrt((v/(2*u))**2-w/u)
    x2 = -v/(2*u) - np.sqrt((v/(2*u))**2-w/u)

    kx1 = np.arccos(x1)/a
    kx2 = np.arccos(x2)/a

    if (np.isnan(kx1)) and (np.isnan(kx2)):
        return (0., get_kyF(0., s_coeff))
    else:
        if np.isnan(kx1): return (kx2, pi/b)
        else: return (kx1, pi/b)

@njit()
def get_kF_lim_2(s_coeff):
    """Compute and return the extreme of the Fermi surface, that is the kx corresponding to ky=0."""

    a, b, c = s_coeff[:3] # lattice constants
    mu, t1, t2, t3, t4, t5 = s_coeff[3:9] # TB parameters

    ## Solve the equation E=0 for ky=0
    u = t3 + t4 + 2*t5
    v = t1/2 + t2 + t4/2
    w = mu + t1/2 - t4/2 - t5

    x1 = -v/(2*u) + np.sqrt((v/(2*u))**2-w/u)
    x2 = -v/(2*u) - np.sqrt((v/(2*u))**2-w/u)

    kx1 = np.arccos(x1)/a
    kx2 = np.arccos(x2)/a

    if (np.isnan(kx1)) and (np.isnan(kx2)):
        return (pi/a, get_kyF(pi/a, s_coeff))
    else:
        if np.isnan(kx1): return (kx2, 0.)
        else: return (kx1, 0.)

@njit()
def get_kF_lim(s_coeff):
    """Get the extremes of the Fermi arc in the first quadrant of the BZ."""
    kx_lim_1, ky_lim_1 = get_kF_lim_1(s_coeff)
    kx_lim_2, ky_lim_2 = get_kF_lim_2(s_coeff)

    return (kx_lim_1, ky_lim_1), (kx_lim_2, ky_lim_2)

@njit()
def get_kF(s_coeff, N_k, in_AFBZ):
    """Compute and return the Fermi wavevectors evenly spaced along the FS."""

    a, b, c = s_coeff[:3] # lattice parameters

    (kx_lim_1, ky_lim_1), (kx_lim_2, ky_lim_2) = get_kF_lim(s_coeff) # the extremes of the Fermi arc
    kxF = np.linspace(kx_lim_1, kx_lim_2, N_k)
    kyF = np.empty_like(kxF)
    for i in range(len(kxF)):
        kyF[i] = get_kyF(kxF[i], s_coeff)
    kyF[0], kyF[-1] = ky_lim_1, ky_lim_2

    ## Evenly space the k-points along the Fermi arc
    dkF = np.sqrt(np.diff(kxF)**2 + np.diff(kyF)**2) # segment lengths
    FS_l = np.zeros_like(kxF)
    FS_l[1:] = np.cumsum(dkF) # integrate FS path
    FS_grid = np.linspace(0, FS_l.max(), len(kxF)) # regular spaced path

    kxF = np.interp(FS_grid, FS_l, kxF)
    kyF = np.empty_like(kxF)
    for i in range(len(kxF)):
        kyF[i] = get_kyF(kxF[i], s_coeff)
    kyF[0], kyF[-1] = ky_lim_1, ky_lim_2

    if in_AFBZ: # map inside the AFBZ
        mask_in_AFBZ = np.abs(kyF)<=(pi/a-np.abs(kxF))
        kxF, kyF = kxF[mask_in_AFBZ], kyF[mask_in_AFBZ]

    ## Negative kxF are obtained by symmetry
    z_ind = 1 if 0. in kxF else 0 # avoid double zero point
    kxF = [-kx for kx in kxF[::-1]] + [kx for kx in kxF[z_ind:]]
    kyF = [ky for ky in kyF[::-1]] + [ky for ky in kyF[z_ind:]]

    return np.array(kxF), np.array(kyF)

@njit()
def get_v0(kx, ky, s_coeff, FLEX = False):
    """Compute and return the velocities computed as the gradient of the dispersion relation in k-space."""

    a, b, c = s_coeff[:3] # lattice constants
    mu, t1, t2, t3, t4, t5 = s_coeff[3:9] # TB parameters

    vx = a*np.sin(kx*a)*(-t1-2*t2*np.cos(ky*b)-t4*np.cos(2*ky*b)\
                                            -4*np.cos(kx*a)*(t3+t4*np.cos(ky*b)+2*t5*np.cos(2*ky*b)))/2/hbar

    vy = b*np.sin(ky*b)*(-t1-2*t2*np.cos(kx*a)-t4*np.cos(2*kx*a)\
                                                -4*np.cos(ky*b)*(t3+t4*np.cos(kx*a)+2*t5*np.cos(2*kx*a)))/2/hbar
    return -vx, -vy

@njit()
def get_DOS(kx, ky, s_coeff):
    """Compute and return the density of states."""

    vx, vy = get_v0(kx, ky, s_coeff)
    return 1/np.sqrt(vx**2+vy**2)/hbar


## Here, cyclotron motion around the FS is implemented
##
@njit()
def map_inside_BZ(kx, ky, s_coeff, in_AFBZ):
    """Return the wavevectors mapped inside the first BZ (or the AFBZ)."""

    a, b, c = s_coeff[:3] # lattice constants

    if not in_AFBZ:
        if ky > pi/b: ky -= 2*pi/b
        if ky < -pi/b: ky += 2*pi/b
        if kx > pi/a: kx -= 2*pi/a
        if kx < -pi/a: kx += 2*pi/a

    if in_AFBZ:
        if ky > (pi/a+kx): kx, ky = kx+pi/a, ky-pi/b
        if ky > (pi/a-kx): kx, ky = kx-pi/a, ky-pi/b
        if ky < (kx-pi/a): kx, ky = kx-pi/a, ky+pi/b
        if ky < (-kx-pi/a): kx, ky = kx+pi/a, ky+pi/b

    return kx, ky

@njit()
def evolve_kF(kx_0, ky_0, B, dt, N_t, s_coeff, in_AFBZ, FLEX):
    """Evolve the Fermi wavevector along the FS according to cyclotron motion
    using Runge-Kutta 4-th order algorithm."""

    a, b, c = s_coeff[:3] # lattice constants

    kx_t, ky_t = np.zeros(N_t), np.zeros(N_t)
    kx_t[0], ky_t[0] = kx_0, ky_0

    dk = e*B/hbar
    # print("Nt :", N_t)
    for t in range(1, N_t):
        vx, vy = get_v0(kx_t[t-1], ky_t[t-1], s_coeff, FLEX) #we use v0 and not v because the CVC doesn't influence the FS (not FSR)
        K1_x, K1_y = dk*vy, -dk*vx

        vx, vy = get_v0(kx_t[t-1]+dt*K1_x/2, ky_t[t-1]+dt*K1_y/2, s_coeff, FLEX)
        K2_x, K2_y = dk*vy, -dk*vx

        vx, vy = get_v0(kx_t[t-1]+dt*K2_x/2, ky_t[t-1]+dt*K2_y/2, s_coeff, FLEX)
        K3_x, K3_y = dk*vy, -dk*vx

        vx, vy = get_v0(kx_t[t-1]+dt*K3_x, ky_t[t-1]+dt*K3_y, s_coeff, FLEX)
        K4_x, K4_y = dk*vy, -dk*vx

        kx = kx_t[t-1] + dt/6*(K1_x + 2*K2_x + 2*K3_x + K4_x)
        ky = ky_t[t-1] + dt/6*(K1_y + 2*K2_y + 2*K3_y + K4_y)

        kx_t[t], ky_t[t] = map_inside_BZ(kx, ky, s_coeff, in_AFBZ) # at every step check if it is inside the BZ

    return kx_t, ky_t


## These routines define the scattering rate
#
#
def get_phi(kx, ky, s_coeff, hole_like):
    """Returns the azimuthal angle φ."""

    a, b, c = s_coeff[:3] # lattice constants

    if not hole_like: # centered at Γ
        return np.arctan2(ky, kx)

    else: # centered at (π,π)
        x = (pi/a-np.abs(kx))
        y = (pi/b-np.abs(ky))
        phi = (-pi/2-np.arctan(x/y)) * ((kx>=0)&(ky>=0))
        phi += -np.arctan(y/x) * ((kx<0)&(ky>=0))
        phi += np.arctan(y/x) * ((kx<0)&(ky<0))
        phi += (pi/2+np.arctan(x/y)) * ((kx>=0)&(ky<0))
        return phi


@njit()
def get_tau_inv(kx, ky, T, s_coeff, tau_choice, tau_coeff, hole_like):
    """Compute and return the scattering rate at a given azimuthal angle phi at temperature T."""

    phi = get_phi(kx, ky, s_coeff, hole_like) # compute phi array (angle on FS, rad) regarding kx, ky

    if tau_choice == 'gamma 0': # tau fully isotropic value : timp (impurity scattering)
        timp = tau_coeff[0]
        return timp + phi*0

    elif tau_choice == 'hot spots elastic': # scattering rate with a peak at the Hot Spots and no T dependence
        timp = tau_coeff[0]
        s_or_c = tau_coeff[-2]
        if int(s_or_c) == 0 :
            mode = 'charge'
        elif int(s_or_c) == 1:
            mode = 'spin'
        else: raise ValueError(f"{tau_coeff[9]} Not a valid value for `mode`; " +\
                         "supported values are `0` for charge and `1` for spin.")

        HS = tau_coeff[-1]
        elastic = timp + compute_tau_HS(mode,phi,HS)
        return elastic

    elif tau_choice == 'hot spots': # peak at hot spots + e-e scattering + anisotropic inelastic scattering peaking at HS

        timp = tau_coeff[0] #impurity scattering coefficient
        s_or_c = tau_coeff[-2] # mode spin or charge
        if int(s_or_c) == 0 :
            mode = 'charge'
        elif int(s_or_c) == 1:
            mode = 'spin'
        else:
            raise ValueError(f"{tau_coeff[-2]} Not a valid value for `mode`; " +\
                         "supported values are `0` for charge and `1` for spin.")

        alpha1 = tau_coeff[1] # T linear coefficient
        alpha2 = tau_coeff[2] # T squared coefficient
        HS = tau_coeff[-1] # location of the hot spots : always the last coeff

        elastic = timp*imp_x_f(s_coeff) + compute_tau_HS(mode,phi,HS) # elastic scattering rate
        inelastic = g(s_coeff)*alpha1*T + alpha2*T**2 # inelastic scattering rate
        e_e_scat = 0.0022e12*(1 + 0.4*np.sin(2*phi)**2)*T**2 # scattering rate coming from the e-e scattering
        return elastic + (compute_beta(mode, phi, HS)*inelastic) + e_e_scat



@njit
def get_HS(s_coeff):
    """compute the angle of hotspot in spin modelisation (i.e. intersection between FS and AFBZ) in [0, pi/4]"""
    a,b = s_coeff[:2]
    kxF, kyF = get_kF(s_coeff, 150, False)
    ky_AFBZ = pi/b - kxF*a/b
    i_HS = np.argmin(np.abs(np.sqrt(kyF**2+kxF**2)-np.sqrt(ky_AFBZ**2+kxF**2))) # find index of closest point on FS to AFBZ
    kxF_HS, kyF_HS = kxF[i_HS], kyF[i_HS]
    return np.arctan(kyF_HS/kxF_HS) #get the angle of the hotspot

@njit
def get_k_HS(s_coeff):
    """compute the angle of hotspot in spin modelisation (i.e. intersection between FS and AFBZ) in [0, pi/4]"""
    a,b = s_coeff[:2]
    kxF, kyF = get_kF(s_coeff, 150, False)
    ky_AFBZ = pi/b - kxF*a/b
    i_HS = np.argmin(np.abs(np.sqrt(kyF**2+kxF**2)-np.sqrt(ky_AFBZ**2+kxF**2))) # find index of closest point on FS to AFBZ
    kxF_HS, kyF_HS = kxF[i_HS], kyF[i_HS]
    return kxF_HS, kyF_HS #get the angle of the hotspot

@njit
def compute_tau_HS(mode,phi,HS) :
    sigma = 0.05 # width of the gaussian peak of the HS
    peak = 22e12 # height of the gaussian peak of the HS
    if mode == 'charge':
        # peak in the charge mode --> peak at the antinodes
        x = np.abs(np.remainder(phi,pi/2))
        for i in range(x.size):
            if x[i] >= pi/4:
                x[i] -= pi/2
        tau_HS = peak*np.exp(-x**2/(2*sigma**2))

    elif mode == 'spin':
        # peak in the spin mode --> peak at the nodes
        x = np.abs(np.remainder(phi,pi/2))
        for i in range(x.size):
            if x[i]> pi/4 :
                x[i] = pi/2 -  x[i]
        tau_HS = peak*np.exp(-(x-HS)**2/(2*sigma**2)) +peak*np.exp(-(x+HS)**2/(2*sigma**2))
    else:
        raise ValueError(f" {mode} Not a valid value for `mode`; " +\
                         "supported values are `charge` or `spin`")
    return tau_HS

@njit
def g(s_coeff):
    """doping dependancy of the inelastic scattering rate of LCCO [0.14, 0.15, 0.16, 0.18]"""
    mu = s_coeff[3]/e
    if mu == -0.0132 : #doping x = 0.14 :
        return 1.9

    elif mu == 0. : #doping x ==0.15
        return 1.

    elif mu == 0.0132 : #doping x ==0.16
        return 0.5

    elif mu == 0.0396 : #doping x ==0.18
        return 0.2
    else :
        return 1 # when other dopings or metals than LCCO [0.14, 0.15, 0.16, 0.18]

@njit
def imp_x_f(s_coeff):
    """doping dependancy of the impurity scattering rate for LCCO [0.14, 0.15, 0.16, 0.18]"""
    mu = s_coeff[3]/e
    if mu == -0.0132 : #doping x = 0.14 :
        return 1.1

    elif mu == 0. : #doping x ==0.15
        return 1.2

    elif mu == 0.0132 : #doping x ==0.16
        return 1.3

    elif mu == 0.0396 : #doping x ==0.18
        return 1.4 # when other dopings or metals than LCCO [0.14, 0.15, 0.16, 0.18]

@njit()
def sig(s_coeff):
    """doping dependancy of the width of the gaussian in the CVC correlation length for LCCO [0.14, 0.15, 0.16, 0.18]"""
    mu = s_coeff[3]/e
    # return 0.14
    if mu == -0.0132 : #doping x = 0.14 :
        return 0.16
    elif mu == 0. : #doping x ==0.15
        return 0.13
    elif mu == 0.0132 : #doping x ==0.16
        return 0.11
    elif mu == 0.0396 : #doping x ==0.18
        return 0.09
    else :
        return 0.13

@njit
def compute_beta(mode,phi,HS) :
    """compute the angle dependance of the inelastic term in the scattering rate"""
    if mode == 'charge' :
        # peak at the antinodes
        return np.cos(2*phi)**2
    elif mode == "spin" :
        # peak at the nodes
        return (np.sin(2*phi + 2*(pi/4-HS))**12 + np.sin(2*phi - 2*(pi/4-HS))**12)

    else:
        raise ValueError("Not a valid value for `mode`; " +\
                         "supported values are `charge` or `spin`")


############################
## FLEX-CVC approximation ##
############################

@njit
def get_epsilon(kx, ky, s_coeff, T=-1, mode = "spin") :
    """function to pondere the CVC (plays the role of the correlation length), computes every HS in 0,pi"""
    hole_like= False # needed to have continuous values in the case of a hole-like FS (we don't use it here : we rather shift the FS by pi)
    phi = get_phi(kx, ky, s_coeff, hole_like)
    x = (np.remainder(phi,pi/2))
    sigma = sig(s_coeff)
    k = 1
    if T != -1 : # If there's a temperature dependance of the CVC
        k = get_CVC_T_dep(T)
        sigma = sig(s_coeff)*get_CVC_T_dep(T)
    # computation in the spin mode
    HS = get_HS(s_coeff)
    eps= np.exp(-(x-HS)**2/(2*sigma**2))+ np.exp(-(x+HS)**2/(2*sigma**2))
    x  = pi/2-x
    eps += np.exp(-(x-HS)**2/(2*sigma**2))+ np.exp(-(x+HS)**2/(2*sigma**2))
    return eps*k

@njit()
def get_epsilon_not_modulo(kx, ky, s_coeff, HS, T=-1, mode = "spin") :
    """Fonction to pondere the CVC computing only 1 HS at a time"""
    # this version is necessary especially in the hole case to keep the CVC continuous
    phi = get_phi(kx, ky, s_coeff, False)
    phi = np.remainder(phi, 2*pi)
    sigma = sig(s_coeff)
    k = 1
    if T != -1 :
        k = get_CVC_T_dep(T)
        sigma = sig(s_coeff)*get_CVC_T_dep(T)
    if T<=0 :T=1
    eps= np.exp(-(phi-HS)**2/(2*sigma**2))
    return eps*k

@njit()
def get_CVC_T_dep(T):
    """Temperature dependance of epsilon"""
    # return 1 +T*0 #in case we want it to be constant
    return  1/(1+np.exp(-0.06*(T-10)))# a = 5

@njit()
def get_v(kx, ky, s_coeff, T = -1, FLEX = False, mode = "spin"):
    """Compute the velocities"""
    return get_v_not_modulo(kx, ky, s_coeff, T, FLEX, mode)

    vx0, vy0 = get_v0(kx, ky, s_coeff)

    if not FLEX:
        return vx0, vy0

    else:
        a, b, c = s_coeff[:3] # lattice constants

        alpha = get_epsilon(kx, ky, s_coeff, T, mode)  # Calculate alpha_k
        # alpha = get_alpha(kx, ky, s_coeff)  # Calculate alpha_k

        vxp1, vyp1 = get_v0(-ky, -kx, s_coeff)
        vxp2, vyp2 = get_v0(ky, kx, s_coeff)
        vxp, vyp = vxp1*((kx*ky)>=0)+vxp2*((kx*ky)<0), vyp1*((kx*ky)>=0)+vyp2*((kx*ky)<0) # points at opposite sides of the BZ
        # kxp, kyp = kx.copy(), ky.copy()
        # kxp[kxp >=0] = kxp[kxp >=0] - np.pi/a
        # kxp[kxp < 0] = kxp[kxp < 0] + np.pi/a
        # kyp[kyp < 0] = kyp[kyp < 0] + np.pi/b
        # kyp[kyp < 0] = kyp[kyp < 0] + np.pi/b
        # vxp, vyp = get_v0(kxp,kyp,s_coeff)
        Jx, Jy =(vx0 + alpha*vxp), (vy0 + alpha*vyp)#*1/(1-alpha**2)
        return Jx, Jy

@njit()
def get_v_not_modulo(kx, ky, s_coeff, T = -1, FLEX = False, mode = "spin"):
    """Compute the velocities. Solves continuity issues encontered in get_v with epsilon (and not epsilon_not_modulo)"""

    vx0, vy0 = get_v0(kx, ky, s_coeff)

    if not FLEX:
        return vx0, vy0

    else:
        # a, b, c = s_coeff[:3] # lattice constants
        HS = get_HS(s_coeff)
        # kx_HS, ky_HS = get_k_HS(s_coeff)

        #compute in different quandrants (0, pi/2,, etc) the value of epsilon, and their influence on points
        # neighboring the hot spots to garanty the continuity of the mean free path in adjacent hotspots not in the same quandrants
        # needed for hole doped cuprates such as LSCO
        #
        eps_q1 = get_epsilon_not_modulo(kx, ky, s_coeff, HS, T, mode)
        eps_q1 += get_epsilon_not_modulo(kx, ky, s_coeff, np.pi/2 - HS, T, mode)
        eps_q2 = get_epsilon_not_modulo(kx, ky, s_coeff, np.pi/2 + HS, T, mode)
        eps_q2 += get_epsilon_not_modulo(kx, ky, s_coeff, np.pi - HS, T, mode)
        eps_q3 = get_epsilon_not_modulo(kx, ky, s_coeff, np.pi + HS, T, mode)
        eps_q3 += get_epsilon_not_modulo(kx, ky, s_coeff, 3*np.pi/2 - HS, T, mode)
        eps_q4 = get_epsilon_not_modulo(kx, ky, s_coeff, 3*np.pi/2 + HS, T, mode)
        eps_q4 += get_epsilon_not_modulo(kx, ky, s_coeff,- HS, T, mode)
        # eps_q2 = 0
        vxp1, vyp1 = get_v0(-ky, -kx, s_coeff) #les vitesses scattered dans le quadrant 1
        vxp2, vyp2 = get_v0(ky, kx, s_coeff) #les vitesses scattered dans le quadrant 2
        # vxp, vyp = vxp1*((kx*ky)>=0)+vxp2*((kx*ky)<0), vyp1*((kx*ky)>=0)+vyp2*((kx*ky)<0) # points at opposite sides of the BZ

        Jx, Jy =(vx0 + (eps_q1 + eps_q3)*vxp1 + (eps_q2 + eps_q4)*vxp2), (vy0 + (eps_q1 + eps_q3)*vyp1 + (eps_q2 + eps_q4)*vyp2) #*1/(1-alpha**2) <-- could be pondered by a factor 1/(1-epsilon²)
        return Jx, Jy
## Eventually, this last set of routines computes the conductivity tensor using SCTIF formalism for the BTE. Then to obtain the resistivity, this tensor is inverted.
##
@njit()
def compute_time_integral(kxF, kyF, B, T, j, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, hole_like, FLEX):
    """Compute the time integral in the Shockley-Chambers formula for the conductivity. The integral
    is evaluated using trapezoidal rule, while the quasi-particle move along the FS in a cyclotron orbit."""

    N_t = int(t_lim/dt)+1

    kx_t, ky_t = evolve_kF(kxF, kyF, -B, dt, N_t, s_coeff, in_AFBZ, FLEX)
    vx_t, vy_t = get_v(kx_t, ky_t, s_coeff, T, FLEX)
    vj = vx_t if j == 'x' else vy_t


    ## The exponential time integral is computed with the trapezoidal rule simultaneously with the other integral
    ## since the cyclotron motion is the same (the timestep is the same).

    tau_inv = get_tau_inv(kx_t, ky_t, T, s_coeff, tau_choice, tau_coeff, hole_like) # scattering rate
    exp_int = np.zeros(N_t)
    for i in range(1,N_t):
        exp_int[i] = exp_int[i-1]+tau_inv[i-1]+tau_inv[i]

    return np.trapezoid(y=vj*np.exp(-exp_int*dt/2), dx=dt)

@njit()
def compute_sigma_ij(B, T, ij, N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    """Compute aget_v0 single term of the conductivity tensor using Shockley-Chambers formula. The FS integral
    is evaluated using trapezoidal rule, while the quasi-particle move along the FS in a cyclotron orbit."""

    a, b, c = s_coeff[:3] # lattice constants

    kxF, kyF = get_kF(s_coeff, N_k, in_AFBZ) # subdivide the FS with evenly spaced points
    vxF, vyF = get_v(kxF, kyF, s_coeff, T, FLEX)
    hole_like = False # 0. not in kxF # `False` if centered at Γ

    vi = vxF if ij[0] == 'x' else vyF
    vF = np.sqrt(vxF**2+vyF**2)

    ## Time integral
    time_int = np.zeros_like(kxF)
    for i in range(len(kxF)):
        time_int[i] = compute_time_integral(kxF[i], kyF[i], B, T, ij[1], dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, hole_like, FLEX)

    ## FS integral
    integrand = (vi/vF)*time_int
    # print(f"L'INTEGRAND EST {integrand} POUR {ij}")
    dkx, dky = np.diff(kxF), np.diff(kyF)
    dkF = np.sqrt(dkx**2+dky**2)

    ## Trapezoidal rule for non-uniform grid spacing
    tot = 0
    if not hole_like:
        for i in range(1, len(kxF)):
            tot += (integrand[i]+integrand[i-1])*dkF[i-1]/2
    elif hole_like:
        for i in range(1, len(kxF[:N_k])): tot += (integrand[:N_k][i]+integrand[:N_k][i-1])*dkF[:N_k][i-1]/2
        for i in range(1, len(kxF[N_k:])): tot += (integrand[N_k:][i]+integrand[N_k:][i-1])*dkF[N_k:][i-1]/2

    prefactor = e**2/(2*pi**2*c*hbar)
    return prefactor * 2 * tot # a factor of two is added to account also for negative ky


@njit()
def compute_resistivity(B, T, N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    """Compute the resistivity (in μΩcm) by inverting the conductivity tensor."""

    ## Only two components are needed, while the others follow from symmetry arguments.
    sigma_xx = compute_sigma_ij(B, T, 'xx', N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    sigma_xy = compute_sigma_ij(B, T, 'xy', N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    det_sigma = sigma_xx**2+sigma_xy**2

    rho_xx = sigma_xx/det_sigma/1e-8 # in μΩcm
    rho_xy = -sigma_xy/det_sigma/1e-8 # in μΩcm

    return rho_xx, rho_xy, -rho_xy, rho_xx


@njit()
def get_tau_iso(T, s_coeff, tau_choice, tau_coeff, in_AFBZ):
    """Return the minimum value of the scattering time. Useful for estimate the t_lim parameters in the computation."""

    kx, ky = get_kF(s_coeff, 11, in_AFBZ)
    # hole_like = 0. not in kx
    hole_like = False
    return np.min(get_tau_inv(kx, ky, T, s_coeff, tau_choice, tau_coeff, hole_like)**(-1))

def compute_B_sweep(B, TB, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    rho_B = np.zeros((4, len(B), len(TB)))
    for i in range(len(B)):
        for j in range(len(TB)):
            tau = get_tau_iso(TB[j], s_coeff, tau_choice, tau_coeff, in_AFBZ)
            t_lim = tau*N_tau
            rho_B[0,i,j], rho_B[1,i,j], rho_B[2,i,j], rho_B[3,i,j] = \
                compute_resistivity(B[i], TB[j], N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    return rho_B

def compute_error_B_sweep(B, TB, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    rho_err_xx, rho_err_RH = np.zeros(6), np.zeros(6)
    tau = get_tau_iso(TB[0], s_coeff, tau_choice, tau_coeff, in_AFBZ)
    t_lim = tau*N_tau

    # absolute error in FS spacing grid
    rho_err_xx[0], rho_err_RH[0], _, _ = \
    compute_resistivity(B[int(B[0]==0)], TB[0], 3*N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    rho_err_xx[3], rho_err_RH[3], _, _ = \
    compute_resistivity(B[-1], TB[0], 3*N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    # absolute error in the timestep
    rho_err_xx[1], rho_err_RH[1], _, _ = \
    compute_resistivity(B[int(B[0]==0)], TB[0], N_k, dt/2, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    rho_err_xx[4], rho_err_RH[4], _, _ = \
    compute_resistivity(B[-1], TB[0], N_k, dt/2, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    # absolute error in the truncation of the time integral
    rho_err_xx[2], rho_err_RH[2], _, _ = \
    compute_resistivity(B[int(B[0]==0)], TB[0], N_k, dt, 5*t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    rho_err_xx[5], rho_err_RH[5], _, _ = \
    compute_resistivity(B[-1], TB[0], N_k, dt, 5*t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    return rho_err_xx, rho_err_RH

def B_sweep(B, TB, N_k, dt, N_tau, out_fileB, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    print('Computing magnetoresistance...')
    t0 = time()
    rho_B = compute_B_sweep(B, TB, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    tot_t = time()-t0
    print(f'  Elapsed time: {strftime("%Hh:%Mm:%Ss", gmtime(tot_t))}\n')

    # print('Estimating the errors...')
    # rho_err_xx, rho_err_RH = compute_error_B_sweep(B, TB, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    # rho_err_xx[0] = abs(rho_err_xx[0]-rho_B[0,int(B[0]==0),0])
    # rho_err_xx[1] = abs(rho_err_xx[1]-rho_B[0,int(B[0]==0),0])
    # rho_err_xx[2] = abs(rho_err_xx[2]-rho_B[0,int(B[0]==0),0])
    # rho_err_xx[3] = abs(rho_err_xx[3]-rho_B[0,-1,0])
    # rho_err_xx[4] = abs(rho_err_xx[4]-rho_B[0,-1,0])
    # rho_err_xx[5] = abs(rho_err_xx[5]-rho_B[0,-1,0])
    # rho_err_RH[0] = abs(rho_err_RH[0]-rho_B[1,int(B[0]==0),0])*10/B[int(B[0]==0)]
    # rho_err_RH[1] = abs(rho_err_RH[1]-rho_B[1,int(B[0]==0),0])*10/B[int(B[0]==0)]
    # rho_err_RH[2] = abs(rho_err_RH[2]-rho_B[1,int(B[0]==0),0])*10/B[int(B[0]==0)]
    # rho_err_RH[3] = abs(rho_err_RH[3]-rho_B[1,-1,0])*10/B[-1]
    # rho_err_RH[4] = abs(rho_err_RH[4]-rho_B[1,-1,0])*10/B[-1]
    # rho_err_RH[5] = abs(rho_err_RH[5]-rho_B[1,-1,0])*10/B[-1]
    # print(f'  ρ_xx absolute error (in μΩcm) at the lowest field at T={TB[0]}K ', end='')
    # print(f'- dk: {rho_err_xx[0]:.2E}, dt: {rho_err_xx[1]:.2E}, tlim: {rho_err_xx[2]:.2E}')
    # print(f'  ρ_xx absolute error (in μΩcm) at the highest field at T={TB[0]}K ', end='')
    # print(f'- dk: {rho_err_xx[3]:.2E}, dt: {rho_err_xx[4]:.2E}, tlim: {rho_err_xx[5]:.2E}\n')
    # print(f'  R_H absolute error (in mm^3/C) at the lowest field at T={TB[0]}K ', end='')
    # print(f'- dk: {rho_err_RH[0]:.2E}, dt: {rho_err_RH[1]:.2E}, tlim: {rho_err_RH[2]:.2E}')
    # print(f'  R_H absolute error (in mm^3/C) at the highest field at T={TB[0]}K ', end='')
    # print(f'- dk: {rho_err_RH[3]:.2E}, dt: {rho_err_RH[4]:.2E}, tlim: {rho_err_RH[5]:.2E}\n')

    np.savez(out_fileB, B=B, rho_B=rho_B, T=TB)
    print(f'The output has been written to the file {out_fileB}.npz\n')
    return rho_B

def compute_T_sweep(T, BT, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    rho_T = np.zeros((4, len(BT), len(T)))
    for j in range(len(T)):
        tau = get_tau_iso(T[j], s_coeff, tau_choice, tau_coeff, in_AFBZ)
        t_lim = tau*N_tau
        for i in range(len(BT)):
            rho_T[0,i,j], rho_T[1,i,j], rho_T[2,i,j], rho_T[3,i,j] = \
                compute_resistivity(BT[i], T[j], N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    return rho_T

def compute_error_T_sweep(T, BT, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    rho_err_xx, rho_err_RH = np.zeros(6), np.zeros(6)
    tau = get_tau_iso(T[0], s_coeff, tau_choice, tau_coeff, in_AFBZ)
    t_lim = tau*N_tau

    # absolute error in FS spacing grid
    rho_err_xx[0], rho_err_RH[0], _, _ = \
    compute_resistivity(BT[-1], T[0], 3*N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    # absolute error in the timestep
    rho_err_xx[1], rho_err_RH[1], _, _ = \
    compute_resistivity(BT[-1], T[0], N_k, dt/2, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    # absolute error in the truncation of the time integral
    rho_err_xx[2], rho_err_RH[2], _, _ = \
    compute_resistivity(BT[-1], T[0], N_k, dt, 5*t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    tau = get_tau_iso(T[-1], s_coeff, tau_choice, tau_coeff, in_AFBZ)
    t_lim = tau*N_tau

    rho_err_xx[3], rho_err_RH[3], _, _ = \
    compute_resistivity(BT[-1], T[-1], 3*N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    rho_err_xx[4], rho_err_RH[4], _, _ = \
    compute_resistivity(BT[-1], T[-1], N_k, dt/2, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    rho_err_xx[5], rho_err_RH[5], _, _ = \
    compute_resistivity(BT[-1], T[-1], N_k, dt, 5*t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    return rho_err_xx, rho_err_RH

def T_sweep(T, BT, N_k, dt, N_tau, out_fileT, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    print('Computing resistivity vs T...')
    t0 = time()
    rho_T = compute_T_sweep(T, BT, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    tot_t = time()-t0
    print(f'  Elapsed time: {strftime("%Hh:%Mm:%Ss", gmtime(tot_t))}\n')

    # if BT[-1] != 0:
        # print('Estimating the errors...')
        # rho_err_xx, rho_err_RH = compute_error_T_sweep(T, BT, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
        # rho_err_xx[0] = abs(rho_err_xx[0]-rho_T[0,-1,0])
        # rho_err_xx[1] = abs(rho_err_xx[1]-rho_T[0,-1,0])
        # rho_err_xx[2] = abs(rho_err_xx[2]-rho_T[0,-1,0])
        # rho_err_xx[3] = abs(rho_err_xx[3]-rho_T[0,-1,-1])
        # rho_err_xx[4] = abs(rho_err_xx[4]-rho_T[0,-1,-1])
        # rho_err_xx[5] = abs(rho_err_xx[5]-rho_T[0,-1,-1])
        # rho_err_RH[0] = abs(rho_err_RH[0]-rho_T[1,-1,0])*10/BT[-1]
        # rho_err_RH[1] = abs(rho_err_RH[1]-rho_T[1,-1,0])*10/BT[-1]
        # rho_err_RH[2] = abs(rho_err_RH[2]-rho_T[1,-1,0])*10/BT[-1]
        # rho_err_RH[3] = abs(rho_err_RH[3]-rho_T[1,-1,-1])*10/BT[-1]
        # rho_err_RH[4] = abs(rho_err_RH[4]-rho_T[1,-1,-1])*10/BT[-1]
        # rho_err_RH[5] = abs(rho_err_RH[5]-rho_T[1,-1,-1])*10/BT[-1]
        # print(f'  ρ_xx absolute error (in μΩcm) at the lowest temperature at B={BT[-1]}T ', end='')
        # print(f'- dk: {rho_err_xx[0]:.2E}, dt: {rho_err_xx[1]:.2E}, tlim: {rho_err_xx[2]:.2E}')
        # print(f'  ρ_xx absolute error (in μΩcm) at the highest temperature at B={BT[-1]}T ', end='')
        # print(f'- dk: {rho_err_xx[3]:.2E}, dt: {rho_err_xx[4]:.2E}, tlim: {rho_err_xx[5]:.2E}\n')
        # print(f'  R_H absolute error (in mm^3/C) at the lowest temperature at B={BT[-1]}T ', end='')
        # print(f'- dk: {rho_err_RH[0]:.2E}, dt: {rho_err_RH[1]:.2E}, tlim: {rho_err_RH[2]:.2E}')
        # print(f'  R_H absolute error (in mm^3/C) at the highest temperature at B={BT[-1]}T ', end='')
        # print(f'- dk: {rho_err_RH[3]:.2E}, dt: {rho_err_RH[4]:.2E}, tlim: {rho_err_RH[5]:.2E}\n')

    np.savez(out_fileT, B=BT, rho_T=rho_T, T=T)
    print(f'The output has been written to the file {out_fileT}\n')
    return rho_T
