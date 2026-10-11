# Tau(phi) plot
from math import gcd
from .BTE_routines import *
import matplotlib.pyplot as plt

#plot datas VS phi with a nice automatised x axis in gradient
def plot_tau_phi (ax,s_coeff, tau_choice, tau_coeff, T,N_k=100, in_AFBZ = False, **kwargs):
    """
    Plot the scattering rate 1/tau (in ps-1) against the azimuthal angle phi.

    PARAMETERS
    ----------
    s_coeff : tight-binding and geometrical parameters of the lattice
    tau_coeff, tau_choice: modelisation of the scattering rate and its parameters
    N_k (opt) : the number of points of the grid in k-space (default 101).
    in_AFBZ : whether the BZ is the Antiferomagnetic
    """

    kxF, kyF = get_kF(s_coeff, N_k, in_AFBZ)
    hole_like = (0. not in kxF)
    phi = get_phi(kxF, kyF, s_coeff, hole_like)
    phi = phi + np.pi*(hole_like)
    tau_phi = get_tau_inv(kxF,kyF,T, s_coeff, tau_choice, tau_coeff, hole_like)

    if hole_like:
        ax.plot(phi[:N_k], tau_phi[:N_k]/1e12, **kwargs)
        ax.plot(phi[N_k:], tau_phi[N_k:]/1e12)
    else :
        ax.plot(phi, tau_phi/1e12, **kwargs)

    set_x_ax(ax, 'phi',lim=1)
    ax.set_ylabel(r'$\tau^{-1}$ $(ps^{-1})$')
    ax.set_ylim(0,150)


def set_x_ax( ax, xl, lim=1):
    frac = [(j//gcd(j,4), 4//gcd(j,4)) for j in range(lim*4+1) ]
    frac = [(str(f[0]),str(f[1])) for f in frac]
    frac = [(""*(f[0]=="1") + f[0]*(f[0]!="1"), ""*(f[1]=="1") + ("/"+f[1])*(f[1]!="1")) for f in frac]
    if xl == 'phi':
        ax.set(xlabel=r'$\phi$', xlim=(0, np.pi),
            xticks=[np.pi/4*n for n in range(4*lim+1)], xticklabels=["0"] + [fr"{n}$\pi${m}" for (n,m) in frac[:-1] if n!= "0"] + [""]) #,r'$5\pi/4$',r'$3\pi/2$',r'$7\pi/4$',r'$2\pi$'])
                # , 5*np.pi/4, 3*np.pi/2,7*np.pi/4, 2*np.pi]

def plot_phi (ax, phi,y, hole_like = False, N_k=100, **kwargs):
    if False:
        ax.plot(phi[:N_k], y[:N_k], **kwargs)
        ax.plot(phi[N_k:], y[N_k:], **kwargs)
    else :
        ax.plot(phi, y, **kwargs)
    set_x_ax(ax, "phi", lim=1)

def make_dotlines(ax, sim, HS= True, zero = True):
    if zero:
        ax.axhline(y=0, ls="--",c="r", lw = 1)
    if HS:

        ax.axvline(x=sim.HS, ls="--", lw = 1, c= "C0")
        ax.axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1, c= "C0")
        ax.axvline(x= np.pi/2 + sim.HS, ls="--", lw = 1, c= "C0")
        ax.axvline(x= np.pi - sim.HS, ls="--", lw = 1, c= "C0")

def make_legend(ax, fig):
    """Draw the legend on the figure."""
    ncol = 1#
    fig.legend(ncol=ncol, loc='outside right center', fontsize=12)
