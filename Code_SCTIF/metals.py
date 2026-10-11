import numpy as np
from .BTE_routines import get_E, get_kF, get_v, get_DOS, get_phi, evolve_kF
from scipy.optimize import curve_fit
from scipy.constants import e, m_e, k as kB, hbar, mu_0
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

pi = np.pi

metals = {'Tl2201': {'a': 3.87, 'b': 3.87, 'c': 11.60,
               'p=0.12' : {'mu': 0.1998, 't1':0.7250, 't2': 0.3020, 't3': 0.0159, 't4': -0.0805, 't5': 0.0034}, #e-lie sign change
               'p=0.25' : {'mu': 0.2382, 't1':0.7250, 't2': 0.3020, 't3': 0.0159, 't4': -0.0805, 't5': 0.0034}},

          'LSCO' : {'a': 3.76, 'b': 3.76, 'c': 6.61,
               'p=0.12' : {'mu': 0.1952, 't1':-1.0000, 't2': 0.1656, 't3':-0.0828, 't4': 0.0800, 't5': 0.0000},
               'p=0.16' : {'mu': 0.2025, 't1':-1.0000, 't2': 0.1500, 't3':-0.0750, 't4': 0.0800, 't5': 0.0000}, #-like sign change
               'p=0.20' : {'mu': 0.2118, 't1':-1.0000, 't2': 0.1376, 't3':-0.0688, 't4': 0.0800, 't5': 0.0000},
               'p=0.23' : {'mu': 0.2201, 't1':-1.0000, 't2': 0.1304, 't3':-0.0652, 't4': 0.0800, 't5': 0.0000},
               'p=0.24' : {'mu': 0.2231, 't1':-1.0000, 't2': 0.1284, 't3':-0.0642, 't4': 0.0800, 't5': 0.0000},
               'p=0.28' : {'mu': 0.2365, 't1':-1.0000, 't2': 0.1224, 't3':-0.0612, 't4': 0.0800, 't5': 0.0000},
               'p=0.32' : {'mu': 0.2519, 't1':-1.0000, 't2': 0.1195, 't3':-0.0597, 't4': 0.0800, 't5': 0.0000}},

          'Nd-LSCO' : {'a': 3.75, 'b': 3.75, 'c': 6.60,
               'p=0.24' : {'mu': 0.1319, 't1':-0.6400, 't2': 0.0873, 't3':-0.0437, 't4': 0.0000, 't5': 0.0000}},

          'Bi2201' : {'a': 3.79, 'b': 3.79, 'c': 12.31,
               'p=0.08' : {'mu': 0.1338, 't1':-0.8800, 't2': 0.1373, 't3':-0.1439, 't4': 0.0573, 't5': 0.0000},
               'p=0.12' : {'mu': 0.1600, 't1':0.8800, 't2': 0.1373, 't3':-0.1439, 't4': -0.0573, 't5': 0.0000}, #e_like sign change
               'p=0.16' : {'mu': 0.1846, 't1':-0.8800, 't2': 0.1373, 't3':-0.1439, 't4': 0.0573, 't5': 0.0000},
               'p=0.20' : {'mu': 0.2074, 't1':0.8800, 't2': 0.1373, 't3':-0.1439, 't4': -0.0573, 't5': 0.0000},#e_like sign change
               'p=0.24' : {'mu': 0.2283, 't1':-0.8800, 't2': 0.1373, 't3':-0.1439, 't4': 0.0573, 't5': 0.0000},
               'p=0.28' : {'mu': 0.2472, 't1':-0.8800, 't2': 0.1373, 't3':-0.1439, 't4': 0.0573, 't5': 0.0000}},
          'LCCO' : {'a': 4.01, 'b': 4.01, 'c': 6.2, #e-like sign change for all of them
            'p=0.100' : {'mu': -0.045,  't1':-1.2, 't2':-0.264, 't3':0.132, 't4': 0.0000, 't5': 0.0000},
            'p=0.125' : {'mu': -0.02059,'t1':-1.2, 't2':-0.27, 't3':0.135, 't4': 0.0000, 't5': 0.0000},
            'p=0.14' : {'mu': -0.0132, 't1':-1.2, 't2': -0.2736, 't3':0.1368, 't4': 0.0000, 't5': 0.0000},
            'p=0.15' : {'mu': 0., 't1':-1.2, 't2': -0.276, 't3':0.1380, 't4': 0.0000, 't5': 0.0000},
            'p=0.16' : {'mu': 0.0132, 't1':-1.2, 't2': -0.2784, 't3':0.1392, 't4': 0.0000, 't5': 0.0000},
            'p=0.175' : {'mu': 0.033,   't1':-1.2, 't2':-0.282, 't3':0.141, 't4': 0.0000, 't5': 0.0000},
            'p=0.18' : {'mu': 0.0396, 't1':-1.2, 't2': -0.2832, 't3':0.1416, 't4': 0.0000, 't5': 0.0000},
            'p=0.200' : {'mu': 0.0615,  't1':-1.2, 't2': -0.288, 't3':0.144, 't4': 0.0000, 't5': 0.0000} },
          'LCCOHole' : {'a': 4.01, 'b': 4.01, 'c': 6.2,
              'p=0.150' : {'mu': 0., 't1':1.2, 't2': -0.276, 't3':0.138, 't4': 0.0000, 't5': 0.0000}}}

class Metal():
    """
    Class used to plot useful quantities of the system, as the energy dispersion,
    the Fermi surface, the velocities and the density of states.
    """

    def __init__(self, met, p):
        try:
            self.met = metals[met]
        except KeyError as error_msg:
            err_str = f"The metal {error_msg} is not implemented yet: possible values are {list(metals.keys())}"
            raise KeyError(err_str)
        self.name = met+"_"+str(p)
        self.p = str(p)
        self.x = p
        self.a, self.b, self.c = self.get_lattice_vectors()
        self.s_coeff = self.get_s_coeff()

    def get_lattice_vectors(self):
        return self.met['a']/1e10, self.met['b']/1e10, self.met['c']/1e10

    def get_TB_params(self):
        try:
            met = self.met['p='+self.p]
        except KeyError as error_msg:
            err_str = f"The doping value {error_msg} is not implemented yet: possible values are {list(self.met.keys())[3:]}"
            raise KeyError(err_str)


        return met['mu']*e, met['t1']*e, met['t2']*e, met['t3']*e, met['t4']*e, met['t5']*e

    def get_s_coeff(self):
        a, b, c = self.get_lattice_vectors()
        mu, t1, t2, t3, t4, t5 = self.get_TB_params()
        s_coeff = np.array([a, b, c, mu, t1, t2, t3, t4, t5])
        return s_coeff

    def get_HS(self):
        """compute the angle of hotspot (i.e. intersection between FS and AFBZ) in [0, pi/4]"""
        a,b = self.s_coeff[:2]
        kxF, kyF = get_kF(self.s_coeff, 1000, False)
        ky_AFBZ = pi/b - kxF*a/b
        i_HS = np.argmin(np.abs(np.sqrt(kyF**2+kxF**2)-np.sqrt(ky_AFBZ**2+kxF**2)))
        kxF_HS, kyF_HS = kxF[i_HS], kyF[i_HS]
        return np.arctan(kyF_HS/kxF_HS)

    def set_x_ax(self, ax, xl):
        if xl == 'phi':
            ax.set(xlabel=r'$\phi$', xlim=(0, np.pi),
                   xticks=[0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi], xticklabels=['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])

        elif xl == 'kx':
            ax.set(xlabel=r'$k_x$', xlim=(-np.pi/self.a, np.pi/self.a),
                   xticks=[-np.pi/self.a, 0, np.pi/self.a], xticklabels=[r'$-\pi/a$', 0, r'$\pi/a$'])

    def plot_AFBZ(self, ax):
        """
        Draw the antiferromagnetic Brillouin zone (AFBZ) in k-space.
        """

        # The antiferromagnetic Brillouin Zone (AFBZ) is half of the standard BZ
        ax.plot([-np.pi/self.a, 0, np.pi/self.a, 0, -np.pi/self.a],
                [0, -np.pi/self.b, 0, np.pi/self.b, 0], '--', color='black', lw=0.6)

        ## PLOT ADJUSTMENT
        self.set_x_ax(ax, 'kx')
        ax.set_ylabel(r'$k_y$', labelpad=-10)
        ax.set(ylim=(-np.pi/self.b, np.pi/self.b), aspect='equal',
               yticks=[-np.pi/self.b, 0, np.pi/self.b], yticklabels=[r'$-\pi/b$', 0, r'$\pi/b$'])

    def plot_FS(self, ax, N_k=100, in_AFBZ=False, full_FS=True, cbar_mode=0, isoE = True, phi_min = 0, phi_max = 0):
        """
        Plot the energy dispersion (in eV) in a contour plot and draw the Fermi surface.

        PARAMETERS
        ----------

        N_k (opt) : the number of points of the grid in k-space (default 101).
        full_FS (opt) : if `True` draw the angle at the Gamma point and put the colorbar (default True).
        cbar_mode (opt): if `0` draw the color bar on the right side, while if `1` draw
                         the color bar on the top
        """

        self.plot_AFBZ(ax)

        kx = np.linspace(-np.pi/self.a, np.pi/self.a, N_k)
        ky = np.linspace(-np.pi/self.b, np.pi/self.b, N_k)

        kx_2D, ky_2D = np.meshgrid(kx, ky)
        kx, ky = np.copy(kx_2D), np.copy(ky_2D)
        if in_AFBZ:
            for i in range(N_k):
                for j in range(N_k):
                    if abs(kx_2D[i,j]) >= np.pi/self.b - abs(ky_2D[i,j]): kx[i,j], ky[i,j] = np.nan, np.nan

        ## PLOT THE DISPERSION RELATION AND 2D FERMI SURFACE (FS)
        energy = get_E(kx, ky, self.s_coeff)/e # compute the dispersion in eV

        if isoE :
            im = ax.contourf(kx, ky, energy, levels=100, cmap='summer')
            if in_AFBZ:
                ax.contourf(kx+np.pi/self.a, ky+np.pi/self.b, energy, levels=100, cmap='summer')
                ax.contourf(kx-np.pi/self.a, ky+np.pi/self.b, energy, levels=100, cmap='summer')
                ax.contourf(kx+np.pi/self.a, ky-np.pi/self.b, energy, levels=100, cmap='summer')
                ax.contourf(kx-np.pi/self.a, ky-np.pi/self.b, energy, levels=100, cmap='summer')

        kxF, kyF = get_kF(self.s_coeff, N_k, in_AFBZ)

        if in_AFBZ:
            ax.plot(kxF[kxF<0], kyF[kxF<0], c='black')
            ax.plot(kxF[kxF<0], -kyF[kxF<0], c='black')
            ax.plot(kxF[kxF>0]-np.pi/self.a, -kyF[kxF>0]+np.pi/self.b, c='black')
            ax.plot(kxF[kxF>0]-np.pi/self.a, kyF[kxF>0]-np.pi/self.b, c='black')

            ax.plot(kxF[kxF>0], kyF[kxF>0], c='black')
            ax.plot(kxF[kxF>0], -kyF[kxF>0], c='black')
            ax.plot(kxF[kxF<0]+np.pi/self.a, -kyF[kxF<0]+np.pi/self.b, c='black')
            ax.plot(kxF[kxF<0]+np.pi/self.a, kyF[kxF<0]-np.pi/self.b, c='black')
        else:
            ax.plot(kxF, kyF, c='black')
            ax.plot(kxF, -kyF, c='black')

        if full_FS :
            if isoE:
                if cbar_mode == 1:
                    cbar = plt.colorbar(im, ax=ax, location='top', pad=0.1, orientation='horizontal', fraction=0.075)
                else:
                    cbar = plt.colorbar(im, ax=ax, orientation='vertical')

                cbar_ticks = list(np.linspace(np.min(energy[np.isfinite(energy)]),
                                            np.max(energy[np.isfinite(energy)]), 5))
                cbar.set_label(label=r'$E\ (eV)$', size=10)
                cbar.set_ticks(cbar_ticks, labels=[f'{tcks:.2f}' for tcks in cbar_ticks], fontsize=10)
            if phi_min != phi_max :
                ax.plot([0, 2*np.pi/self.a*np.cos(phi_min)],[0, 2*np.pi/self.a*np.sin(phi_min)], "--")
                ax.plot([0, 2*np.pi/self.a*np.cos(phi_max)],[0, 2*np.pi/self.a*np.sin(phi_max)], "--")

            kx_phi = np.pi/self.a/3
            phi = np.pi/6 # 30°
            ky_phi = kx_phi*np.tan(phi)
            ax.plot([0,kx_phi], [0,ky_phi], c='black', lw=1)
            ax.plot([0,kx_phi*1.2], [0,0], '--', c='black', lw=1)
            arc_angles = np.linspace(0, phi, 20)
            arc_xs = kx_phi*0.7*np.cos(arc_angles)
            arc_ys = kx_phi*0.7*np.sin(arc_angles)
            ax.plot(arc_xs, arc_ys, color='black', lw=1)
            ax.text(kx_phi*0.8, ky_phi*0.3, r'$\phi$', fontsize=10)
            ax.text(-np.pi/self.a*0.4, np.pi/self.b*0.3, 'FS', fontsize=10)
            ax.text(np.pi/self.a*0.55, np.pi/self.b*0.55, 'AFBZ', fontsize=10)

        ax.scatter(0, 0, s=2, c='black')
        ax.text(0.45, 0.43, r'$\Gamma$', fontsize=10, transform=ax.transAxes)

    def plot_vels(self, ax, N_k=100, in_AFBZ=False, FLEX=False):
        """
        Plot the Fermi velocities (in km/s) against the azimuthal angle phi.

        PARAMETERS
        ----------

        N_k (opt) : the number of points of the grid in k-space (default 101).
        """

        kxF, kyF = get_kF(self.s_coeff, N_k, in_AFBZ)
        hole_like = (0. not in kxF)
        phi = get_phi(kxF, kyF, self.s_coeff, hole_like)
        phi = phi + np.pi*(hole_like)

        vxF, vyF = get_v(kxF, kyF, self.s_coeff, FLEX)

        if hole_like:
            ax.plot(phi[:N_k], vxF[:N_k]/1e3, '--', c='blue', lw=1.25, label=r'$v_x$')
            ax.plot(phi[:N_k], vyF[:N_k]/1e3, '--', c='red', lw=1.25, label=r'$v_y$')
            ax.plot(phi[:N_k], np.sqrt(vxF[:N_k]**2 + vyF[:N_k]**2)/1e3, c='green', label=r'$v_F$')

            ax.plot(phi[N_k:], vxF[N_k:]/1e3, '--', c='blue', lw=1.25)
            ax.plot(phi[N_k:], vyF[N_k:]/1e3, '--', c='red', lw=1.25)
            ax.plot(phi[N_k:], np.sqrt(vxF[N_k:]**2 + vyF[N_k:]**2)/1e3, c='green')
        else:
            ax.plot(phi, vxF/1e3, '--', c='blue', lw=1.25, label=r'$v_x$')
            ax.plot(phi, vyF/1e3, '--', c='red', lw=1.25, label=r'$ky_F$')
            ax.plot(phi, np.sqrt(vxF**2 + vyF**2)/1e3, c='green', label=r'$v_F$')

        # ax.set_ylim(400,650)
        self.set_x_ax(ax, 'phi')
        ax.set_ylabel(r'$v_F \ (km/s)$')
        ax.legend()

    def plot_kf(self, ax, N_k=100, in_AFBZ = False, FLEX = False):
        """
        Plot the Fermi wavevectord (in nm-1) against the azimuthal angle phi.

        PARAMETERS
        ----------

        N_k (opt) : the number of points of the grid in k-space (default 101).
        """

        kxF, kyF = get_kF(self.s_coeff, N_k, in_AFBZ)
        hole_like = (0. not in kxF)
        phi = get_phi(kxF, kyF, self.s_coeff, hole_like)
        phi = phi + np.pi*(hole_like)

        if hole_like:
            ax.plot(phi[:N_k], kxF[:N_k]/1e9, '--', c='blue', lw=1.25, label=r'$kxF$')
            ax.plot(phi[:N_k], kyF[:N_k]/1e9, '--', c='red', lw=1.25, label=r'$kyF$')
            ax.plot(phi[:N_k], np.sqrt(kxF[:N_k]**2 + kyF[:N_k]**2)/1e9, c='green', label=r'$k_F$')

            ax.plot(phi[N_k:], kxF[N_k:]/1e9, '--', c='blue', lw=1.25)
            ax.plot(phi[N_k:], kyF[N_k:]/1e9, '--', c='red', lw=1.25)
            ax.plot(phi[N_k:], np.sqrt(kxF[N_k:]**2 + kyF[N_k:]**2)/1e9, c='green')
        else:
            ax.plot(phi, kxF/1e9, '--', c='blue', lw=1.25, label=r'$kx_F$')
            ax.plot(phi, kyF/1e9, '--', c='red', lw=1.25, label=r'$ky_F$')
            ax.plot(phi, np.sqrt(kxF**2 + kyF**2)/1e9, c='green', label=r'$k_F$')

        self.set_x_ax(ax, 'phi')
        ax.set_ylabel(r'$k_F \ (nm)$')
        ax.legend()

    def plot_DOS(self, ax, N_k=100, in_AFBZ=False, FLEX=False):
        """
        Plot the density of state (in 1/eV*Å) against the azimuthal angle phi.

        PARAMETERS
        ----------

        N_k (opt) : the number of points of the grid in k-space (default 101).
        """

        kxF, kyF = get_kF(self.s_coeff, N_k, in_AFBZ)
        hole_like = (0. not in kxF)
        phi = get_phi(kxF, kyF, self.s_coeff, hole_like)
        phi = phi + np.pi*(hole_like)

        dos = get_DOS(kxF, kyF, self.s_coeff)*e*1e-10 # 1/eV*Å

        if hole_like:
            ax.plot(phi[:N_k], dos[:N_k], c='black')
            ax.plot(phi[N_k:], dos[N_k:], c='black')
        else:
            ax.plot(phi, dos, c='black')

        self.set_x_ax(ax, 'phi')
        ax.set_ylabel(r'$1/|\nabla\varepsilon_k|\ (eV^{-1}Å^{-1})$')

    def plot_cycl_orbit(self, fig, ax, N_k=100, B=10, dt=2e-15, N_t=20000, in_AFBZ=False, repeat=False):
        self.plot_FS(ax, N_k=N_k, full_FS=False)
        kxF, kyF = get_kF(self.s_coeff, N_k=5, in_AFBZ=in_AFBZ)
        kx_t, ky_t = evolve_kF(kxF[0], kyF[0], B, dt, N_t, self.s_coeff, in_AFBZ=in_AFBZ)

        end_orbit = np.where(np.abs(kx_t-kx_t[0])/np.abs(kx_t[0])+np.array(ky_t>0)<=1e-3)[0][-1]

        pt, = ax.plot([], [], '^', c='black', ms=6.)
        nstep = 20
        def kx_evol(t):
            pt.set_data([kx_t[nstep*t]], [ky_t[nstep*t]])
            return pt,

        return FuncAnimation(fig, kx_evol, interval=5, blit=True, repeat=repeat, frames=end_orbit//nstep)

def get_TB_params_LSCO(p):
    a, b, c = 3.76e-10, 3.76e-10, 6.61e-10

    dopings = [0.16, 0.185, 0.1914, 0.21, 0.26]
    mu = [0.2025, 0.2080, 0.2095, 0.2145, 0.2295]
    t2 = [0.1500, 0.1420, 0.1400, 0.1350, 0.1250]

    plt.close('all')
    fig, ax = plt.subplots(1, 2, figsize=(9,3.5), constrained_layout=True)
    fig.canvas.header_visible = False

    ax[0].plot(dopings, mu, 'D', c='red', label=r'$\mu$')
    ax[0].plot(dopings, t2, 'D', c='blue', label=r'$t_2$')
    ax[0].legend(shadow=True, fancybox=True, fontsize=10)

    def quad_fit(x, a0, a1, a2):
        return a0 + a1*x + a2*x**2

    popt_mu, _ = curve_fit(quad_fit, dopings, mu)
    popt_t2, _ = curve_fit(quad_fit, dopings, t2)

    dopings = sorted(dopings+[p])
    dopings = np.linspace(dopings[0], dopings[-1], 100)

    ax[0].plot(dopings, quad_fit(dopings, *popt_mu), '--', lw=1., c='red')
    ax[0].plot(dopings, quad_fit(dopings, *popt_t2), '--', lw=1., c='blue')

    mu = quad_fit(p, *popt_mu)
    t1 = -1
    t2 = quad_fit(p, *popt_t2)
    t3 = -0.5*t2
    t4 = 0.08

    ax[0].plot(p, mu, 'D', c='white', mec='red')
    ax[0].plot(p, t2, 'D', c='white', mec='blue')

    ax[0].axvline(0.1914, ls='-.', lw=1., c='black')
    ax[0].text(0.195, 0.227, 'vHs', fontsize=10)

    kx = np.linspace(-np.pi/a, np.pi/a, 100)
    kx, ky = np.meshgrid(kx, kx)

    s_coeff = np.array([a, b, c, mu, t1, t2, t3, t4, 0.])
    energy = get_E(kx, ky, s_coeff)

    im = ax[1].contourf(kx, ky, energy, levels=100, cmap='summer')
    FS = ax[1].contour(kx, ky, energy, levels=[0], colors=['black'])

    cbar = plt.colorbar(im, ax=ax[1], location='right')
    cbar.set_label(label=r'$E\ (eV)$', size=10, labelpad=10)
    cbar_ticks = list(np.linspace(np.min(energy), np.max(energy), 5))
    cbar.set_ticks(cbar_ticks, labels=[f'{tcks:.2f}' for tcks in cbar_ticks], fontsize=10)

    ax[1].plot([-np.pi/a, 0, np.pi/a, 0, -np.pi/a], [0, -np.pi/b, 0, np.pi/b, 0], '--', color='black', lw=0.8)
    ax[1].text(0.76, 0.76, 'AFBZ', fontsize=11, transform=ax[1].transAxes,
           horizontalalignment='left', verticalalignment='bottom')

    ax[1].scatter(0, 0, s=2, c='black')
    ax[1].text(0.51, 0.51, r'$\Gamma$', fontsize=11, transform=ax[1].transAxes,
           horizontalalignment='left', verticalalignment='bottom')
    ax[1].text(0.05, 0.9, f'p={p}', fontsize=11, transform=ax[1].transAxes)

    ax[0].set(xlabel='p', ylabel='hopping energy (eV)')
    ax[1].set_ylabel(r'$k_y$', labelpad=-10)
    ax[1].set(xlabel=r'$k_x$', xlim=(-np.pi/a, np.pi/a), ylim=(-np.pi/b, np.pi/b), aspect='equal',
              xticks=[-np.pi/a, 0, np.pi/a], xticklabels=[r'$-\pi/a$', 0, r'$\pi/a$'],
              yticks=[-np.pi/b, 0, np.pi/b], yticklabels=[r'$-\pi/b$', 0, r'$\pi/b$'])
    plt.show()

    print(f'Tight-Binding parameterization (in eV) at {p} doping: μ={mu:.4f}, t1={t1:.4f}, t2={t2:.4f}, t3={t3:.4f}, t4={t4:.4f}')
    return mu, t1, t2, t3, t4

def Luttinger_count_Bi2201(mu):
    a, b, c = 3.79e-10, 3.79e-10, 12.31e-10
    t1, t2, t3, t4, t5 = -0.88, 0.1373, -0.1439, 0.0573, 0.

    s_coeff = np.array([a, b, c, mu*e, t1*e, t2*e, t3*e, t4*e, t5*e])

    A_BZ = (2*np.pi/a)*(2*np.pi/b)
    kxF, kyF = get_kF(s_coeff, N_k=501, in_AFBZ=False)
    if 0. not in kxF:
        s_coeff = np.array([a, b, c, mu*e, -t1*e, t2*e, t3*e, -t4*e, t5*e])
        kxF, kyF = get_kF(s_coeff, N_k=501, in_AFBZ=False)

    A_holes = 2*(np.trapz(kyF, x=kxF))
    p_lutt = (2*A_holes/A_BZ) - 1

    print(f'At μ={mu:.4f} eV, the inferred doping is p={p_lutt:.4f}')
    return p_lutt
