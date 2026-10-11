"""Fichier contenant la classe qui permet de simuler plus facilement
notamment en stockant les constantes pour ne les modifier qu'une seule fois en exécutant"""

import Code_SCTIF.BTE_routines as btr
import Code_SCTIF.metals as met
import numpy as np
import scipy.constants as sc

pi = np.pi
metal = met.Metal('LCCO',p='0.18')
tau_coeff_sans_HS = {"timp" : 27.e12, "alpha1" : 0.65e12, "alpha2" : 0.0052e12, "s_or_c": 1}
tau_choice = "hot spots"
FLEX = True
B0 = np.linspace(0, 35, 40)
TB0 = np.array( [20, 30, 40, 50])

class Simulation():
    """class used to set the simulation parameters in order
    to modify them just one before running a simulation"""
    def __init__(self,
        metal = metal,
        tau_coeff_sans_HS = tau_coeff_sans_HS,
        tau_choice = tau_choice,
        in_AFBZ = False,
        FLEX = FLEX):

        self.metal = metal # type of metal, object from class Metal
        self.HS = metal.get_HS() # location of the hot spots
        self.tau_coeff = np.array([param for param in tau_coeff_sans_HS.values()] + [self.HS]) #hot spots computed after the dictionnary containing the tau_coeff values (for easier readability) and adding at last the loation of the hot spots
        self.tau_choice = tau_choice #defined in the preambule
        self.in_AFBZ = in_AFBZ #idem
        self.FLEX = FLEX #idem

    def sweep_B(self,
        B = B0,
        TB = TB0,
        N_k = 200,
        dt = 1e-15,
        N_tau = 70):

        """ B : values of magnetic field sweep
        TB : values of temperature for each sweep
        N_k = 50 number of k-points of a single Fermi arc
        dt = 1e-15 timestep (in s)
        N_tau = 70 the upper limit in the time integral is N_tau times the scattering time"""

        s_coeff = self.metal.get_s_coeff()
        out_file_B = self.get_file_name("B",B,"T",TB)
        rho_xx, rho_xy, rho_yx, rho_yy = btr.B_sweep(
            B, TB, N_k, dt, N_tau, out_file_B, s_coeff, self.tau_choice, self.tau_coeff, self.in_AFBZ, self.FLEX
        )
        return rho_xx, rho_xy, rho_yx, rho_yy

    def sweep_T(self,
        T = np.logspace(0, 2, 50),
        BT = np.array([0,5,10,14]),
        N_k = 50,
        dt = 1e-15,
        N_tau = 70,
        out_file_T = ''):

        """ B : values of magnetic field sweep
        TB : values of temperature for each sweep
        N_k = 50 number of k-points of a single Fermi arc
        dt = 1e-15 timestep (in s)
        N_tau = 70 the upper limit in the time integral is N_tau times the scattering time"""

        s_coeff = self.metal.get_s_coeff()
        if out_file_T == "" : out_file_T = self.get_file_name("T",T,"B",BT)
        rho_xx, rho_xy, rho_yx, rho_yy = btr.T_sweep(
            T, BT, N_k, dt, N_tau, out_file_T, s_coeff, self.tau_choice, self.tau_coeff, self.in_AFBZ, self.FLEX
        )
        return rho_xx, rho_xy, rho_yx, rho_yy

    def compute_T_sweep_conductivity(self,
        ij,
        B =np.array([1,2,3,5,10,15]),
        T = np.linspace(1, 100, 50),
        N_k = 50,
        dt = 1e-15,
        N_tau = 70
        ):
        sigma = np.zeros((len(B), len(T)))
        for t in range(len(T)):
            tau = btr.get_tau_iso(T[t], self.metal.s_coeff, self.tau_choice, self.tau_coeff, self.in_AFBZ)
            t_lim = tau*N_tau
            for b in range(len(B)):
                sigma[b,t] = btr.compute_sigma_ij(B[b], T[t], ij, N_k, dt, t_lim, self.metal.s_coeff, self.tau_choice, self.tau_coeff, self.in_AFBZ, self.FLEX)
        return sigma

    def get_file_name(self,sweep,sweep_v,param,param_v):
        # param_v = list(param_v)
        if int(self.tau_coeff[-2]) == 0 :
            mode = "charge"
        else :
            mode = "spin"
        tag_mesure = f"{sweep}={int(np.min(sweep_v))}_to_{int(np.max(sweep_v))}_at_{param}={param_v}_FLEX={self.FLEX}"
        return f"./SCTIF/npz/{self.metal.name}_{self.tau_choice}_{mode}_"+tag_mesure
