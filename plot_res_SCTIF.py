from Code_SCTIF.BTE_routines_SCTIF import *
from Code_SCTIF.metals_SCTIF import *
from Code_SCTIF.plot_results_SCTIF import *
import Code_SCTIF.plot_tau_phi_SCTIF as pt
import matplotlib.pyplot  as plt
from matplotlib.ticker import NullLocator
from simul_SCTIF import Simulation
import matplotlib as mpl
import pandas as pd

plt.style.use("pacoty.mplstyle")

LCCOH = Metal('LCCOHole', p='0.150')
LCCO = Metal('LCCO', p='0.150')
pi = np.pi
sim = Simulation()
out_LCCO_T = Results(f'./SCTIF/npz/{sim.metal.name}_{sim.tau_choice}_{int(sim.tau_coeff[-2])*"spin" + (1-int(sim.tau_coeff[-2]))*"charge"}_T=1_to_100_at_B=[ 1  2  3  5 10 15]_FLEX={sim.FLEX}.npz')
out_LCCO_B = Results(f'./SCTIF/npz/{sim.metal.name}_{sim.tau_choice}_{int(sim.tau_coeff[-2])*"spin" + (1-int(sim.tau_coeff[-2]))*"charge"}_B=0_to_35_at_T=[20 30 40 50]_FLEX={sim.FLEX}.npz')
# out_tau_HS_0 = Results(f'./tau_HS_est_nul.npz')
# out_tau_HS_1 = Results(f'./tau_HS_est_un.npz')
# out_tau_HS_T_dep = Results(f'./tau_HS_depend_T.npz')

# out_tau_HS_0_NO_CVC = Results(f'./tau_HS_est_nul_NO_CVC.npz')
# out_tau_HS_1_NO_CVC = Results(f'./tau_HS_est_un_NO_CVC.npz')
# out_tau_HS_T_dep_NO_CVC = Results(f'./tau_HS_T_dependant_NO_CVC.npz')

# out_tau_hs_0_cvc_not_t= Results('./tau_HS_est_nul_CVC_not_T_dependant.npz')
# out_tau_hs_1_cvc_not_t= Results('./tau_HS_est_un_CVC_not_T_dependant.npz')

# out_no_gamma_N_cvc_not_t_no_tau_hs = Results("./tau_HS_est_nul_gamma_N_0_CVC_not_T_dependant.npz")

# out_tau_hs_dans_beta = Results("./tau_HS_dans_beta_w_cvc.npz")
N=15
figs = []
names = []
for i in range(N) :
    figs.append(plt.subplots(figsize=(6, 4.5), constrained_layout=True))
    names.append("")
ax = [f[1] for f in figs]
fig = [f[0] for f in figs]

kx,ky = get_kF(sim.metal.s_coeff,200,False)
Jx, Jy = get_v(kx, ky, sim.metal.s_coeff, FLEX=sim.FLEX)
vx, vy = get_v0(kx, ky, sim.metal.s_coeff)

hole_like = False #(0. not in kx)
phi = get_phi(kx,ky,sim.metal.s_coeff, hole_like)
phi = phi + np.pi*hole_like
phi_min = -1
phi_max = 15


cmap = plt.cm.Reds
T=[0,10,20,30,40,50,60,70,80]
col = [cmap(i) for i in np.linspace(0, cmap.N, len(T)+2, dtype=np.int32)][2:]
i = 0
for Ti in T:
    gamma = get_tau_inv(kx,ky,Ti, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, hole_like = False)
    Jx, Jy = get_v(kx, ky, sim.metal.s_coeff, T = Ti, FLEX=sim.FLEX)
    lx, ly = Jx/gamma, Jy/gamma
    dly = np.gradient(ly,phi)
    pt.plot_phi(ax[0], phi, lx*dly, hole_like, c=col[i], label = f"${Ti}$ K")
    ax[0].set_ylabel(r"$l_x*dl_y/d\phi$")
    ax[0].axhline(y=np.sum(lx*dly/(2*np.pi)), ls="--",c=col[i], lw = 1)
    pt.plot_phi(ax[14], phi, gamma, hole_like, c=col[i], label = f"${Ti}$ K")
    ax[14].legend()
    ax[0].legend()
    i+=1

i=0
T = np.linspace(1,200, 50)
B = np.array([2]) #1,2,3,5,10,15

# dsigma_dphi = compute_dsigma_dphi(B[0], 100, "xy", 200, 1e-15, 10e-12, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, sim.in_AFBZ, sim.FLEX, phi_min=phi_min, phi_max=phi_max)
# pt.plot_phi(ax[0], phi, dsigma_dphi)

# ax[0].set_ylabel(r"$\mathrm{d}\sigma_{xy}/\mathrm{d}\phi$")
ax[0].axvline(x=sim.HS, ls="--", lw = 1)
ax[0].axhline(y=0, ls="--",c="r", lw = 1)
ax[0].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
ax[2].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
ax[2].axvline(x= sim.HS, ls="--", lw = 1)
ax[2].axvline(x= np.pi/2 + sim.HS, ls="--", lw = 1)
ax[2].axvline(x= np.pi - sim.HS, ls="--", lw = 1)
ax[0].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
ax[0].axvline(x=np.pi/2 + sim.HS, ls="--", lw = 1)
ax[0].axvline(x=np.pi - sim.HS, ls="--", lw = 1)

dJy = np.gradient(Jy/gamma, phi)
norm = mpl.colors.CenteredNorm()#Normalize(vmin = -1, vmax = 1) # = np.min(Jx*dJy), vmax = -np.min(Jx*dJy))
norm2 = mpl.colors.Normalize(vmin = -np.abs(np.min(Jx/gamma*dJy)), vmax = -np.abs(np.min(Jx/gamma*dJy))) # = np.min(Jx*dJy), vmax = -np.min(Jx*dJy))
color = Jx/gamma*dJy/np.abs(Jx/gamma*dJy)
pt.plot_phi(ax[2], phi, Jx, False, label = "$J_x$")
pt.plot_phi(ax[2], phi, vx, False, label = "$v_x$")
ax[2].set_ylabel(r"$v_x$")
ax[2].axhline(y=0, ls="--", lw=1, c="r")
# ax[2].scatter(phi, Jx, c =  color, norm= norm, cmap = mpl.colormaps["coolwarm"])

lx = Jy/gamma

# pt.plot_phi(ax[3], phi, lx, False, label = "$J_x$")
# pt.plot_phi(ax[3].twinx(), phi, dJy, False, label = "$v_x$")
# ax[3].set_ylabel(r"$l_x$")
# ax[3].axhline(y=0, ls="--", lw=1, c="r")
# ax[3].scatter(phi, lx, c =  color, norm= norm, cmap = mpl.colormaps["coolwarm"])


# pt.plot_phi(ax[4], phi, Jx/gamma*dJy, color, False, label = "$J_x$")

# ax[4].scatter(phi, Jx/gamma*dJy,c = color, norm= norm, cmap = mpl.colormaps["coolwarm"])
# ax[4].set_ylabel(r"$l_x*dl_y/\phi$")
# ax[4].axhline(y=0, ls="--", lw=1, c="r")


a = sim.metal.a
ax[1].plot([0, 2*np.pi/a*np.cos(sim.HS)],[0, 2*np.pi/a*np.sin(sim.HS)], "--", lw = 1, c="r")
ax[1].plot([0, 2*np.pi/a*np.cos(pi/2 - sim.HS)],[0, 2*np.pi/a*np.sin(np.pi/2 - sim.HS)], "--", lw = 1, c="r")

# ax[0].scatter(phi, dsigma_dphi, c =  color, norm= norm, cmap = mpl.colormaps["coolwarm"])

cmap = mpl.colormaps["twilight"]
cbar2 = plt.colorbar(mpl.cm.ScalarMappable(norm= norm,cmap = mpl.colormaps["coolwarm"]),ax=ax[2], orientation='vertical', label = r"$J_x\mathrm{d}J_y/\mathrm{d}\phi$")
# cbar0 = plt.colorbar(mpl.cm.ScalarMappable(norm= norm,cmap = mpl.colormaps["coolwarm"]),ax=ax[0], orientation='vertical', label = r"$J_x\mathrm{d}J_y/\mathrm{d}\phi$")
# cbar3 = plt.colorbar(mpl.cm.ScalarMappable(norm= norm,cmap = mpl.colormaps["coolwarm"]),ax=ax[3], orientation='vertical', label = r"$J_x\mathrm{d}J_y/\mathrm{d}\phi$")
# sign_sigma = dsigma_dphi

sim.metal.plot_FS(ax[1], isoE = False, phi_min=phi_min, phi_max=phi_min)
ax[1].quiver(kx,ky ,Jx/gamma, Jy/gamma, Jx*ky - Jy*kx, cmap = mpl.colormaps["coolwarm"], width = 0.002, scale =  0.0000001) #scale = 0.00000005 np.zeros_like(phi) phi, cmap=cmap
ax[1].quiver(kx,-ky ,Jx/gamma, -Jy/gamma, Jx*ky - Jy*kx, cmap = mpl.colormaps["coolwarm"], width = 0.002, scale =  0.0000001) #scale = 0.00000005 np.zeros_like(phi)
# ax[5].plot(kxt, kyt, marker = '.')

# sim.metal.plot_FS(ax0)4

# pt.make_dotlines(ax[6], sim, HS =False)
# eps = get_epsilon(kx, ky, sim.metal.s_coeff)
# pt.plot_phi(ax[5], phi, eps)
# # ax[5].set_ylabel(r"$\varepsilon(k)$")
# ax[5].axhline(y = 0, ls = '--', c='r', lw = '1')
# ax[5].axhline(y = 1, ls = '--', c='r', lw = '1')

og_data_rh = pd.read_csv('r_h_exp_data_x_0.159.csv')
T_exp = og_data_rh["Temperature"]
Rh_exp = og_data_rh["Rh"]

# out_tau_HS_T_dep.plot_RH_vs_T(ax[6], marker = '.')
# names[6] = "rho_xx_tau_HS_T_dep"
# ax[6].plot(T_exp, Rh_exp, c= 'C3', label = "exp results")
# ax[6].set(xscale = 'log', xlim = 1,title = r"LCCO with $\tau{HS}$ T dependant")
# pt.make_dotlines(ax[6],sim, HS= False)
# names[10] = out_LCCO_B.plot_RH_vs_B(ax[10])
# names[8] = out_LCCO_B.plot_MR(ax[8])
# out_LCCO.plot_MR(ax2)
# out_tau_HS_T_dep.plot_rho_ij_vs_T("xx",ax[7])
# names[7] = "rho_xx_tau_T_dep"
# ax[7].set(title = r"$\rho_{xx}$ with $\tau{HS}$ T dependant")
# names[9] = out_LCCO_B.plot_drho_db_vs_b(ax[9])

# out_tau_HS_T_dep.make_legend(fig[7])
# out_tau_HS_T_dep.make_legend(fig[6])

# out_LCCO_B.make_legend(fig[8])
# out_LCCO_B.make_legend(fig[9])
# out_LCCO_B.make_legend(fig[10])

# out_tau_hs_0_cvc_not_t.plot_RH_vs_T(ax[8], marker = '.')
# names[8] = "R_h_tau_HS_0"
# ax[8].plot(T_exp, Rh_exp, c= 'C3', label = "exp results")
# ax[8].set(xscale = 'log', xlim = 1, title = r"LCCO with $\tau{HS}=0$ CVC not T dep")
# pt.make_dotlines(ax[8],sim, HS= False)

# out_tau_hs_0_cvc_not_t.plot_rho_ij_vs_T("xx",ax[3])
# names[3] = "rho_xx_tau_HS_is_0_cvc_T_indep"
# ax[3].set(title = r"$\rho_{xxcuprates}$ with $\tau{HS}=0$ and CVC not T-dependant")

out_tau_hs_dans_beta.plot_RH_vs_T(ax[9], marker = '.')
names[9] = "R_h_tau_HS_dans_beta"
ax[9].plot(T_exp, Rh_exp, c= 'C3', label = "exp results")
ax[9].set(xscale = 'log', xlim = 1, title = r"LCCO with $\tau{HS}$ dans $\beta$ and CVC T-Dep")
pt.make_dotlines(ax[9],sim, HS= False)

out_tau_hs_dans_beta.plot_rho_ij_vs_T("xx",ax[10])
names[10] = "rho_xx_tau_HS_dans_beta"
ax[10].set(title = r"$\rho_{xx}$ with $\tau{HS}$, dans $\beta$ CVC T-dependant")

out_tau_hs_dans_beta.plot_drho_dT_vs_b(ax[7])
names[7] = "drho_dT tau_hs dans beta cvc T dep"
ax[7].set(title = r"$d\rho_{xx}/dT$ with $\tau{HS}$ dans $\beta$ CVC T-dependant")# tau_hs 1 cvc not T dep") dans $\beta$ CVC T-dependant")
out_tau_hs_dans_beta.make_legend(fig[7])
ax[7].set(xlim =0)
# ax[4].plot(out_tau_HS_0.T, out_tau_HS_0.rho_T[0,2,:], c="C0", label = "tau_hs_0")
# ax[4].plot(out_tau_HS_1.T, out_tau_HS_1.rho_T[0,2,:], c="C1", label = "tau_hs_1")
# ax[4].plot(out_tau_HS_T_dep.T, out_tau_HS_T_dep.rho_T[0,2,:], c="C2", label = "tau_hs T dep")

# ax[4].plot(out_tau_HS_0_NO_CVC.T, out_tau_HS_0_NO_CVC.rho_T[0,2,:], c="C3", label = "no_cvc_tau_0")
# ax[4].plot(out_tau_HS_1_NO_CVC.T, out_tau_HS_1_NO_CVC.rho_T[0,2,:], c="C4", label = "tau_hs_1_no_cvc")
# ax[4].plot(out_tau_HS_T_dep_NO_CVC.T, out_tau_HS_T_dep_NO_CVC.rho_T[0,2,:], c="C5", label = "tau_hs T dep no cvc")

# ax[4].plot(out_tau_hs_dans_beta.T, out_tau_hs_dans_beta.rho_T[0,0,:], c="C6", label = "tau_hs dans beta")
# ax[4].plot(out_tau_hs_1_cvc_not_t.T, out_tau_hs_1_cvc_not_t.rho_T[0,2,:], c="C7", label = "tau_hs 1 cvc not T dep")

ax[4].set(xlabel = "$T$ (K)", ylabel = r"$\rho_{xx}$", xlim=1)
ax[4].legend()

out_LCCO_B.plot_rho_ij_vs_B('xx',ax[5])
out_LCCO_B.make_legend(fig[5])

i = 0
for n in names:
    if n!="":
        n+= f"{sim.metal.name}_{sim.tau_choice}_FLEX={sim.FLEX}_SCTIF"
        fig[i].savefig("./SCTIF/results pdf/" + n + ".pdf")
    i+=1

out_LCCO_B.plot_MR(ax[12])
out_LCCO_B.plot_drho_db_vs_b(ax[13])

out_LCCO_B.make_legend(fig[12])
out_LCCO_B.make_legend(fig[13])

T=[20,30,40,50] #,75,110,150,190]
i=0
fg,axss = plt.subplots(1,7, sharex = True, sharey = True)
for Ti in T :
    gamma = get_tau_inv(kx,ky, Ti, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, False)
    cmap = plt.cm.Reds
    axss[i].plot(out_LCCO_B.B, (out_LCCO_B.rho_B[0,:,i]-out_LCCO_B.rho_B[0,0,i])/out_LCCO_B.rho_B[0,0,i], c=out_LCCO_B.color[i], label = f"${Ti}$ K")
    axss[i].set(xlabel = "$B$ (T)", ylabel = r'$\Delta\rho/\rho_0$')
    axss[i].legend()
    pt.plot_phi(ax[14],phi, gamma, label = fr'$T={Ti}$ K', c=out_LCCO_B.color[i])

    i+=1

# ## Calcul de l'effet Hall à 5K et 5T
# dt = 1e-15
# N_tau = 500
# B = 5
# #%% LATTICE
#
# rho_xx, rho_xy, rho_yx, rho_yy = compute_resistivity_JZ(B, 5, 100, dt, dt*N_tau, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, sim.in_AFBZ, sim.FLEX)
# print(f"L'EFFET HALL À {B}T ET 5K EST {rho_xy/B}")
# print(rf"LA $\rho_xx$ À {B}T ET 5K EST {rho_xx}")
# rho_xx, rho_xy, rho_yx, rho_yy = compute_resistivity(B, 5, 100, dt, dt*N_tau, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, sim.in_AFBZ, sim.FLEX)
# print(f"L'EFFET HALL À {B}T ET 5K EST {rho_xy/B}")
# print(rf"LA $\rho_xx$ À {B}T ET 5K EST {rho_xx}")
# plt.show()

# sim.tau_coeff
#%% Résistivity
#
# ij = "xx"
# out_LCCO.plot_rho_ij_vs_B(ij,ax3, marker='.',xlim = (1,20))
# out_LCCO.plot(ax2, marker='.', xlim = (1,20))
# ax2.set_title(fr'$R_H$ vs T, metal = {sim.metal.name}')
# ax3.set_title(fr'$\rho_{ij}$ vs T, metal = {sim.metal.name}')
# out_LCCO.make_legend(fig2)
# out_LCCO.make_legend(fig3)

############
# ARCHIVES #
############
#

# J_f = np.sqrt(Jx**2 + Jy**2)
# vx, vy = get_v0(kx, ky, sim.metal.s_coeff)
# ax0.yaxis.set_major_locator(NullLocator())
# ax0.xaxis.set_major_locator(NullLocator())
# ax1.set_title(r"$v_x(\phi)$ and $J_x(\phi)$")
# pt.plot_phi(ax4, phi, Jx,hole_like, label = f"$J_x$ {sim.metal.name}")
# pt.plot_phi(ax4, phi, vy,hole_like,c='red', label = f"$v_y$ {sim.metal.name}")
# pt.plot_phi(ax4, phi, Jy,hole_like,c='blue', label = f"$J_y$ {sim.metal.name}")
# # pp.set_param(fig1,ax1)
# ax1.set_xlim(0)
# ax4.set_xlim(0)
# pt.set_x_ax(ax1,"phi")
# pt.set_x_ax(ax4,"phi")
# ax1.legend()
# ax4.legend()
# tau_phi = get_tau_inv(kx,ky, 15, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, hole_like)
# lx_phi =Jx/tau_phi #Jx
# ly_phi = Jy/tau_phi
# l = np.sqrt(lx_phi**2 + ly_phi**2)
# # pt.plot_phi(ax4, p + hi, lx_phi, hole_like )
# ax5.plot(phi,l)
# ax5.yaxis.set_major_locator(NullLocator())
# HS= sim.HS
# ax5.vlines([HS, pi/2 - sim.HS,pi/2 + HS, pi - HS , pi + HS, 3*pi/2 - HS, 3*pi/2 + HS, 2*pi - HS], np.min(l), np.max(l), colors='r', label = 'HS') #, transform=ax5.get_xaxis_transform())
# ax5.plot(phi + pi,l, c= "C0")
# ax5.legend()
# ax5.set_title(r"$|l(\phi)| = \tau(\phi)*v(\phi)$")

#
# # sim_hole = Simulation(LCCOH)

# LSCO = Metal('LSCO',p='0.32')
# LSCO2 = Metal('LSCO',p='0.12')
# LSCO3 = Metal('LSCO',p='0.20')
# LSCO4 = Metal('LSCO',p='0.24')
#
# s_2 = sim_hole.metal.s_coeff
# kx2,ky2 = get_kF(s_2,100, False)
# hole_like2 = (0. not in kx2)
# print("COEFF 2", sim_hole.metal.s_coeff)
# phi2 = get_phi(kx2,ky2,s_2, hole_like2)
# phi2 = phi2 + np.pi*(hole_like2)
# vx2, vy2 = get_v0(kx2, ky2, s_2)
# Jx2, Jy2 = get_v(kx2, ky2, s_2, FLEX=sim.FLEX)
# sim_hole.metal.plot_FS(ax3)
# pt.plot_phi(ax4, phi2, vx2,hole_like2,c='green', label = f"$v_x$ {"hole"*hole_like2 + "electron"*(1-(hole_like2))}")
# pt.plot_phi(ax4, phi2, Jx2,hole_like2,c='purple', label = f"$J_x$ {"hole"*hole_like2 + "electron"*(1-(hole_like2))}")
# ax4.legend()

# list_p = [p[2:] for p in list(metals['LCCO'].keys())[3:]]
# list_LCCO = [Metal('LCCO',p=p) for p in list_p]
# ax1_t = ax0
# sim.metal.plot_FS(ax1)

# ax1.set_xscale('log')
# out_LCCO.plot_MR(ax1, ylim= (0,0.08), marker='.')
# alpha = get_epsilon(kx,ky,sim.metal.s_coeff)#, sim.HS)
# pt.set_x_ax(ax0,"phi")

# pt.plot_phi(ax1,phi, alpha,hole_like )



# out_LCCO.make_legend(fig0)
# out_LCCO.make_legend(fig1)

plt.show()
