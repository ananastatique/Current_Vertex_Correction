from numpy._core.umath import sign
from Code_SCTIF.BTE_routines_SCTIF import *
from Code_SCTIF.metals_SCTIF import *
from Code_SCTIF.plot_results_SCTIF import *
import Code_SCTIF.plot_tau_phi_SCTIF as pt
from simul_SCTIF import Simulation
import matplotlib.pyplot  as plt
from matplotlib.ticker import NullLocator
import matplotlib as mpl

plt.style.use("pacoty.mplstyle")

LCCOH = Metal('LCCOHole', p='0.150')
LCCO = Metal('LCCO', p='0.15')
pi = np.pi
sim = Simulation()

N=1
figs = []
for i in range(N) :
    figs.append(plt.subplots(figsize=(6, 4.5), constrained_layout=True))
ax = [f[1] for f in figs]
fig = [f[0] for f in figs]

kx,ky = get_kF(sim.metal.s_coeff,200,False)
Jx, Jy = get_v(kx, ky, sim.metal.s_coeff, FLEX=sim.FLEX)
vx, vy = get_v0(kx, ky, sim.metal.s_coeff)

hole_like = False #(0. not in kx)
phi = get_phi(kx,ky,sim.metal.s_coeff, hole_like)
phi = phi + np.pi*hole_like
phi_min = -1#sim.HS #0#sim.HS
phi_max = 15#np.pi/2 - sim.HS #np.pi #sim.HS
mask_phi = phi < phi_max
mask_phi_min = phi[mask_phi] > phi_min

gamma = get_tau_inv(kx,ky,40, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, hole_like = False)
norm = mpl.colors.Normalize(vmin=-np.max(np.abs(Jx*ky - Jy*kx)), vmax=np.max(np.abs(Jx*ky - Jy*kx)))
col = np.abs(Jx*ky- kx*Jy)
sim.metal.plot_FS(ax[0], isoE = False, phi_min=phi_min, phi_max=phi_min)
ax[0].plot(kx+(gamma - np.min(gamma))*np.cos(phi)*1e-5*2, ky+ (gamma - np.min(gamma))*np.sin(phi)*1e-5*2, c="C1")
ax[0].plot(-kx-(gamma - np.min(gamma))*np.cos(phi)*1e-5*2, -ky- (gamma - np.min(gamma))*np.sin(phi)*1e-5*2, c="C1")
# ax[0].quiver(kx,ky ,Jx/gamma,Jy/gamma, col, cmap = mpl.colormaps["coolwarm"], width = 0.002, scale =  0.00000005, norm = norm) #scale = 0.00000005 np.zeros_like(phi) phi, cmap=cmap
# ax[0].quiver(kx,-ky ,Jx/gamma, -Jy/gamma, col, cmap = mpl.colormaps["coolwarm"], width = 0.002, scale =  0.00000005, norm = norm) #scale = 0.00000005 np.zeros_like(phi)a
# ax[0].plot([0,15*1e9], [0, sim.HS*15*1e9])
# for T in [20,30,40,50]:
# #     # pt.plot_tau_phi(ax0, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, T)
#     # gamma = get_tau_inv(kx,ky,T, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, hole_like = False)
#     # lx = vx/gamma
#     # ly = vy/gamma
#     # dly = np.gradient(ly, phi)
#     # pt.plot_phi(ax[4], phi, lx*dly, label = f"$T=${T}K")
#     # ax[4].scatter(phi, lx*dly, c = (lx*dly)/np.abs(lx*dly), cmap = mpl.colormaps["coolwarm"], label = f"$T=${T}K")
#     # ax[4].set_ylabel(r"$dl_y/d\phi*l_x$")
#     # pt.plot_phi(ax[5], phi, gamma, hole_like, label = f"$l$ {sim.metal.name}")
#     pass

# T = np.linspace(1,200)
# B = np.linspace(0,15,20)
# cmap = plt.cm.Reds
# T = np.linspace(0,200,50)
# color = [cmap(i) for i in np.linspace(0, cmap.N, len(T)+2, dtype=np.int32)][2:]
# i=0
# for Ti in np.linspace(0,200,50): #[20,30,40,50]:#
#     # print(Ti)
#     # RHO = [compute_resistivity_JZ(Bi, Ti, 200, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, sim.in_AFBZ, sim.FLEX) for Bi in B]
#     # rho_xx = np.array([R[0] for R in RHO])
#     # rho_xy = np.array([R[1] for R in RHO])
#     # ax[6].plot(B, (rho_xx - rho_xx[0])/rho_xx[0], marker = '.', label = f'{Ti}K')
#     gamma = get_tau_inv(kx,ky,Ti, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, hole_like = False)
#     pt.plot_phi(ax[5],phi,gamma, c = color[i])# label = f'{Ti}K',
#     ax[5].set_ylabel(r"$\Gamma$ ($ps^{-1}$)")
#     i+=1
# i=0
# T = np.linspace(0,50,10)
# color = [cmap(i) for i in np.linspace(0, cmap.N, len(T)+2, dtype=np.int32)][2:]
# gamma = get_tau_inv(kx,ky,T[0], sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, hole_like = False)
# # pt.plot_phi(ax[1],phi,gamma, c = color[i], label="0K")# label = f'{Ti}K',
# # pt.plot_phi(ax[3].twinx(),phi,vx, c = color[8])# label = f'{Ti}K',
# # pt.plot_phi(ax[3],phi,vx, c = "C0")# label = f'{Ti}K',
# #
# for Ti in np.linspace(0,50,10): #[20,30,40,50]:#
#     # print(Ti)
#     # RHO = [compute_resistivity_JZ(Bi, Ti, 200, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, sim.in_AFBZ, sim.FLEX) for Bi in B]
#     # rho_xx = np.array([R[0] for R in RHO])
#     # rho_xy = np.array([R[1] for R in RHO])
#     # ax[6].plot(B, (rho_xx - rho_xx[0])/rho_xx[0], marker = '.', label = f'{Ti}K')
#     # gamma = get_tau_inv(kx,ky,Ti, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, hole_like = False)
#     # Jx, Jy = get_v(kx, ky, sim.metal.s_coeff, Ti, sim.FLEX, mod = True)
#     # Jx_not_modulo, Jy_not_mode = get_v(kx, ky, sim.metal.s_coeff, Ti, sim.FLEX, mod = False)

#     # pt.plot_phi(ax[3],phi,Jx, c = color[i])# label = f'{Ti}K',

#     # ax[6].set_ylabel('eps CVC turning on')
#     # pt.plot_phi(ax[1],phi,Jx, c = color[i])# label = f'{Ti}K',

#     # ax[1].set_ylabel(r"$\Gamma$ ($ps^{-1}$)")
#     i+=1

# Jx, Jy = get_v(kx, ky, sim.metal.s_coeff, Ti, sim.FLEX)
# Jx_not_modulo, Jy_not_mode = get_v_not_modulo(kx, ky, sim.metal.s_coeff, Ti, sim.FLEX)

# pt.plot_phi(ax[6],phi,get_epsilon(kx,ky, sim.metal.s_coeff, Ti), c = color[i-1])# label = f'{Ti}K',
# pt.plot_phi(ax[6],phi,get_epsilon_not_modulo(kx,ky, sim.metal.s_coeff, sim.HS, Ti), c = "C0")# label = f'{Ti}K',

# gamma = get_tau_inv(kx,ky,T[-1], sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, hole_like = False)
# pt.plot_phi(ax[1],phi,Jy, c = color[i-1], label = f'Jx')
# pt.plot_phi(ax[1],phi,Jy_not_mode, c = "C0", label = f'Jx not mod')

# ax[1].set_ylabel(r"Jx")
# ax[1].legend()

# sigmaxy_B = [compute_sigma_ij(B[1],Ti, "xy", 200, 1e-15, 10e-12, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, sim.in_AFBZ, sim.FLEX, phi_min = -1, phi_max = 15) for Ti in T]
# sigmaxx_B = [compute_sigma_ij(B[1],Ti, "xx", 200, 1e-15, 10e-12, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, sim.in_AFBZ, sim.FLEX, phi_min = -1, phi_max = 15) for Ti in T]

# sigmaxy_B = np.array(sigmaxy_B)
# sigmaxx_B = np.array(sigmaxx_B)

# det_sigma = sigmaxx_B**2 + sigmaxy_B**2

# pt.make_dotlines(ax[1], sim)#, zero =False)
# ax[1].legend()

# ax[2].plot(T, -2*sigmaxy_B/(sigmaxx_B**2 +sigmaxy_B**2),marker = '.')
# ax[2].set_ylabel( r"$dRh/dT$")
# ax[2].set_xlabel(r"$B$ (T)")

# ax[3].plot(T, np.gradient(sigmaxy_B,T) + 2*sigmaxy_B/sigmaxx_B*np.gradient(sigmaxx_B,T))
# ax[3].twinx().plot(T, np.gradient(sigmaxy_B,T), color = 'r', label= "approx")# + 2*sigmaxy_B/sigmaxx_B*np.gradient(sigmaxx_B,T))
# ax[3].set(xlabel = "T", ylabel= "signe de dsigma/dT")
# ax[3].legend()
#
#
# ax[4].plot(T, 2*sigmaxy_B/sigmaxx_B*np.gradient(sigmaxx_B,T))
# ax[4].set(xlabel = 'T', ylabel = 'dsigma/dT*sxy/sxx')



# ax[0].set_ylabel(r"$\mathrm{d}\sigma_{xy}/\mathrm{d}\phi$")
# ax[0].axvline(x=sim.HS, ls="--", lw = 1)
# ax[0].axhline(y=0, ls="--",c="r", lw = 1)
# ax[0].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
# ax[4].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
# ax[4].axvline(x= sim.HS, ls="--", lw = 1)
# ax[4].axvline(x= np.pi/2 + sim.HS, ls="--", lw = 1)
# ax[4].axvline(x= np.pi - sim.HS, ls="--", lw = 1)
# ax[2].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
# ax[2].axvline(x= sim.HS, ls="--", lw = 1)
# ax[2].axvline(x= np.pi/2 + sim.HS, ls="--", lw = 1)
# ax[2].axvline(x= np.pi - sim.HS, ls="--", lw = 1)
# ax[0].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
# ax[0].axvline(x=np.pi/2 + sim.HS, ls="--", lw = 1)
# ax[0].axvline(x=np.pi - sim.HS, ls="--", lw = 1)
    # pt.plot_phi(ax[5].twinx(), phi, Jy/gamma,hole_like, label = f"$l$ {sim.metal.name}")
#     pt.plot_phi(ax[1], phi, np.sqrt(Jx**2+Jy**2)/gamma,hole_like, label = f"$l$ {sim.metal.name}")
    # pt.plot_phi(ax[3], phi, np.sqrt(Jx**2+Jy**2)/gamma, False, label = "$|l|$")


# ax[4].plot(kxt, kyt)
# ax[4].scatter(kxt, kyt,c=t, norm = norm, cmap = mpl.colormaps["viridis"])
# ax[4].plot(kxt_F, kyt_F)
# ax[4].scatter(kxt_F, kyt_F, c=t, norm = norm, cmap = mpl.colormaps["plasma"])
# ax[4].axis("equal")

T = np.linspace(1,200, 50)
B = np.array([10,2,3,5,10,15])

Tsw = np.linspace(1,200,30)


# dsigma_dphi = compute_dsigma_dphi(B[0], 150, "xy", 200, 1e-15, 10e-12, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, sim.in_AFBZ, sim.FLEX, phi_min=phi_min, phi_max=phi_max)
# pt.plot_phi(ax[0], phi, dsigma_dphi)
# for Bi in [0,1,2,3,5,10,15]:
#     sigma_xx_B = [compute_sigma_ij(Ti, Bi, "xx", 200, 1e-15, 10e-12, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, sim.in_AFBZ, sim.FLEX, phi_min=phi_min, phi_max=phi_max) for Ti in Tsw]
#     sigma_xx_B = np.array(sigma_xx_B)
#     ax[6].plot(Tsw,sigma_xx_B, marker = '.', label = f"$B={Bi}$T") #np.gradient(sigma_xx_B, Tsw), marker = '.')
#     ax[6].set(xlabel = "$T$ (K)", ylabel = r"$\sigma{xx}$")
#     ax[6].legend()
# ax[0].set_ylabel(r"$\mathrm{d}\sigma_{xy}/\mathrm{d}\phi$")
# ax[0].axvline(x=sim.HS, ls="--", lw = 1)
# ax[0].axhline(y=0, ls="--",c="r", lw = 1)
# ax[0].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
# ax[4].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
# ax[4].axvline(x= sim.HS, ls="--", lw = 1)
# ax[4].axvline(x= np.pi/2 + sim.HS, ls="--", lw = 1)
# ax[4].axvline(x= np.pi - sim.HS, ls="--", lw = 1)
# ax[2].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
# ax[2].axvline(x= sim.HS, ls="--", lw = 1)
# ax[2].axvline(x= np.pi/2 + sim.HS, ls="--", lw = 1)
# ax[2].axvline(x= np.pi - sim.HS, ls="--", lw = 1)
# ax[0].axvline(x=np.pi/2 - sim.HS, ls="--", lw = 1)
# ax[0].axvline(x=np.pi/2 + sim.HS, ls="--", lw = 1)
# ax[0].axvline(x=np.pi - sim.HS, ls="--", lw = 1)

# sigma_xy =   sim.compute_T_sweep_conductivity("xy")
# cmap = plt.cm.Blues
# color = [cmap(i) for i in np.linspace(0, cmap.N, B.size+2, dtype=np.int32)][2:]

# for i in range(B.size):
#     print( sigma_xy[i,:].shape)
#     ax[0].plot(T, sigma_xy[i,:], label = f"$B= {B[i]}$ T", c= color[i])

# ax[0].set_xlabel("$T$ (K)")
# ax[0].set_ylabel(r"$\sigma_{xy}$")
# ax[0].legend()

# dJy = np.gradient(Jy/gamma, phi)
# norm = mpl.colors.CenteredNorm()#Normalize(vmin = -1, vmax = 1) # = np.min(Jx*dJy), vmax = -np.min(Jx*dJy))
# norm2 = mpl.colors.Normalize(vmin = -np.abs(np.min(Jx/gamma*dJy)), vmax = -np.abs(np.min(Jx/gamma*dJy))) # = np.min(Jx*dJy), vmax = -np.min(Jx*dJy))
# color = Jx/gamma*dJy/np.abs(Jx/gamma*dJy)
# pt.plot_phi(ax[2], phi, Jx, False, label = "$J_x$")
# pt.plot_phi(ax[2], phi, vx, False, label = "$v_x$")
# ax[2].set_ylabel(r"$v_x$")
# ax[2].axhline(y=0, ls="--", lw=1, c="r")
# ax[2].scatter(phi, Jx, c =  color, norm= norm, cmap = mpl.colormaps["coolwarm"])

# lx = Jy/gamma

# pt.plot_phi(ax[3], phi, lx, False, label = "$J_x$")
# pt.plot_phi(ax[3].twinx(), phi, dJy, False, label = "$v_x$")
# ax[3].set_ylabel(r"$l_y$")
# ax[3].axhline(y=0, ls="--", lw=1, c="r")
# ax[3].scatter(phi, lx, c =  color, norm= norm, cmap = mpl.colormaps["coolwarm"])


# pt.plot_phi(ax[4], phi, Jx/gamma*dJy, color, False, label = "$J_x$")
# # pt.plot_phi(ax[4], phi, vx, False, label = "$v_x$")
# ax[4].scatter(phi, Jx/gamma*dJy,c = color, norm= norm, cmap = mpl.colormaps["coolwarm"])
# ax[4].set_ylabel(r"$l_x*dl_y/\phi$")
# ax[4].axhline(y=0, ls="--", lw=1, c="r")


# a = sim.metal.a
# ax[1].plot([0, 2*np.pi/a*np.cos(sim.HS)],[0, 2*np.pi/a*np.sin(sim.HS)], "--", lw = 1, c="r")
# ax[1].plot([0, 2*np.pi/a*np.cos(pi/2 - sim.HS)],[0, 2*np.pi/a*np.sin(np.pi/2 - sim.HS)], "--", lw = 1, c="r")

# ax[0].scatter(phi, dsigma_dphi, c =  color, norm= norm, cmap = mpl.colormaps["coolwarm"])

# cmap = mpl.colormaps["twilight"]
# cbar2 = plt.colorbar(mpl.cm.ScalarMappable(norm= norm,cmap = mpl.colormaps["coolwarm"]),ax=ax[2], orientation='vertical', label = r"$J_x\mathrm{d}J_y/\mathrm{d}\phi$")
# cbar0 = plt.colorbar(mpl.cm.ScalarMappable(norm= norm,cmap = mpl.colormaps["coolwarm"]),ax=ax[0], orientation='vertical', label = r"$J_x\mathrm{d}J_y/\mathrm{d}\phi$")
# cbar3 = plt.colorbar(mpl.cm.ScalarMappable(norm= norm,cmap = mpl.colormaps["coolwarm"]),ax=ax[3], orientation='vertical', label = r"$J_x\mathrm{d}J_y/\mathrm{d}\phi$")
# sign_sigma = dsigma_dphi


# # ax[5].plot(kxt, kyt, marker = '.')

# # sim.metal.plot_FS(ax0)
# out_LCCO_T.plot_RH_vs_T(ax[6])
# out_LCCO_B.plot_MR(ax[8])
# # out_LCCO.plot_MR(ax2)
# out_LCCO_T.plot_rho_ij_vs_T("xy",ax[7])
# out_LCCO_B.plot_drho_db_vs_b(ax[9])

# out_LCCO_T.make_legend(fig[7])
# out_LCCO_T.make_legend(fig[6])

# out_LCCO_B.make_legend(fig[8])
# out_LCCO_B.make_legend(fig[9])



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
# phi  = phi[mask_phi][mask_phi_min]
# Jx = Jx[mask_phi][mask_phi_min]
# Jy = Jy[mask_phi][mask_phi_min]
# kx = kx [mask_phi][mask_phi_min]
# ky = ky[mask_phi][mask_phi_min]
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
# B = 60
# i0 = 25# sim.metal.plot_FS(ax0)

# norm = mpl.colors.Normalize (vmin = 0, vmax = 1000)
# t = np.linspace(0, 1000, 1000)
# kxt, kyt = evolve_kF(kx[i0], ky[i0], B, 1e-15, 1000, sim.metal.s_coeff, False, True)
# kxt_F, kyt_F = evolve_kF(kx[i0], ky[i0], B, 1e-15, 1000, sim.metal.s_coeff, False, False)

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
