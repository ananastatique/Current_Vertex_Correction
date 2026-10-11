from Code_SCTIF.BTE_routines_SCTIF import *
from Code_SCTIF.metals_SCTIF import *
from Code_SCTIF.plot_results_SCTIF import *
import Code_SCTIF.plot_tau_phi_SCTIF as pt
import matplotlib.pyplot  as plt
from matplotlib.ticker import NullLocator
from simul_SCTIF import Simulation
import matplotlib as mpl
import pandas as pd
import os

plt.style.use("pacoty.mplstyle")

#SIMULATION
sim = Simulation()

## FILES
out_LCCO_T = Results(f'./SCTIF/npz/{sim.metal.name}_{sim.tau_choice}_{int(sim.tau_coeff[-2])*"spin" + (1-int(sim.tau_coeff[-2]))*"charge"}_T=1_to_100_at_B=[ 1  2  5 10]_FLEX={sim.FLEX}.npz')
out_LCCO_B = Results(f'./SCTIF/npz/{sim.metal.name}_{sim.tau_choice}_{int(sim.tau_coeff[-2])*"spin" + (1-int(sim.tau_coeff[-2]))*"charge"}_B=0_to_35_at_T=[20 30 40 50]_FLEX={sim.FLEX}.npz')
# out_LCCO_B_15 = Results("./SCTIF/npz/LCCO_0.150_hot spots_spin_T=1_to_100_at_B=[ 1  2  3  5 10 15]_FLEX=True.npz")
## INITIALIZING FIGURES and labels
N=13
figs = []
names = []
for i in range(N) :
    figs.append(plt.subplots())
    names.append("")
ax = [f[1] for f in figs]
fig = [f[0] for f in figs]

label_B = r"$B$ (T)"
label_T = r"$T$ (K)"
met_name_wo_doping = sim.metal.name.split("_")[0]
title_set_up = rf"$-$ {met_name_wo_doping} {sim.metal.p} \texttt{{FLEX={sim.FLEX}}}"
s_or_c = sim.tau_coeff[-2]
#COMPUTING INITIAL VALUES OF VARIABLES

kx,ky = get_kF(sim.metal.s_coeff,200,False)
Jx, Jy = get_v(kx, ky, sim.metal.s_coeff, FLEX=sim.FLEX)#, mode=s_or_c)#), mode=sim.tau_coeff[-2])
vx, vy = get_v0(kx, ky, sim.metal.s_coeff)

hole_like = False #(0. not in kx)
phi = get_phi(kx,ky,sim.metal.s_coeff, hole_like)
phi = phi + np.pi*hole_like
phi_min = -1
phi_max = 15
#COLORS
cmap = plt.cm.Reds
T=[0,5,10,20,30,40,50,60,80,100]
reds = [cmap(i) for i in np.linspace(0, cmap.N, len(T)+2, dtype=np.int32)][2:]

#EXEPERIMENTAL DATA

og_data_rh_14 = pd.read_csv('r_h_exp_data_x_0.14.csv')
og_data_rh_15 = pd.read_csv('r_h_exp_data_x_0.15.csv')
og_data_rh_16 = pd.read_csv('r_h_exp_data_x_0.16.csv')
og_data_rh_18 = pd.read_csv('r_h_exp_data_x_0.18.csv')
T_exp_14 = og_data_rh_14["Temperature"]
Rh_exp_14 = og_data_rh_14["Rh"]
T_exp_15 = og_data_rh_15["Temperature"]
Rh_exp_15 = og_data_rh_15["Rh"]
T_exp_16 = og_data_rh_16["Temperature"]
Rh_exp_16 = og_data_rh_16["Rh"]
T_exp_18 = og_data_rh_18["Temperature"]
Rh_exp_18 = og_data_rh_18["Rh"]

ls = "--",
og_data_rho_xx = pd.read_csv('rho_T_0.159_exp.csv')
T_rho_exp = og_data_rho_xx["Temperature"]
Rho_exp = og_data_rho_xx["Rho"]
drho_compute = np.gradient(Rho_exp, T_rho_exp)

og_data_drho_xx = pd.read_csv('dho_dT_exp_0.15.csv')
T_drho_exp = og_data_drho_xx["Temperature"]
dRho_exp = og_data_drho_xx["drho"]

og_data_linear_fit_exp_drho_xx = pd.read_csv('dashed_drho_dxx_T_exp.csv')
T_drho_fit_exp = og_data_linear_fit_exp_drho_xx["Temperature"]
dRho_fit_exp = og_data_linear_fit_exp_drho_xx["drho"]

#PLOTTING MR
out_LCCO_B.plot_MR(ax[0])
names[0] = "MR"
out_LCCO_B.make_legend(fig[0])
ax[0].set(ylim = (0), title = "MR "+title_set_up)

out_LCCO_B.plot_drho_db_vs_b(ax[1])
T = [20, 30, 40, 50]
names[1] = "dMR_dB"
out_LCCO_B.make_legend(fig[1])
ax[1].set(ylim = (0), title = r"$d\rho_{xx}/dB$" + title_set_up)

#plot Hall Number and exp data
out_LCCO_T.plot_RH_vs_T(ax[2])
# out_LCCO_B_15.plot_RH_vs_T(ax[2])
names[2] = "RH"
pt.make_dotlines(ax[2], sim, HS = False)
ax[2].plot(T_exp_14 , Rh_exp_14, c = "C3", label = "exp results 0.14") # -Rh_exp_14[0]
# ax[2].plot(T_exp_15, Rh_exp_15, c = "C3", label = "exp results 0.15") # - Rh_exp_15[0]
# ax[2].plot(T_exp_16, Rh_exp_16, c = "C0", label = "exp results 0.16") # - Rh_exp_16[0]
# ax[2].plot(T_exp_18, Rh_exp_18, c = "C2", label = "exp results 0.18") # - Rh_exp_18[0]
ax[2].set(xscale = "log", title = r"$R_H$" + title_set_up)
out_LCCO_T.make_legend(fig[2])

#Plot velocity along x and y
pt.plot_phi(ax[3], phi, vx, c=reds[-1], label="$v_x$")
ax[3].set(ylabel = r"$v_x$ ($\mathrm{ms^{-1}}$)", title=rf"$v_x$ {title_set_up}")
names[3] = "v_x"
pt.make_dotlines(ax[3], sim)

pt.plot_phi(ax[4], phi, vy, c=reds[-1],  label="$v_y$")
ax[4].set(ylabel = r"$v_y$ ($\mathrm{ms^{-1}}$)", title=rf"$v_y$ {title_set_up}")
names[4] = "v_y"
pt.make_dotlines(ax[4], sim)

#PLOT Jx and Jy CVC
pt.plot_phi(ax[3], phi, Jx, c="C0", label="$v^*_x$")
ax[3].set(ylabel = r"$v_x$ ($\mathrm{ms^{-1}}$)", title=rf"$v_x$ {title_set_up}")
names[3] = "J_x"

pt.plot_phi(ax[4], phi, Jy, c="C0", label = "$v^*_y$")
ax[4].set(ylabel = r"$v_y$ ($\mathrm{ms^{-1}}$)", title=rf"$v^*_y$ {title_set_up}")
names[4] = "J_y"

#Plot GAMMA, JX and EPSILON varying with temperature :
T=[0,5,10,20,30,40,50,60,80,100]

i=0
for Ti in T:
    gamma = get_tau_inv(kx,ky,Ti, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, hole_like = False)
    Jx, Jy = get_v(kx, ky, sim.metal.s_coeff, Ti, FLEX=sim.FLEX)#, mode = s_or_c)
    eps_T = get_epsilon(kx, ky, sim.metal.s_coeff, Ti)#, mode=s_or_c)
    pt.plot_phi(ax[5], phi, gamma, label = f"$T=${Ti}K", c= reds[i])
    pt.plot_phi(ax[10], phi, Jx, label = f"$T=${Ti}K", c= reds[i])
    pt.plot_phi(ax[11], phi, eps_T, label = f"$T=${Ti}K", c= reds[i])

    i+=1

ax[5].set(ylabel = r"$\Gamma \ (\mathrm{ps^{-1}})$", title = rf"$\Gamma$ {title_set_up}")
names[5] =  "Gamma"
ax[10].set(ylabel = r"$J_x \ (\mathrm{ms^{-1}})$", title = rf"$J_x$ {title_set_up}", ylim = (-500000,500000))
names[10] =  "Jx_cvc_turning_on"
ax[11].set(ylabel = r"$\varepsilon(k)$", title = rf"$\varepsilon$ {title_set_up}", ylim = 0)
names[11] =  "epsilon_T_cvc"
# ax[5].legend()
pt.make_dotlines(ax[5], sim)
pt.make_dotlines(ax[10], sim)
pt.make_dotlines(ax[11], sim)
pt.make_legend(ax[5], fig[5])
pt.make_legend(ax[10], fig[10])
pt.make_legend(ax[11], fig[11])

#PLOT FERMI SURFACE
sim.metal.plot_FS(ax[6], isoE = False, phi_min=-1, phi_max=-1)
norm = mpl.colors.Normalize(vmin=-np.max(np.abs(Jx*ky - Jy*kx)), vmax=np.max(np.abs(Jx*ky - Jy*kx)))
ax[6].quiver(kx, ky, Jx, Jy, np.abs(Jx*ky - Jy*kx), cmap = mpl.colormaps["coolwarm"], norm = norm, width = 0.004, scale =  5000000) # np.abs(Jx*ky - Jy*kx), cmap = mpl.colormaps["coolwarm"], #scale = 0.00000005 np.zeros_like(phi) phi, cmap=cmap
ax[6].quiver(kx, -ky ,Jx, -Jy,np.abs(Jx*ky - Jy*kx), cmap = mpl.colormaps["coolwarm"], norm = norm, width = 0.004, scale =  5000000) #  np.abs(Jx*ky - Jy*kx), cmap = mpl.colormaps["coolwarm"] #scale = 0.00000005 np.zeros_like(phi)
ax[6].legend()
names[6] = "FS"

#plot \epsilon(k)
epsilon = get_epsilon(kx, ky, sim.metal.s_coeff)#, mode= s_or_c)
pt.plot_phi(ax[7], phi, epsilon)
names[7] = "epsilon_cvc_vs_phi"
ax[7].set(ylabel = r"$\varepsilon (\phi)$ ([1])", title = r"$\varepsilon (\phi)$"+ title_set_up )
pt.make_dotlines(ax[7],sim)

#plot rho xx
out_LCCO_T.plot_rho_ij_vs_T("xx",ax[8])
names[8] = "rho_xx_vs_T"
ax[8].set(title = r"$\rho_{xx}$ " + title_set_up, xlim = (0,100), ylim = 0)
ax[8].plot(T_rho_exp, Rho_exp, c = "C3", label = "exp results")
out_LCCO_T.make_legend(fig[8])

#plot drho_xx/dT
out_LCCO_T.plot_drho_dT_vs_b(ax[9])
names[9] = "drho_xxdT_vs_T"
ax[9].set(title = r"$d\rho_{xx}/dT$ " + title_set_up, xlim = (0,100), ylim = 0)
ax[9].plot(T_drho_exp, dRho_exp, c = "C3", label = "exp results")
# ax[9].plot(T_drho_fit_exp, dRho_fit_exp, c = "C4", ls = "--", label = "exp results linear fit")
# ax[9].plot(T_rho_exp, drho_compute, c = "C1", ls = "--", label = "exp results co")
out_LCCO_T.make_legend(fig[9])

# Plot k(T) for the Turning on of the CVC

T_k = np.linspace(1,100,100)
k = get_CVC_T_dep(T_k)
ax[12].plot(T_k,k)
ax[12].set(xlabel = label_T, ylabel = r"$\alpha (T)$ ([1])", title = r"$\alpha(T)$"+ title_set_up )
pt.make_dotlines(ax[12],sim,HS=False)
names[12] = "alpha_T_turning_on_CVC"

#plot MR T
#
T=[20,30,40,50]
i=0
fg,axss = plt.subplots(1,7, sharex = True, sharey = True)
for Ti in T :
    gamma = get_tau_inv(kx,ky, Ti, sim.metal.s_coeff, sim.tau_choice, sim.tau_coeff, False)
    cmap = plt.cm.Reds
    axss[i].plot(out_LCCO_B.B, (out_LCCO_B.rho_B[0,:,i]-out_LCCO_B.rho_B[0,0,i])/out_LCCO_B.rho_B[0,0,i], c=out_LCCO_B.color[i], label = f"${Ti}$ K")
    axss[i].set(xlabel = "$B$ (T)", ylabel = r'$\Delta\rho/\rho_0$')
    axss[i].legend()
    i+=1


# Save Figures

os.makedirs("./SCTIF/results propres/png/", exist_ok=True)
os.makedirs("./SCTIF/results propres/svg/", exist_ok=True)

n = "MR_vs_T_diff_axes"
fg.savefig("./SCTIF/results propres/png/" + n + ".png", transparent=False)
fg.savefig("./SCTIF/results propres/" + n + ".pdf", transparent=False)
fg.savefig("./SCTIF/results propres/svg/" + n + ".svg", transparent=False)

s_or_c = int(sim.tau_coeff[-2])
mode = s_or_c*"spin" + (1-s_or_c)*"charge"
i = 0
for n in names:
    if n!="":
        n = f"{sim.metal.name}_FLEX={sim.FLEX}_SCTIF_" + n
        fig[i].savefig("./SCTIF/results propres/" + mode + "_" + n + ".pdf")
        fig[i].savefig("./SCTIF/results propres/png/"+ mode + "_"+ n + ".png", transparent=False)
        fig[i].savefig("./SCTIF/results propres/svg/"+ mode + "_"+ n + ".svg", transparent=False)
    i+=1

plt.show()
