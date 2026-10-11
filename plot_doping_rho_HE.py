from Code_SCTIF.BTE_routines_SCTIF import *
from Code_SCTIF.metals_SCTIF import Metal
from Code_SCTIF.plot_results_SCTIF import *
import Code_SCTIF.plot_tau_phi_SCTIF as pt
import matplotlib.pyplot  as plt
from matplotlib.ticker import NullLocator
from simul_SCTIF import Simulation
import matplotlib as mpl
import pandas as pd
import os
# *.pdf
# *.npz
# *.json
# *test*
# *.pyc

plt.style.use("pacoty.mplstyle")

#SIMULATION
sim = Simulation()
list_dopings = ["0.14","0.15","0.16","0.18"]
cuprate = "LCCO"
list_met = []
sims = []
## FILES
out_LCCO_T = []
out_LCCO_B = []
for i, x in enumerate(list_dopings):
    metal = Metal(cuprate, p=x)
    list_met.append(metal)
    sims.append(Simulation(metal = metal))
    out_LCCO_T.append(Results(f'./SCTIF/npz/{sims[i].metal.name}_{sims[i].tau_choice}_{int(sims[i].tau_coeff[-2])*"spin" + (1-int(sims[i].tau_coeff[-2]))*"charge"}_T=1_to_100_at_B=[ 1  2  5 10]_FLEX={sims[i].FLEX}.npz'))
    out_LCCO_B.append(Results(f'./SCTIF/npz/{sims[i].metal.name}_{sims[i].tau_choice}_{int(sims[i].tau_coeff[-2])*"spin" + (1-int(sims[i].tau_coeff[-2]))*"charge"}_B=0_to_35_at_T=[20 30 40 50]_FLEX={sims[i].FLEX}.npz'))
# out_LCCO_B_15 = Results("./SCTIF/npz/LCCO_0.150_hot spots_spin_T=1_to_100_at_B=[ 1  2  3  5 10 15]_FLEX=True.npz")
## INITIALIZING FIGURES and labels
N=3
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
T=[20, 30, 40, 50]
B = [1, 2 , 5, 10]
reds = [cmap(i) for i in np.linspace(0, cmap.N, len(T)+2, dtype=np.int32)][2:]


fg,axss = plt.subplots(1,4, figsize=(12, 6), sharex = True, sharey = True)

#PLOT HE FOR EACH DOPING
for i, sim in enumerate(sims):
    j = 2
    for (k, Ti) in enumerate(T) :
        if k ==0 :
            axss[k].plot(out_LCCO_B[i].B, (out_LCCO_B[i].rho_B[0,:,k]-out_LCCO_B[i].rho_B[0,0,k])/out_LCCO_B[i].rho_B[0,0,k], c=out_LCCO_B[i].color[i], label = f"$x = {float(list_dopings[i])}$")
        else:
            axss[k].plot(out_LCCO_B[i].B, (out_LCCO_B[i].rho_B[0,:,k]-out_LCCO_B[i].rho_B[0,0,k])/out_LCCO_B[i].rho_B[0,0,k], c=out_LCCO_B[i].color[i])
        axss[k].set(xlabel = "$B$ (T)", ylabel = r'$\Delta\rho/\rho_0$', title = f"$T = {Ti}K$", xlim = (0, 35), ylim = 0)
    ax[0].plot(out_LCCO_T[i].T, out_LCCO_T[i].rho_T[1,j,:]*100/out_LCCO_T[i].B[j], c=out_LCCO_T[i].color[i], label = f"$x = {float(list_dopings[i])}$")
    ax[1].plot(out_LCCO_T[i].T, out_LCCO_T[i].rho_T[0,0,:], c= out_LCCO_T[i].color[i], label = f"$x = {float(list_dopings[i])}$")
    drho_dT = np.gradient(out_LCCO_T[i].rho_T[0,j,:], out_LCCO_T[i].T)
    ax[2].plot(out_LCCO_T[i].T, drho_dT, c= out_LCCO_T[i].color[i], label = f"$x = {float(list_dopings[i])}$")



ax[0].set(xlabel=label_T, ylabel=r"$R_H (10^{-10} \Omega m/T)$", title=f"$R_H$ vs $T$ for $x = {[float(x) for x in list_dopings]}$", xscale='log')
ax[1].set(xlabel=label_T, ylabel=r"$\rho_{xx}$", title=rf"$\rho_{{xx}}$ vs $T$ for $B=0$ and $ x = {[float(x) for x in list_dopings]}$", ylim= 0, xlim = 0)
ax[2].set(xlabel=label_T, ylabel=r"$\partial\rho/\partial T$", title=rf"$\partial\rho/\partial T$ vs $T$ for $B={B[j]}T$ and $ x = {[float(x) for x in list_dopings]}$", ylim= 0, xlim = 0)
pt.make_legend(ax[0], fig[0])
pt.make_dotlines(ax[0],sim,HS=False)
pt.make_legend(ax[1], fig[1])
pt.make_legend(axss[-1], fg)

names[0] = "R_H vs T"
names[1] = "rho0 vs T"
names[2] = "drho_dT vs T"

n = "MR_vs_T_diff_axes"

os.makedirs("./SCTIF/results propres/doping studies/png/", exist_ok=True)
os.makedirs("./SCTIF/results propres/doping studies/svg/", exist_ok=True)

fg.savefig("./SCTIF/results propres/doping studies/png/" + n + ".png", transparent=False)
fg.savefig("./SCTIF/results propres/doping studies/" + n + ".pdf", transparent=False)
fg.savefig("./SCTIF/results propres/doping studies/svg/" + n + ".svg", transparent=False)

i=0
# Save Figures


for n in names:
    if n!="":
        n = f"LCCO_FLEX={sims[0].FLEX}_SCTIF_" + n
        fig[i].savefig("./SCTIF/results propres/doping studies/" + "doping_" + n + ".pdf")
        fig[i].savefig("./SCTIF/results propres/doping studies/png/"+ "doping_"+ n + ".png", transparent=False)
        fig[i].savefig("./SCTIF/results propres/doping studies/svg/"+ "doping_"+ n + ".svg", transparent=False)
    i+=1

plt.show()
