from bayes_opt import BayesianOptimization
from Code_SCTIF.BTE_routines_SCTIF import *
import Code_SCTIF.metals_SCTIF as met
from Code_SCTIF.plot_results_SCTIF import *
from simul_SCTIF import Simulation
import matplotlib.pyplot  as plt
import Code_SCTIF.plot_tau_phi_SCTIF as pt
import pandas as pd
import os

plt.style.use("pacoty.mplstyle")
list_dopings = ["0.14","0.15","0.16","0.18"]
cuprate = "LCCO"
list_met = []
sims = []

og_data_rho0 = pd.read_csv('rho_0_vs_x_40K_exp.csv')
dopings_rho0_exp = og_data_rho0["x"]
rho0_exp = og_data_rho0["rho0"]

og_data_MR = pd.read_csv('MR_20T_40K_vs_x_exp.csv')
dopings_exp = og_data_MR["x"]
MR_exp = og_data_MR["MR"]

for x in list_dopings:
    metal = met.Metal(cuprate, p=x)
    list_met.append(metal)
    sims.append(Simulation(metal = metal))

def get_resistivity(sim, B=20, T = 40) :
    s_coeff = sim.metal.s_coeff
    tau_coeff = sim.tau_coeff
    tau_choice = sim.tau_choice
    N_k = 200
    dt = 1e-15
    N_tau = 70
    tau = get_tau_iso(T, s_coeff, tau_choice, tau_coeff, False)
    t_lim = tau*N_tau

    rho_xx, rho_xy, rho_yx, rho_yy = compute_resistivity(B, T, N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, False, sim.FLEX)

    return rho_xx, rho_xy, rho_yx, rho_yy

def compute_MR(B: float,T: float, sim: Simulation) :
    """compute a single MR value @ B and T for sim (= metal + doping + scattering rate)"""
    rho_xx, rho_xy, rho_yx, rho_yy= get_resistivity(sim,B, T)
    rho_xx_0, rho_xy, rho_yx, rho_yy = get_resistivity(sim,0, T)
    return (rho_xx - rho_xx_0)/rho_xx_0

def plot_MR_vs_x(ax, dopings, B=20, T=40):
    sims = []
    for x in dopings:
        metal = met.Metal(cuprate, p=x)
        sims.append(Simulation(metal = metal))
    MR_x = []
    for i in range(len(sims)):
        MR_x.append(compute_MR(B,T, sims[i]))
        dopings[i] = float(dopings[i])
    ax.scatter(dopings, MR_x)
    ax.set(xlabel = "Doping, $x$", ylabel = r'$\Delta\rho/\rho_0$',
        title = rf"$\Delta\rho/\rho_0$(20T, 40K) $-$ LCCO {sims[i].metal.p} \texttt{{FLEX={sims[i].FLEX}}}")


def plot_rho_0_vs_x(ax, dopings, T=40):
    sims = []
    for x in dopings:
        metal = met.Metal(cuprate, p=x)
        sims.append(Simulation(metal = metal))
    rho_0 = []
    for i in range(len(sims)):
        rho_0.append(get_resistivity(sims[i], 0, T)[0])
        dopings[i] = float(dopings[i])
    ax.scatter(dopings, rho_0)
    ax.set(xlabel = "Doping, $x$", ylabel = r'$\rho_{xx}\ (\mu\Omega cm)$',
        title = rf"$\rho_{{xx}}$(0T, 40K) $-$ LCCO {sims[i].metal.p} \texttt{{FLEX={sims[i].FLEX}}}")

opt = False

if not opt:
    N=2
    figs = []
    names = []
    for i in range(N) :
        figs.append(plt.subplots())
        names.append("")
    ax = [f[1] for f in figs]
    fig = [f[0] for f in figs]

    plot_MR_vs_x(ax[0], list_dopings)
    ax[0].scatter(dopings_exp, MR_exp, label = "exp results")
    pt.make_legend(ax[0],fig[0])
    names[0] = "MR_vs_x"

    plot_rho_0_vs_x(ax[1], list_dopings)
    ax[1].scatter(dopings_rho0_exp, rho0_exp, label = "exp results")
    pt.make_legend(ax[1],fig[1])
    names[1] = "rho0_vs_x"
    plt.show()
    i=0

    os.makedirs("./SCTIF/results propres/doping studies/png/", exist_ok=True)
    os.makedirs("./SCTIF/results propres/doping studies/svg/", exist_ok=True)

    for n in names:
        if n!="":
            n = f"LCCO_FLEX={sims[0].FLEX}_SCTIF_" + n
            fig[i].savefig("./SCTIF/results propres/doping studies/" + "doping_" + n + ".pdf")
            fig[i].savefig("./SCTIF/results propres/doping studies/png/"+ "doping_"+ n + ".png", transparent=False)
            fig[i].savefig("./SCTIF/results propres/doping studies/svg/"+ "doping_"+ n + ".svg", transparent=False)
        i+=1
else :
    ## BAYESIAN OPT OF DOPING PARAMETERS
    def black_box(g_x, imp_x):
        B = 20
        T = 40
        metal = met.Metal(cuprate, p="0.14")
        sim= Simulation(metal = metal, tau_choice = "doping opt", tau_coeff_sans_HS = np.array([27.e12, g_x, imp_x, 0., 0.75e12, 0.0112e12, 12, 0, 0, 1] ))
        MR = compute_MR(B,T, sim)
        rho = get_resistivity(sim, B, T)[0]
        MR_exp = 0.0315
        rho_exp = 28.42
        return -np.abs(MR - MR_exp) - np.abs(rho - rho_exp)
    min_max = [(0.05,2), (0.7,1)]
    # Bounded region of parameter space
    name_args = ["g_x","imp_x"]
    pbounds = {name_args[i]:min_max[i] for i in range (len(name_args))}

    optimizer = BayesianOptimization(
        f=black_box,
        pbounds=pbounds,
        verbose=2, # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state=1,
    )

    optimizer.maximize(
        init_points=6,
        n_iter=15,
    )

    print(optimizer.max)

plt.show()
