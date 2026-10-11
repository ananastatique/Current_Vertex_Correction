from simul import Simulation
import Code_SCTIF.metals as met
import numpy as np


sim_spin = Simulation()
sim_spin.sweep_B()
sim_spin.sweep_T(BT = np.array([1, 2, 5, 10]))
