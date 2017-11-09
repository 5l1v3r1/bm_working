# front for box-model
import numpy as np
import matplotlib.pyplot as plt

from main_time_loop_swe import main_time_loop as mtl
from gas_phase_reac_ode_make_swe import gas_phase_reac_ode_make as gprom

# get ode information for gas-phase
[stor, reac_rate] = gprom()
# times to output results at (s)
res_tim = np.linspace(0, 3.6e3, 101)
# call on main time loop for chamber processes
[C_g_all, C_p_all] = mtl(res_tim, stor, reac_rate)

# plot change in gas-phase concentration (ug/m3) 
# of components over time
p1, = plt.plot(res_tim, C_g_all[:,0], '--g')
p2, = plt.plot(res_tim, C_g_all[:,1], '--b')
p3, = plt.plot(res_tim, np.sum(C_p_all,1), '--r')


plt.legend([p1, p2, p3],["VOC","O","PM"],loc='upper left')
plt.ylim(ymin = -0.1)
plt.xlabel(r'$t$ (s)')
plt.ylabel(r'$C\, (\mu g\, m^{-3})$')
plt.show()