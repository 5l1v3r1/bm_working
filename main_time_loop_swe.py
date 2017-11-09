import numpy as np
from scipy.integrate import odeint as int

# function for looping through time steps to 
# solve chamber processes

def main_time_loop(res_tim, stor, reac_rate):

	# -----------------------------------------------------
	# inputs:
	# res_tim - times to output results at (s)
	# stor - stoichiometry of gas-phase reactions
	# reac_rate - reaction rate of gas-phase reactions
	# -----------------------------------------------------
	
	# initial gas-phase concentrations (ug/m3 (air))
	C0 = np.array((3.0e0, 1.0e-4, 0.0, 0.0))
	# initial particle-phase concentrations (ug/m3 (air))
	C_p0 =  np.array((0.0, 0.0, 0.0, 0.0))
	# record results
	C_g_all = np.zeros((np.shape(res_tim)[0], (np.shape(stor))[0]))
	C_p_all = np.zeros((np.shape(res_tim)[0], (np.shape(stor))[0]))
	C_g_all[0, :] = C0
	C_p_all[0, :] = C_p0
	
	# ode function for gas-phase reaction
	def gas_reac(C0, tspan, reac_rate, stor):
	
		dcdt = np.zeros((4))
		# solve new gas-phase concentrations (ug/m3)
		dcdt[0] = np.sum(reac_rate[0,:]*(C0[:]**stor[0,:]))
		dcdt[1] = 0.0 # keep oxidant concentration constant
		dcdt[2] = np.sum(reac_rate[2,:]*(C0[:]**stor[2,:]))
		dcdt[3] = np.sum(reac_rate[3,:]*(C0[:]**stor[3,:]))
		return dcdt
	
	# ode function for gas-particle partitioning
	def part_ode(C0, tspan, k_qmt, Csta_imt):
		
		dcdt = -k_qmt*(C0-Csta_imt)
		return dcdt
	
	# loop through output result times
	for itn in range(1, np.shape(res_tim)[0]):

		# time span to integrate over
		tspan = np.array((0, (res_tim[itn]-res_tim[itn-1])))
		
		# use integration to estimate new gas-phase 
		# concentrations (ug/m3)
		sol = int(gas_reac, C0, tspan, args=(reac_rate, stor))
		Cn = sol[-1,:]
		
		# mass transfer coefficient for partitioning (m3/s)
		k_qmt = 1.0
		# effective saturation concentrations (ug/m3)
		Csta_imt = np.array((1.0e6, 1.0e6, 1.0e-6, 1.0e-4))
		# solve partitioning equation
		sol = int(part_ode, Cn, tspan, args=(k_qmt, Csta_imt))
		Cnn = sol[-1,:]
		
		# prepare particle-phase concentration for change
		C_pn = C_p0
		# ensure particle-phase concentration not gone negative
		neg_ind = (C_p0-(Cnn-Cn)<0.0)
		
		C_pn[neg_ind] = 0.0
		Cnn[neg_ind] = Cn+C_p0
		
		# index of particle-phase components not exhausted
		pos_ind = (C_p0-(Cnn-Cn)>=0.0) 
		# change particle-phase concentration (ug/m3 (air))
		C_pn[pos_ind] = C_pn[pos_ind]-(Cnn[pos_ind]-Cn[pos_ind])
		
		# reset initial variables
		C0 = Cnn
		C_p0 = C_pn
		C_g_all[itn, :] = Cnn
		C_p_all[itn, :] = C_pn
	return C_g_all, C_p_all	
	