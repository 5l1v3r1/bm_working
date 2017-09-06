# function to track particle number size distribution using Full-Stationary 
# size structure (p. 462 of Jacobson 2005).

import numpy as np

def full_stat_main(n0, s0, m0, rho, Vref):

	# input:
	# n0 - initial particle number concentration per size bin 
	# (# particle m^{-3} (air))
	# s0 - initial volume bounds per size bin (m^3) (1st dim.)
	# m0 - initial particle phase mass per size bin (1st dim.) 
	# (g m^{-3} (air))
	# rho - particle phase component (2nd dim.) densities (g m^{-3}), 
	# repeated across size bins (1st dim.)
	# Vref - reference volume of single particles in each fixed 
	# size bin (m^3)
	# output:
	# n1 - end of time step particle number concentration per size bin
	# (# particle m^{-3} (air))
	# m1 - end of time step mass per size bin (g m^{-3} (air)) 

	# new volume of total particles per original size bin 
	# (m^3 (particle) m^{-3} (air)) (row array)
	VT = m0*(1.0/rho)
		
	# break and display error if volume array and number concentration
	# of particles array not in agreement
	if np.sum(VT>0)!=np.sum(n0[n0[:, 0]>0, 0]):
		print 'disagreement between particle volumes and particle 				number concentration, line 29 of full_stat_main2.py'
		return
	# new volume of single particles per original size bin 
	# (m^3 (particle) m^{-3} (air)) (row array)
	V = np.zeros((VT.shape[0]))
	V[VT>0] = VT[VT>0]/n0[n0[:, 0]>0, 0]
	 	
	# repeat fixed bin volume bounds array over variable bins 
	# (fixed bin bounds in 1st dim. and variable bins in 2nd dim.)
	smat0 = np.transpose(np.tile(s0, (V.shape[0], 1)))
	
	# indices of bin bounds less than or equal to current volumes
	ib = smat0<V	
	# index of fixed size bin particle volumes fall inside
	ib = np.sum(ib, axis=0)-1
	# repeat bin indices over fixed size bins (size bins in 1st dim. and
	# bin indices in 2nd dim.)
	ib = np.tile(ib, (Vref.shape[0], 1))
	
	# reference size bin index (1st dim.), repeated over variable volume 
	# array (2nd dim.)
	refi = np.arange(0, Vref.shape[0], 1)
	refi = np.transpose(np.tile(refi, (V.shape[0], 1)))
	
	# find where bin indices for volumes and reference indices are in 
	# agreement
	agree = (refi == ib)
	del refi, ib, smat0
	# repeat particle number (# m^{-3} (air)) (2nd dim.) over fixed 
	# size bins (1st dim.)
	nrep = np.transpose(np.tile(n0, (1, Vref.shape[0])))
	# multiply number concentration by truth matrix and sum numbers per 
	# size bin (# m^{-3})
	n1 = np.sum(agree*nrep, axis=1)	
	del agree, nrep
	print n1	
	# particle volume per size bin (m^3 (particle) m^{-3} (air))
	V1 = n1*Vref
	# particle mass per size bin (g m^{-3} (air))
	m1 = V1*rho 
		
	return n1, m1
