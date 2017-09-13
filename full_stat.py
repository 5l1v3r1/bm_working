# Copyright Notice: This code is in Copyright.  Any use leading to publication 
# or financial gain is prohibited without the permission of the authors 
# Simon O'Meara and # David Topping: simon.omeara@manchester.ac.uk.  First 
# published 2017.

# This file is part of box_model

# box_model is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# box_model is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with box_model (see COPYING.txt).  If not, see 
# <http://www.gnu.org/licenses/>.

# ------------------------------------------------------------------------------# function to track particle number size distribution using Full-Stationary 
# size structure (p. 462 Jacobson 2005).

import numpy as np

def full_stat(n0, s0, m0, rho, Vref):

	# -----------------------------------------------------------------
	# input:
	# n0 - initial particle number concentration per size bin 
	# (# particle /m3 (air))
	# s0 - initial volume bounds per size bin (m3) (1st dim.)
	# m0 - initial particle phase mass per size bin (1st dim.) per
	# component (2nd dim) (g/m3 (air))
	# rho - particle phase component (2nd dim.) densities (g/m3), 
	# repeated across size bins (1st dim.)
	# Vref - reference volume of single particles in each fixed 
	# size bin (volume at centre of bin) (m3)
	
	# output:
	# n1 - end of time step particle number concentration per size bin
	# (# particle/m3 (air))
	# m1 - end of time step mass per component (1st dim.) per size bin 
	# (2nd dim.) (g m/3 (air)) 
	
	# notes:
	# The full-stationary size structure keeps the volume of single 
	# particles in each size bin constant. It works by finding the new 
	# volume of single particles based on 
	# their mass and density, then allocating them to the size bins who's
	# volume bounds they fall within.  The allocated particles are then 
	# assumed to have the fixed volume set for that size bin.  Therefore, 
	# number is conserved but not mass.     
	
	# -----------------------------------------------------------------
	# mass fractions per component per variable size bin
	
	# tile the total mass of each size bin (1st dim.) over components
	# (2nd dim.)
	mTmat = np.tile(np.sum(m0,1), (m0.shape[1],1))
	# mass fraction per component (1st dim.) per variable size bin
	# (2nd dim.)
	mTmat[mTmat==0.0] = 1.0e6 # prevent error in division by setting high
	mcrat = np.transpose(m0)/mTmat
	# tile mcrate over fixed size bins in 3rd dim.
	mcrat = np.repeat(mcrat[:,:,np.newaxis],Vref.shape[0],2)	

	# -----------------------------------------------------------------
	# get new volume of single particles in size bins:

	# new total volume of particles per original size bin 
	# (m^3(particle)/m3(air)) (row array) (sum across volume of individual
	# components)
	VT = np.sum(m0*(1.0/rho),1)	
	
	# new volume of single particles per original size bin 
	# (m^3(particle)/m3(air)) (row array)
	V = np.zeros((VT.shape[0]))
	V[VT>0] = VT[VT>0]/n0[n0[:, 0]>0, 0] # (only where n>0 to avoid error)
	
	# ----------------------------------------------------------------	 	# find size bins where particles belong:

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
	
	# reference size bin index (1st dim.), repeated over variable size bins
	#  (2nd dim.)
	refi = np.arange(0, Vref.shape[0], 1)
	refi = np.transpose(np.tile(refi, (V.shape[0], 1)))
	
	# find where bin indices for volumes and reference indices are in 
	# agreement
	agree = (refi == ib)
	del refi, ib, smat0
	
	# -----------------------------------------------------------------
	# allocate numbers of particles to their appropriate bin and find
	# new mass concentration (g/m3 (air)) 

	# repeat particle number (# /m3 (air)) (2nd dim.) over fixed 
	# size bins (1st dim.)
	nrep = np.transpose(np.tile(n0, (1, Vref.shape[0])))
	# multiply number concentration by truth matrix and sum numbers per 
	# size bin (# /m3)
	n1 = np.sum(agree*nrep, axis=1)	
	
	# -----------------------------------------------------------------
	# number fraction from each variable bin (2nd dim.) moving into each
	# fixed bin (1st dim.)
	n2 = n1.reshape(n1.shape[0],1) # 2 dim. version rather than 1 dim.	
	n2[n2==0] = 1.0e6 # set very high to prevent error when dividing
	nfpb = np.transpose(agree*nrep/(n2))
	# repeat over components in 1st dim. (variable size bins in 2nd dim.
	# and fixed size bins in 3rd dim.)
	nfpb = (np.repeat(nfpb[np.newaxis,:,:],m0.shape[1],0))
	
	# multiply number fraction by component mass fraction and sum across
	# variable size bins to get the new component mass fractions in fixed
	# size bins		
	mf1= np.sum(nfpb*mcrat,1)		

	# ----------------------------------------------------------------	
	
	# new total particle volume per size bin (m3(particle)/m3(air))
	VT = n1*Vref
	# average density (g/m3(particle)) per size bin
	rho1 = rho.reshape(rho.shape[0],1) # change from 1 dim. to 2 dim.
	av_rho = np.sum(mf1*rho1,0) # size bins in 2nd dim.
	# new total mass per fixed size bin (g(particle)/m3(air))
	mT = VT*av_rho
	# new mass per components (1st dim.) per size bin (2nd dim.)
	m1 = mT*mf1	

	print m1
	return 
		
	return n1, m1
