# Copyright Notice: This code is in Copyright.  Any use leading to publication 
# or financial gain is prohibited without the permission of the authors 
# Simon O'Meara and David Topping: simon.omeara@manchester.ac.uk.  First 
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

# ------------------------------------------------------------------------------
# function to track particle number size distribution using Moving-centre 
# size structure (p. 466 Jacobson 2005).

import numpy as np

def move_cent(n0, s0, m0, rho, Vref):

	# ----------------------------------------------------------------
	# input:
	# n0 - initial particle number concentration per size bin 
	# (# particle/m3 (air))
	# s0 - initial volume bounds per size bin (m3) (1st dim.)
	# m0 - initial particle phase mass per size bin (1st dim.) per
	# component (2nd dim.) (g/m3 (air))
	# rho - particle phase component (2nd dim.) densities (g/m3), 
	# repeated across size bins (1st dim.)
	# Vref - reference volume of single particles in each fixed 
	# size bin (m^3)

	# output:
	# n1 - end of time step particle number concentration per size bin
	# (# particle/m3 (air))
	# m1 - end of time step mass per component (1st dim.) per size bin 
	# (2nd dim.) (g/m3 (air)) 
	
	# notes:
	# in the moving-centre size structure, particles are allocated to the
	# fixed size bin who's bounds they fall within.  The average volume of
	# single particles in a size bin is allowed to vary, which is the
	# difference to the full-stationary size structure 

	# ----------------------------------------------------------------
	# get new volume of single particles per size bin:

	# new volume of total particles per original size bin 
	# (m3(particle)/m3 (air)) (row array) (sum over components 
	# in 1st dim.)
	VT = np.sum(m0*(1.0/rho),1)
	VT.shape = (VT.shape[0],1) # ensure 2D not 1D 
	
	# new volume of single particles per original size bin 
	# (m^3 (particle) m^{-3} (air)) (row array)
	V = np.zeros((VT.shape[0], 1))
	V[VT>0] = VT[VT>0]/n0[n0[:, 0]>0, 0]
	
	# ------------------------------------------------------------------
	# find size bins where particles belong:
 	
	# repeat fixed bin volume bounds array over variable bins 
	# (fixed bin bounds in 1st dim. and variable bins in 2nd dim.)
	smat0 = np.transpose(np.tile(s0, (V.shape[0], 1)))
	
	# indices of bin bounds less than or equal to current volumes
	ib = smat0<np.transpose(V)	
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
	
	# -----------------------------------------------------------------
	# allocate numbers of particles to their appropriate bin and find
	# new mass concentration (g/m3 (air)):
	
	# repeat particle number (#/m3(air)) (2nd dim.) over fixed 
	# size bins (1st dim.)
	nrep = np.transpose(np.tile(n0, (1, Vref.shape[0])))	
	# multiply number concentration by truth matrix and sum numbers per 
	# size bin (#/m3(air))
	n1 = np.sum(agree*nrep, axis=1)	
	
	# repeat truth matrix over components (1st dim), and have variable
	# size bins in 2nd dim. and fixed size bins in 3rd dim.	
	agree = np.repeat(np.transpose(agree)[np.newaxis,:,:],rho.shape[0],0)
	m0 = np.repeat(np.transpose(m0)[:,:,np.newaxis],Vref.shape[0],2)
	# get new component masses by summing the truth and component mass 
	# array and sum over variable size bins to leave a matrix with
	# components in 1st dim. and fixed size bins in 2nd dim.
	m1 = np.sum(agree*m0,1)
	
			
	return n1, m1
