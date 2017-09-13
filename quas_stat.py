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

# ------------------------------------------------------------------------------# function to track particle number size distribution using Quasistationary 
# size structure (p. 465 Jacobson 2005).

import numpy as np

def quas_stat(n0, s0, m0, rho, Vref):

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
	# size bin (volume at centre of bin) (m3) (1st dim.)
	
	# output:
	# n1 - end of time step particle number concentration per size bin
	# (# particle/m3 (air))
	# m1 - end of time step mass per component (1st dim.) per size bin 
	# (2nd dim.) (g m/3 (air)) 
	
	# notes:
	# The quasistationary size structure keeps the volume of single 
	# particles in each size bin constant. It works by finding the new 
	# volume of single particles based on 
	# their mass and density.  The new volumes (vnew) are used to partition
 	# particles between adjoining size bins (those with fixed 
	# volumes below vnew (vj) and those above (vk)) when vj<=vnew<vk.
	# The allocated particles are then 
	# assumed to have the fixed volume set for that size bin (vj or vk).
	# Using eq. 13.34 of Jacobson 2005 ensures the partitioning conserves
	# number and volume (and therefore mass).  Like the full-stationary 
	# structure it is numerically diffusive.
	
	# -----------------------------------------------------------------
	# mass fractions per component per variable size bin	

	# tile the total mass of each size bin (1st dim.) over components 
	# (2nd dim.)
	mTmat = np.tile(np.sum(m0,1), (m0.shape[1], 1))	
	# mass fraction per variable size bin (1st dim.) per component 
	# (2nd dim.)
	mTmat[mTmat==0.0]=1.0e6 # prevent error by setting unrealistically high
	# components in 1st dim, variable size bin in 2nd dim.
	mcrat = np.transpose(m0)/mTmat 
	# tile mcrat over fixed size bins in 3rd dim.
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

	# repeat fixed bin centre volumes array over variable size bins 
	# (fixed bin bounds in 1st dim. and variable size bins in 2nd dim.)
	Vrefmat0 = np.transpose(np.tile(Vref, (V.shape[0], 1)))
		
	# indices of fixed bin centres less than or equal to current volumes
	ib = (Vrefmat0<=V)	
	# index of fixed size bin particle volumes are greater than or 
	# equal to (1D array)
	ib = np.sum(ib, axis=0)-1
	# ensure that the largest reference size bin is not being indexed (this
	# prevents error when 1 is added to indexes below to reference the next
	# size bin up.  The subtraction fo one from the index in this case, 
	# then leads to a zero in the corresponding delnj element, so no 
	# partitioning occurs to the lower size bin)
	ibi=(ib==Vref.shape[0]-1)
	ib[ibi]=ib[ibi]-1  
			
	# ----------------------------------------------------------------
	# get numbers to partition and do partition	
		
	# number to partition to lower size bin
	delnj=n0[:,0]*(Vref[ib+1]-V)/(Vref[ib+1]-Vref[ib])	
	# and to higher size bin
	delnk=n0[:,0]*(V-Vref[ib])/(Vref[ib+1]-Vref[ib])	
	# index set on a matrix, with set size bins in 1st dim. and variable
	# bins in 2nd dim.	
	ibmat = (np.tile(ib, (Vref.shape[0],1)))
	# index array for set size bins, spread over variable size bins
	ias = np.transpose(np.tile(np.arange(0,Vref.shape[0],1), 
			(V.shape[0],1)))
	# truth array for lower size bins
	tal = ibmat==ias
	# and for upper size bins
	tau = ibmat+1==ias
	
	# number to partition spread across fixed bins (1st dim.)
	delnj = np.tile(delnj, (Vref.shape[0],1))
	delnk = np.tile(delnk, (Vref.shape[0],1))

		
	# multiply by corresponding truth arrays and sum for each fixed size 
	# bin to get new numbers
	n1 = np.sum(tal*delnj,1)
	n1 = n1+np.sum(tau*delnk,1)

	# -------------------------------------------------------------------
	# get new mass fractions per component:

	# number fraction from each variable bin (2nd dim.) moving into each 
	# fixed bin (1st dim.)
	# total number being moved to lower bins
	tnpbl = np.transpose(np.tile((np.sum(tal*delnj,1)),(V.shape[0],1)))	
	# and to higher bins
	tnpbu = np.transpose(np.tile((np.sum(tau*delnk,1)),(V.shape[0],1)))
	
	# number fraction
	ig0 = tnpbl==0.0 # index where delnj equals zero	
	# where zero set very high to prevent error when dividing
	tnpbl[ig0]=1.0e6 	
	nfpbl = ((tal*delnj))/tnpbl # number fraction in lower bins 
	ig0 = tnpbu==0.0 # where equal zero
	tnpbu[ig0]=1.0e6 # same as above to prevent error
	nfpbu = ((tau*delnk))/tnpbu # number fraction in upper bins
	# empty array for overall number fraction
	nfpb = np.zeros((V.shape[0],Vref.shape[0]))
	# combine lower and upper bin results
	nfpb[:,:] = np.transpose(nfpbl+nfpbu)
	# repeat over components in 1st dim. (variable size bins in 2nd dim.
	# and fixed size bins in 3rd dim.)
	nfpb = np.repeat(nfpb[np.newaxis,:,:],m0.shape[1],0)

	# multiply number fraction by component mass fraction and sum across 
	# variable size bins to get the new component mass fractions in fixed 
	# size bins
	mf1 = np.sum(nfpb*mcrat,1)
	
	# ------------------------------------------------------------------
	# get new masses per component:

	# new total volume per size bin (m3(particle)/m3(air))
	VT = n1*Vref 
	# average density (g/m3(particle)) per size bin
	rho1 = rho.reshape(rho.shape[0],1) # change from 1 dim. to 2 dim.	
	av_rho = np.sum(mf1*rho1,0) # size bins spread across 2nd dim.
	# new total mass per size bin (g(particle)/m3(air))
	mT = VT*av_rho
	# new mass per components (1st dim.) per size bin (2nd dim.)
	m1 = mT*mf1
	
	return n1, m1
