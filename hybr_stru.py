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

# ------------------------------------------------------------------------------# function to track particle number size distribution using hybrid 
# size structure (p. 420 Jacobson 2000).

import numpy as np
from full_stat import full_stat
from full_move import full_move

def hybrid_main(n0, s0, mnv0, msv0, rhonv, rhosv, sbn, Vref, Vsv0):

	# ------------------------------------------------------------------
	# input:
	# n0 - initial particle number concentration per size bin 
	# (# particle m/3 (air))
	# s0 - initial volume bounds per size bin (m^3) (1st dim.)
	# mnv0 - initial nonvolatile particle phase mass per size bin 
	# (1st dim.) (g/m3 (air))
	# msv0 - initial semi-volatile particle phase mass per size bin
	# (1st dim.) (g/m3 (air))
	# rho - particle phase component (2nd dim.) densities (g/m3), 
	# repeated across size bins (1st dim.)
	# sbn - number of size bins
	# Vref - reference volume of single particles in each fixed 
	# size bin (m3)
	# Vsv0 - initial volumes of semi-volatile material per fully moving
	# size bin (m3) (summed across all particles in that bin)

	# output:
	# n1 - end of time step particle number concentration per size bin
	# (# particle/m3 (air))
	# mnv1 - end of time step mass of involatile components per size bin 
	# (g/m3 (air)) 
	# msv1 - end of time step mass of semi-volatile components per size bin 
	# (g/m3 (air)) 

	# notes:
	# core (involatile part) of particles treated with full-stationary 
	# whilst shell (volatile part) treated with full-moving.  Numerical
	# diffusion a disadvantage.  The size bin order of sv corresponds to 
	# that of the nv, so that mnv1 in bin #n has the msv1 in bin #n 
	# condensed onto it
	
	# ------------------------------------------------------------------
	# call on full-stationary method to find the number concentration of
	# particles per size bin (/m3(air)) and mass concentration 
	# (g/m3 (air)) of involatile components per size bin 
	[nnv1, mnv1] = full_stat(n0, s0, mnv0, rhonv, Vref)	

	# call on full-moving method to calculate average semi-volatile volume 
	# per particle per size bin (m3(particle)/m3(air)) and average 
	# mass of semi-volatile per particle per size bin 
	# (g(condensed phase)/m3(air)).  The size bin order of sv is aligned
	# with that of the nv, so that mnv1 in bin #n has the msv1 in bin #n 
	# condensed onto it
	[msv1, nsv1] = full_move(n1, msv0, rhosv)
	
	return n1, mnv1, msv1
