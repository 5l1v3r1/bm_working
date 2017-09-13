# Copyright Notice: This code is in Copyright.  Any use leading to publication or 
# financial gain is prohibited without the permission of the authors Simon O'Meara and # David Topping: simon.omeara@manchester.ac.uk.  First published 2017.

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
# along with box_model (see COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

# -----------------------------------------------------------------------------------
# function to treat a dynamic number size distribution by the Full-Moving 
# Size Structure method (p. 465 Jacobson 2005).

def full_move_main(n0, m0, rho):
	
	# -------------------------------------------------------------------	
	# input:
	# n0 - number of particles per size bin (#/m3(air))
	# m0 - mass of each component (2nd dim.) per size bin (1st dim.) 
	# (g/m3(air)) 
	# rho - component (2nd dim.) density (g/m3 (condensed phase)) 
	# repeated across size bins (1st dim.)

	# output:
	# n1 - end of time step particle number concentration per size bin 
	# (# particle/m3(air))
	# m1 - end of time step mass per size bin per particle (g/m3 (air))	

	# notes:
	# Full-moving is when particles do not transfer between size bins.
	# Instead they stay in one size bin and the volume of particles
	# in the size bin changes.  Although mass of initial particles is
	# conserved in this way, it can mean there is no size bin available
	# when new particles enter the system or existing particles coagulate.
	# Since no rebinning is required, this code could be useless, it is
	# only provided in case the volume (total or
	# individual) of particles per size bin  is required 
	
	# ------------------------------------------------------------------
	# get new volume of single particles in size bins:

	# new total volume of particles per size bin (m3(particle)/m3(air))
	# (sum volume of individual components)
	VT = np.sum(m0*(1.0/rho),0)
	# volume of indivual particles per size bin 
	# (m3(particle)/m3(air))
	V1 = VT/n1
	

	return m1, n1
