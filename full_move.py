# function to treat a dynamic number size distribution by the Full-Moving Size Structure method (p. 465 Jacobson 2005).

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
