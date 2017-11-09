# function to create the matrices required for gas-phase chemistry ode solution
import numpy as np
from Parse_eqn_file_SPO_test_swe import Parsing as Pars

def gas_phase_reac_ode_make():

	# call parsing script
	filename = 'test_eqn_swe.eqn'
	[stor_mat, reac_mat, stop_mat, prod_mat, kine_mat] = Pars(filename)
	
	# unique components
	unicr = np.unique(reac_mat)
	unicp = np.unique(prod_mat)
	unic = np.append(unicr,unicp)
	
	nuni = len(unic)
	
	# matrix for storing reactant stoichiometries
	stor = np.zeros((nuni,nuni))
	
	# empty results matrix for total reaction rate coefficient
	reac_rate = np.zeros((nuni,nuni))
	
	# number of equations
	eqnn = (np.shape(reac_mat))[0]
	
	for icom1 in range(0, nuni): # first component loop
	
		for icom2 in range(0, nuni): # second component loop
				
			# equation where first component reacts
			indcom1 = np.sum((reac_mat==str(unic[icom1])),1)
			# equation where second component reacts
			indcom2 = np.sum((reac_mat==str(unic[icom2])),1)
			
			indcom1[indcom1>0] = 1 
			indcom2[indcom2>0] = 1 
			match_eqs = ((indcom1+indcom2)==2)
			
			# equation where first component produced
			indcom3 = np.sum((prod_mat==str(unic[icom1])),1)
			indcom3[indcom3>0] = 1 
			match_eqs2 = ((indcom3+indcom2)==2)
			
			for eqli in range(0, eqnn): # equation loop
				
				if match_eqs[eqli]==1: # relevant equation
					# reaction rate coefficient
					reac_rate[icom1,icom2] = reac_rate[icom1,icom2]+-kine_mat[eqli,0]
					# stoichiometry of products in this equation
					stor[icom1, icom2] = stor[icom1, icom2]+stor_mat[eqli,icom2]
					
				if match_eqs2[eqli]==1: # relevant equation
					# reaction rate coefficient
					reac_rate[icom1,icom2] = reac_rate[icom1,icom2]+kine_mat[eqli,0]
					# stoichiometry of products in this equation
					stor[icom1, icom2] = stor[icom1, icom2]+stor_mat[eqli,icom2]

	return stor, reac_rate	
