# function to parse .eqn files output from kpp
import numpy as np

# define function
def Parsing(filename):

	# inputs
	# filename - string containing name of file to parse

	# open file in read mode
	with open(filename,'r') as file:
			
		# emptry matrices
		stor_mat = np.ones((1,1))	
		stop_mat = np.ones((1,1))
		reac_mat = np.zeros((1,1), dtype='U25')	
		prod_mat = np.zeros((1,1), dtype='U25')	
		kine_mat = np.zeros((1,1))	
		# equation line count
		eqlc = -1
		# count on LHS important equation elements
		elcl = -1
		# count on RHS important equation elements
		elcr = -1
		# count on kinetics important equation elements
		elck = -1
		
		# loop through lines
		for line in file:
			si = line.find('>')
			if si>-1:
			
				
				# count on whether RHS reached
				RHSc = 0
					
				# split string into components 
				# separated by white space
				line1 = str.split(line)
				
				# equation part loop
				for stri in range(1, len(line1)):
					cstr = line1[stri]
					
					
					
					# mark for RHS of equation
					if cstr=='=':
						RHSc = 1
						continue
					# mark for kinetic information
					if cstr==':':
						RHSc = 2
						continue
					
					if str.isdigit(cstr[0]) or str.isalpha(cstr[0]):
						if RHSc==0:
							elcl = elcl+1
					if str.isdigit(cstr[0]) or str.isalpha(cstr[0]):
						if RHSc==1:
							elcr = elcr+1
					if RHSc==2:
						elck = elck+1
					
					dig_sea=0 # count on characters in strings
					
					# character loop
					for cind in range(0, len(cstr)):
						
						if cind+dig_sea>=len(cstr):
							break
						
						
						
						# reactants and their stoichometry
						if RHSc==0:
							
							if eqlc>=stor_mat.shape[0]:
								stor_mat = np.append(stor_mat,np.ones((1,stor_mat.shape[1])),0)
								reac_mat = np.append(reac_mat,np.zeros((1,reac_mat.shape[1]),dtype='U25'),0)
							if elcl>=stor_mat.shape[1]:
								stor_mat = np.append(stor_mat,np.ones((stor_mat.shape[0],(elcl-stor_mat.shape[1]+1))),1)				
								reac_mat = np.append(reac_mat,np.zeros((reac_mat.shape[0],(elcl-reac_mat.shape[1]+1)),dtype='U25'),1)				
							
							if str.isdigit(cstr[cind]):
								
								dig_sea=0
								# look forwards to see how many digits are in this number
								if (cind+dig_sea+1)<(len(cstr)):
									while str.isdigit(cstr[cind+dig_sea+1]) or cstr[cind+dig_sea+1]=='.' or cstr[cind+dig_sea+1]=='E' or cstr[cind+dig_sea+1]=='-':
										dig_sea = dig_sea+1
										break
								
								stor_mat[eqlc,elcl] = cstr[cind:cind+dig_sea+1]
							
							if str.isalpha(cstr[cind]) or cstr[cind]=='(':
								dig_sea=0
								# look forwards to see how many letters are in this reactant
								if (cind+dig_sea)<(len(cstr)):
									
									while (cind+dig_sea+1) < (len(cstr)):
										dig_sea = dig_sea+1
										
							
								reac_mat[eqlc,elcl] = str(cstr[cind:cind+dig_sea+1])
								
						# products and their stoichometry
						if RHSc==1 and cstr!='=':		
							
							if eqlc>=stop_mat.shape[0]:
								stop_mat = np.append(stop_mat,np.ones((1,stop_mat.shape[1])),0)
								prod_mat = np.append(prod_mat,np.zeros((1,prod_mat.shape[1]),dtype='U25'),0)
							if elcr>=stop_mat.shape[1]:
								stop_mat = np.append(stop_mat,np.ones((stop_mat.shape[0],(elcr-stop_mat.shape[1]+1))),1)				
								prod_mat = np.append(prod_mat,np.zeros((prod_mat.shape[0],(elcr-prod_mat.shape[1]+1)),dtype='U25'),1)				
							
							if str.isdigit(cstr[cind]):
								
								dig_sea=0
								# look forwards to see how many digits are in this number
								if (cind+dig_sea+1)<(len(cstr)):
									while str.isdigit(cstr[cind+dig_sea+1]) or cstr[cind+dig_sea+1]=='.' or cstr[cind+dig_sea+1]=='E' or cstr[cind+dig_sea+1]=='-':
										dig_sea = dig_sea+1
										break
								
								stop_mat[eqlc,elcr] = cstr[cind:cind+dig_sea+1]
							
							if str.isalpha(cstr[cind]) or cstr[cind]=='(':
								dig_sea=0
								# look forwards to see how many letters are in this reactant
								if (cind+dig_sea)<(len(cstr)):
									
									while (cind+dig_sea+1) < (len(cstr)):
										dig_sea = dig_sea+1
										
								
								prod_mat[eqlc,elcr] = str(cstr[cind:cind+dig_sea+1])
								break # from character loop
						# kinetics
						if RHSc==2:
							
							if eqlc>=kine_mat.shape[0]:
								kine_mat = np.append(kine_mat,np.ones((1,kine_mat.shape[1])),0)
							if elck>=kine_mat.shape[1]:
								kine_mat = np.append(kine_mat,np.ones((kine_mat.shape[0],(elck-kine_mat.shape[1]+1))),1)				
							
							dig_sea=0
							
							while cstr[cind+dig_sea]!=')':
								dig_sea = dig_sea+1
							
							kine_mat[eqlc,elck] =  cstr[cind+1:cind+dig_sea]
							break # onto next equation component
			
	return stor_mat, reac_mat, stop_mat, prod_mat, kine_mat



		
