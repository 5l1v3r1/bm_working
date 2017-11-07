import numpy as np

# open file in read mode
with open('small_strato2.eqn','r') as file:
		
	# emptry matrices
	stoi_mat = np.ones((1,1))		

	# equation line count
	eqlc = -1
	# count on important equation elements
	elc = -1
	# loop through lines
	for line in file:
		si = line.find('>')
		if si>-1:
			eqlc = eqlc+1
			# count on whether RHS reached
			RHSc = 0
			# count on whether stoichometry or formula seen
			stoi = 1
			# split string into components 
			# separated by white space
			line1 = str.split(line)
			# equation part loop
			for stri in range(1, len(line1)):
				cstr = line1[stri]
				if str.isdigit(cstr[0]) or str.isalpha(cstr[0]):
					elc=elc+1
				
				# character loop
				for cind in range(0, len(cstr)):
					
					if cstr=='=':
						RHSc=1
					if str.isalpha(cstr[cind]):
						stoi=1
							 
					# stoichometry of reactants
					if RHSc==0 and str.isdigit(cstr[cind]) and stoi==1:
						dig_sea=0
						# look forwards to see how many digits are in this number
						if (cind+dig_sea)<(len(cstr)-1):
							while str.isdigit(cstr[cind+dig_sea+1]) or cstr[cind+dig_sea+1]=='.' or cstr[cind+dig_sea+1]=='E' or cstr[cind+dig_sea+1]=='-':
								dig_sea = dig_sea+1
						if eqlc>=stoi_mat.shape[0]:
							stoi_mat = np.append(stoi_mat,np.ones((1,stoi_mat.shape[1])),0)
						if elc>=stoi_mat.shape[1]:
							stoi_mat = np.append(stoi_mat,np.ones((stoi_mat.shape[0],(elc-stoi_mat.shape[1]+1))),1)				
						
						print cstr[cind+dig_sea]
						stoi_mat[eqlc,elc] = cstr[cind+dig_sea]
			
print stoi_mat


		
