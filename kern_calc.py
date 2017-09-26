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

# ------------------------------------------------------------------------------


# code to check coagulation kernel calculations

import numpy as np
import scipy.constants as si
from matplotlib import pyplot as plt
from Knud_calc import Knud_calc
from Reyn_numb import Reyn_numb as Reyn_numb
from vanW_fact import vanW_fact as vanW_fact

# ---------------------------------------------------------
# notes:
# aim is to reproduce Figs. 15.7 and 15.8 Jacobson 2005


# ---------------------------------------------------------
# constants:

# radii of particles (m)
radi = np.logspace(-8.0,-5.0)
# number of size bins
sbn = radi.shape[0]
radi.shape = (sbn,1) # 2D not 1D

# temperature (K)
T = 300.0
# Relative Humidity
RH = 0.01
# particle density (g/m3)
rho_p = 1.0e6 


# --------------------------------------------------------------------
# Brownian motion:

# call on function to determine the Knudsen no. and therefore flow 
# regime
[Kni, eta_a, rho_a, kin_visc] = Knud_calc(RH, T, radi)
# Cunningham slip-flow correction (15.30) with constant taken
# from text below 15.30 (dimensionless)
G = 1.0+Kni*(1.249+0.42*(np.exp(-0.87/Kni)))
# particle diffusion coefficient (15.29) (m2/s) (note the 
# Boltzmann constant has units (kg m2)/(s2 K), so *1e3 to
# convert to g from kg)
Dp = (((si.k*1.0e3)*T)/(6.0*np.pi*radi*eta_a))*G
# single particle mass (g) (volume and density product)
Mp = ((4.0/3.0)*np.pi*(radi**3.0))*rho_p
# thermal speed of particle (15.32) (m/s)
nu_p = ((8.0*(si.k*1.0e3)*T)/(np.pi*Mp))**0.5
# particle mean free path (15.34) (m)
lam_p = (8.0*Dp)/(np.pi*nu_p)

i = Kni<1.0 # continuum regime

# repeat over size bins
radim = np.repeat(radi,sbn,1)
Dpm = np.repeat(Dp,sbn,1)
im = np.repeat(i,sbn,1)
K_B = np.zeros((sbn,sbn))
# transposed version
radimt = np.transpose(radim)
Dpmt = np.transpose(Dpm)
# collision kernel (15.28) (m3/s)
K_B[im] = 4.0*np.pi*(radim[im]+radimt[im])*(Dpm[im]+Dpmt[im])


i = Kni<=10.0 
j = Kni>=1.0 # transition regime
i = i*j

# mean distance from centre of a sphere travelled by particles
# leaving sphere's surface and travelling lam_p (m) (15.34)
sig_p = ((2.0*radi+lam_p)**3.0-(4.0*radi**2.0+
	lam_p**2.0)**(1.5))/(6.0*radi*lam_p)-2.0*radi	
# repeat over size bins
im = np.repeat(i,sbn,1)
sig_pm = np.repeat(sig_p,sbn,1)
nu_pm = np.repeat(nu_p,sbn,1)
# transposed version
sig_pmt = np.transpose(sig_pm)
nu_pmt = np.transpose(nu_pm)

# kernel numerator
K_Bnum = np.zeros((sbn,sbn))
K_Bnum[im] = 4.0*np.pi*(radim[im]+radimt[im])*(Dpm[im]+Dpmt[im])
# left term kernel denominator
K_Blden = np.zeros((sbn,sbn))
K_Blden[im] = (radim[im]+radimt[im])/((radim[im]+radimt[im])+
			((sig_pm[im]**2.0+sig_pmt[im]**2.0)**0.5))
# right term kernel denominator
K_Brden = np.zeros((sbn,sbn))
K_Brden[im] = (4.0*(Dpm[im]+Dpmt[im]))/(((nu_pm[im]**2.0+nu_pmt[im]**2.0)**0.5)*(radim[im]+radimt[im]))
# coagulation kernel (15.33) (m3/(particle s))	
K_B[im] = K_Bnum[im]/(K_Blden[im]+K_Brden[im])


i = Kni>10.0 # free molecular regime
# repeat over size bins
im = np.repeat(i,sbn,1)
# coagulation kernel (15.31) (m3/(particle s))
K_B[im] = np.pi*((radim[im]+radimt[im])**2.0)*((nu_pm[im]+nu_pmt[im])**0.5)


# -------------------------------------------------------------	
# Convective Brownian Diffusion Enhancement kernel:

# Reynold numbers and terminal fall speeds (m/s)
[Rei, Vf] = Reyn_numb(radi, eta_a, rho_a, kin_visc, rho_p, Kni)

# particle Schmidt number (dimensionless) (15.36)
Scp = kin_visc/Dp
Scpm = np.repeat(Scp,sbn,1) # repeat over size bins
# index of particles to coagulate to with Re<=1 and with Re>1
i = Rei<=1.0
i2 = Rei>1.0
# repeat over size bins
im = np.repeat(i,sbn,1)
i2m = np.repeat(i2,sbn,1)
Reim = np.repeat(Rei,sbn,1)
# transpose
Reimt =  np.transpose(Reim)
# empty results matrix
K_DE = np.zeros((sbn,sbn))
# convective Brownian diffusion enhancement kernel (15.35) (m3/(particle s))
K_DE[im] = (K_B[im]*0.45*(Reimt[im]**(1.0/3.0))*(Scpm[im]**(1.0/3.0)))
K_DE[i2m] = (K_B[i2m]*0.45*(Reimt[i2m]**(0.5))*(Scpm[i2m]**(1.0/3.0)))
print Reimt[0,:]

# --------------------------------------------------------------
# Gravitational Collection kernel (could not get agreement with
# fig. 15.7 as of 21/09/2017):

# indices where rj>=ri
im = radimt>=radim
Vfm = np.repeat(Vf,sbn,1)
Vfmt = np.transpose(Vfm)
# dimensionless Stokes number (in text just below eq. 15.39) (needs 
# to be for rj>=ri)
Stm = np.zeros((sbn,sbn))
Stm[im] = (Vfm[im]*(np.abs(Vfmt[im]-Vfm[im])))/(radim[im]*si.g)

# alternative eq from eq 2.32 of Ludlum (1980 (Clouds and Storms))
#St = (2.0*rho_p*(radi**2.0)*(Vf-Vf[0]))/(9.0*eta_a*radj)

# eq. 15.39
EA = np.zeros((sbn,sbn))
EA[im] = (Stm[im]**(2.0))/((Stm[im]+0.5)**2.0)	
# EV matrix for eq. 15.38
EV = np.zeros((sbn,sbn))
# eq. 15.38
im = (Stm>1.214)
EV[im] = (1.0+(0.75*np.log(2.0*Stm[im]))/(Stm[im]-1.214))**(-2.0)
im = (Stm<=1.214)
EV[im] = 0	
# collision efficiency (15.38)
Ecoll = (60.0*EV+EA*Reim)/(60.0+Reim)
im = radimt>=radim
# Gravitational collection kernel (15.37) (m3/(particle s))
K_GC = np.zeros((sbn,sbn))
K_GC[im] = (Ecoll[im]*np.pi*((radim[im]+radimt[im])**2.0)*
		(np.abs(Vfm[im]-Vfmt[im])))	


# ---------------------------------------------------------
# Turbulent Inertial Motion (15.40):
epsilon = 5.0e-4

# note kinematic viscosity is a scalar as it applies to 
# the air, not particles
K_TI = ((np.pi*epsilon**(3.0/4.0))/(si.g*
	kin_visc**(1.0/4.0)))*(radim+radimt)**2.0*(np.abs(Vfm-Vfmt))

# Turbulent Shear (15.41)
K_TS = ((8.0*np.pi*epsilon)/(15.0*kin_visc))**0.5*((radim+radimt)**3.0)

# plot kernels (*1e6 to transform to cm3/(particle s) from 
# m3/(particle s)) against radius of j
f, axarr = plt.subplots(2)
axarr[0].plot(radi,K_B[0,:]*1.0e6,'-b')
axarr[0].plot(radi,K_DE[0,:]*1.0e6,'-g')
axarr[0].plot(radi,K_GC[0,:]*1.0e6,'-r')
axarr[0].plot(radi,K_TI[0,:]*1.0e6,'-m')
axarr[0].plot(radi,K_TS[0,:]*1.0e6,'-k')

axarr[1].plot(radi,K_B[:,49]*1.0e6,'-b')
axarr[1].plot(radi,K_DE[:,49]*1.0e6,'-g')
axarr[1].plot(radi,K_GC[:,49]*1.0e6,'-r')
axarr[1].plot(radi,K_TI[:,49]*1.0e6,'-m')
axarr[1].plot(radi,K_TS[:,49]*1.0e6,'-k')

axarr[0].set_xscale('log')
axarr[0].set_yscale('log')
axarr[0].set_ylim((1.0e-17,1.0e-5))
axarr[1].set_xscale('log')
axarr[1].set_yscale('log')
axarr[1].set_ylim((1.0e-17,1.0e-5))
plt.show()
# ---------------------------------------------------------
# van der Waals and viscous forces p512 - try to replicate
# fig. 15.8:

radi = np.logspace(-11,-5) # radii array (m)
# factor results array and particle pair 
# Knudsen number array: radii of particle i
# in 1st dim., radii of particle j in 2nd
res_all = np.zeros((radi.shape[0],3))
res_Knp = np.zeros((radi.shape[0],3))
# ratios of second particle radius to first 
# particle
jrat = np.array(([1.0,5.0,50.0]))
T = 298.15 # temperature (K)
RH = 0.5 # relative humidity (fraction)
A_H = 200.0*(si.k*1.0e3*T) # Hamaker constant

# j particle radius loop
for j in range(0,jrat.shape[0]): 
	# i particle radius loop
	for i in range(0,radi.shape[0]):

		ri = radi[i] # radius of i (m)
		rj = ri*jrat[j] # radius of j (m)
		# get the correction factors
		[W_c, W_k] = vanW_fact(A_H, T, ri, rj)
		# -----------------------------------------
		# get the particle pair Knudsen number (15.47):

		[Kni, eta_ai, rho_ai, kin_visci] = Knud_calc(RH, T, ri)
		[Knj, eta_aj, rho_aj, kin_viscj] = Knud_calc(RH, T, rj)
		Gi = 1.0+Kni*(1.249+0.42*(np.exp(-0.87/Kni)))
		Gj = 1.0+Knj*(1.249+0.42*(np.exp(-0.87/Knj)))
		# particle diffusion coefficient (15.29) (m2/s) (note the 
		# Boltzmann constant has units (kg m2)/(s2 K), so *1e3 to
		# convert to g from kg)
		Dpi = (((si.k*1.0e3)*T)/(6.0*np.pi*ri*eta_ai))*Gi
		Dpj = (((si.k*1.0e3)*T)/(6.0*np.pi*rj*eta_aj))*Gj	
		Mi = ((4.0/3.0)*np.pi*ri**3.0)*1.0e6
		Mj = ((4.0/3.0)*np.pi*rj**3.0)*1.0e6
		vbari = ((8.0*si.k*1.0e3*T)/(np.pi*Mi))**0.5 #15.32
		vbarj = ((8.0*si.k*1.0e3*T)/(np.pi*Mj))**0.5 #15.32
		# mean free path (15.34)
		lami = (8.0*Dpi)/(np.pi*vbari)
		lamj = (8.0*Dpj)/(np.pi*vbarj)

		Knp = ((lami**2.0+lamj**2.0)**0.5)/(ri+rj)

		# -----------------------------------------
		# get the van der Waals/viscous collision correction 
		# factor (15.42):		
		fac = 4.0*(Dpi+Dpj)/(((vbari**2.0+vbarj**2.0)**0.5)*(ri+rj))
		V_E = (W_c*(1.0+fac))/(1.0+(W_c/W_k)*fac)
		
		# fill in results
		res_Knp[i,j] = Knp
		res_all[i,j] = V_E
		
#axarr[1].plot(res_Knp[:,0],res_all[:,0],'--r')
#axarr[1].plot(res_Knp[:,1],res_all[:,1],'--g')
#axarr[1].plot(res_Knp[:,2],res_all[:,2],'--b')
#axarr[1].set_xscale('log')

#plt.show()



