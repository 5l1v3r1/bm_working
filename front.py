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

# ---------------------------------------------------------
# code to initialse the box model, eq numbers refer to
# Jacobson 2005

import numpy as np
from matplotlib import pyplot as plt


#----------------------------------------------------------
#basic information:

# number of size bins
sbn = np.int(10.0)
# lower and upper size bin bound radii (m)
lsbr = 1.0e-8
usbr = 1.0e-5
# select size bin distribution type:
# 1 for volume-ratio size distribution (p450)
sbdt = 1

# number of components
cn = 3
# component densities (g/m3)
rho = np.array([1.0e6,1.0e6,1.0e6])
# mass fractions per size bin (1st dim.) per component 
# (2nd dim.)
mf = np.zeros((sbn,cn))
mf[:,0] = 0.2
mf[:,1] = 0.3
mf[:,2] = 0.5

# select number size ditribution form:
# 1 for normal (not lognormal)
nsdf=1
if nsdf == 1:
	# mean for probability distribution function
	mu = (sbn-1.0)/2.0 
	# standard deviation for distribution function
	sig = mu/2.0
	# multiplication factor for distribution function 
	# (integral beneath probability distribution 
	# curve is 1)
	fac = 100.0

# ---------------------------------------------------------
# distribute volumes per size bin geometrically over the 
# size range of interest, volume-ratio size distribution:

# volume constant
Vcon = (4.0/3.0)*np.pi
# upper and lower size bin bounds volumes (m3)
usbV = Vcon*(usbr**3.0)
lsbV = Vcon*(lsbr**3.0)

if sbdt == 1: # volume-ratio size distribution (p450)
	# volume ratio between neighbouring larger bin and 
	# smaller bin (13.3)
	Vrat = (usbV/lsbV)**(1.0/(sbn))
	# array of size bin bound volumes (m)
	Vbou = np.arange(0,sbn+1)
	Vbou = lsbV*(Vrat**Vbou)

# volume at size bin centre (m3) (13.5):
Varr = 0.5*(Vbou[0:sbn]+Vbou[1::])


# ---------------------------------------------------------
# populate size bins with number concentration of 
# particles (#/m3(air)):

if nsdf == 1:
	 
	# base for distribution function
	x = np.arange(0,sbn) 	
	# number concentration per size bin (#/m3(air))
	num = ((1.0/((2.0*np.pi*sig**2.0)**0.5))*
		np.exp(-((x-mu)**2.0)/(2.0*sig**2.0)))*fac


# ---------------------------------------------------------
# mass concentration of components per size bin (g/m3(air))

# average density (g/m3(particle)) per size bin
rhobar = np.sum(mf*rho,1)
# total volume of particles per size bin 
# (m3(particle)/m3(air))
VT = Varr*num
# total mass of particles per size bin (g/m3(air))
mT = np.zeros((sbn,1))
mT[:,0] = VT*rhobar
# mass per size bin (1st dim.) per component (2nd dim.)
# (g/m3(air))
m0 = mT*mf

# verify above calculations:
#plt.plot(Varr,'-xg')
#plt.plot(Vbou,'--or')
#plt.yscale('log')
#plt.show()
