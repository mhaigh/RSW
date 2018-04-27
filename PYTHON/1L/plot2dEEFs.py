# plot2dEEFs.py

#========================================================

import numpy as np
import matplotlib.pyplot as plt

from core import BG_state

from inputFile import *

#========================================================

Ns = 151

sigma_set = np.linspace(0.015,0.045,Ns) * 3840000.0
ds = sigma_set[1] - sigma_set[0]

y0_min = y[0] + L/3;					# We want to keep the forcing at least one gridpoint away from the boundary
y0_max = y[N-1] - L/3;
y0_set = [];							# Initialise an empty set of forcing latitudes
y0_index_set = [];
for j in range(0,N):
	if y0_min <= y[j] <= y0_max:
		y0_set.append(y[j]);			# Build the set of forcing locations, all at least 1 gridpoint away from the boundary.	
		y0_index_set.append(j);
y0_set = np.array(y0_set) / L;
y0_index_set = np.array(y0_index_set);
nn = np.shape(y0_set)[0]

sigma_mesh, y0_mesh = np.mgrid[slice(sigma_set.min()/1000.0,(sigma_set.max()+ds)/1000.,ds/1000.),slice(y0_set.min(),y0_set.max()+dy_nd,dy_nd)];

#========================================================


Qy = np.zeros((Ns,nn))

for ui in range(0,Ns):

	sigma = sigma_set[ui]
	U0, H0 = BG_state.BG_Gaussian(Umag,sigma,JET_POS,Hflat,f0,beta,g,y,L,N)
	U0 = U0 / U; H0 = H0 / chi
	U0 = U0[y0_index_set]; H0 = H0[y0_index_set]

	U0y =  diff(U0,2,0,dy_nd)
	Q = (f_nd[y0_index_set] / Ro - U0y) / H0
	Qy[ui,:] = diff(Q,2,0,dy_nd)
	

fs = 18

CS = plt.contour(y0_set,sigma_set/1000.,Qy,1,colors='k')
plt.clabel(CS, fontsize=9, inline=1)
plt.pcolormesh(y0_mesh,sigma_mesh,Qy)

plt.xlabel('y0',fontsize=fs)
plt.ylabel('Jet width (km)',fontsize=fs)
plt.xticks(fontsize=fs-2)
plt.yticks(fontsize=fs-2)
plt.show()


