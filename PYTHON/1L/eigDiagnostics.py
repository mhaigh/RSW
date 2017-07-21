# eigDiagnostics.py
# File containing functions to be called by the master script EIG.py.
#====================================================

import sys

import numpy as np
import matplotlib.pyplot as plt

from diagnostics import diff
from diagnostics import extend

#====================================================

# The potential vorticity (both of background state and normal modes)
# PV
def PV(u,v,eta,u_full,eta_full,H_nd,U_nd,f_nd,dx_nd,dy_nd,N):
# Takes input u, v, eta which are the eigenvectors of the 1-L SW system, projected onto two dimensions, 
# i.e. u = u_vec(y) * exp(i*k*x). 

	# Relative vorticity
	RV_full = diff(v,1,1,dx_nd) - diff(u_full,0,0,dy_nd);
	RV_prime = diff(v,1,1,dx_nd) - diff(u,0,0,dy_nd);
	RV_BG = - diff(U_nd,2,0,dy_nd);

	# Potential vorticity
	PV_BG = (f_nd + RV_BG) / H_nd; 
	PV_full = np.zeros((N,N));
	PV_prime = np.zeros((N,N));
	for i in range(0,N):
		PV_full[:,i] = (f_nd[:] + RV_full[:,i]) / eta_full[:,i];
		PV_prime[:,i] = PV_full[:,i] - PV_BG[:];

	return PV_full, PV_prime;

#====================================================

# footprint
# A function that calculates the PV footprint of the 1L SW solution as produced by RSW_1L.py
def footprint(u_full,v,PV_full,x_nd,dx_nd,dy_nd,N):
# This code calculates the PV footprint, i.e. the PV flux convergence defined by
# P = -(div(u*q,v*q)-div((u*q)_av,(v*q)_av)), where _av denotees the time average.
# We will take the time average over one whole forcing period, and subtract this PV flux
# convergence from the PV flux convergence at the end of one period.
# To calculate footprints, we use the full PV and velocity profiles.


	qu = PV_full * u_full;		# Zonal PV flux
	qv = PV_full * v;		# Meridional PV flux

	# Next step: taking appropriate derivatives of the fluxes. To save space, we won't define new variables, but instead overwrite the old ones.
	# From these derivatives we can calculate the 
	P = - diff(qu,1,1,dx_nd) - diff(qv,0,0,dy_nd);		
	
	# We are interested in the zonal average of the footprint
	P_xav = np.trapz(P,x_nd[:N],dx_nd,axis=1);
	
	P[0,:] = 0; P[N-1,:] = 0;

	P = extend(P);

	return P, P_xav


#====================================================

# eigPlot
def eigPlots(u_proj,v_proj,eta_proj,u_nd,v_nd,eta_nd,x_nd,y_nd,sol):
	
	u_proj = extend(u_proj);
	v_proj = extend(v_proj);
	eta_proj = extend(eta_proj);

	if sol:

		u_nd = extend(u_nd);
		v_nd = extend(v_nd);
		eta_nd = extend(eta_nd);

		plt.figure(1,figsize=[21,6]);
		plt.subplot(231)
		plt.contourf(x_nd,y_nd,u_proj);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.subplot(232)
		plt.contourf(x_nd,y_nd,v_proj);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.subplot(233)
		plt.contourf(x_nd,y_nd,eta_proj);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
	
		plt.subplot(234)
		plt.contourf(x_nd,y_nd,u_nd);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.subplot(235)
		plt.contourf(x_nd,y_nd,v_nd);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.subplot(236)
		plt.contourf(x_nd,y_nd,eta_nd);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.tight_layout();
		plt.show();

	else:	

		plt.figure(1,figsize=[21,6]);
		plt.subplot(131)
		plt.contourf(x_nd,y_nd,u_proj);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.subplot(132)
		plt.contourf(x_nd,y_nd,v_proj);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.subplot(133)
		plt.contourf(x_nd,y_nd,eta_proj);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.tight_layout();
		plt.show();

#====================================================

# orderEigenmodes
def orderEigenmodes(vec,val,N,VECS):
# A function that takes the set of eigenmodes, given by vec, and orders them according to the number of zero crossings.
# When two or more eigenmodes cross zeros the same amount of times, they are ordered by their frequency, smallest first.
	
	dim = np.size(val);

	if VECS:
		vec = np.array(vec);
		u_vec = vec[0,:,:];
		v_vec = vec[1,:,:];
		eta_vec = vec[2,:,:];

	else:
		u_vec = vec[0:N,:];		# Extract the eigenmodes.
		v_vec = vec[N:2*N,:];
		eta_vec = vec[2*N:3*N,:];		

	# Initialise a counter for the number of zero crossings. 
	count = np.zeros((dim),dtype=int);
	for wi in range(0,dim):
		for j in range(1,N):
			if (u_vec[j-1,wi] >= 0 and u_vec[j,wi] <= 0) or (u_vec[j-1,wi] <= 0 and u_vec[j,wi] >= 0):
				count[wi] = count[wi] + 1;

	i_count = np.argsort(count);

	return count, i_count;

#====================================================

# vec2vecs
def vec2vecs(vec,N,dim,BC):
# A function to take the full eigenvector, vec = (u,v,eta), and 
# extract u, v, eta, depending on the boundary condition.

	if BC == 'FREE-SLIP':
		u_vec = vec[0:N,:];
		v_vec = np.zeros((N,dim),dtype=complex);
		v_vec[1:N-1,:] = vec[N:2*N-2,:];
		eta_vec = vec[2*N-2:3*N-2,:];
	elif BC == 'NO-SLIP':
		u_vec = np.zeros((N,dim),dtype=complex);
		u_vec[1:N-1,:] = vec[0:N-2,:];
		v_vec = np.zeros((N,dim),dtype=complex);
		v_vec[1:N-1,:] = vec[N:2*N-2,:];
		eta_vec = vec[2*N-4:3*N-4,:];	
		return val, u_vec, v_vec, eta_vec;		
	else:
		sys.exit('ERROR: choose valid BC');

	
	return u_vec, v_vec, eta_vec;


	
	
