# EQUIV_EDDY_FLUXES.py
#=======================================================
# This is an executable code that solves the 1L shallow water system a number of times, each time storing the equivalent eddy flux.
#=======================================================

import os
import numpy as np
import matplotlib.pyplot as plt
import time

import diagnostics
import PV
import momentum
import forcing_1L
import solver
#import output

from inputFile import *

#=======================================================

start = time.time();

filename = 'EEF_' + str(int(period_days));
filename_u = 'EEF_u_' + str(int(period_days));
filename_v = 'EEF_v_' + str(int(period_days));

# Can test against U0 or y0, or find the buoyancy vs U0 or y0.
TEST = 'y0';
nn = 5;

if TEST == 'U0':
	U0_set = np.linspace(0.2,0.25,nn);
	U_ones = np.ones(N);
	U0_nd_set = U0_set / U; 
	if FORCE_TYPE == 'CTS':
		F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_cts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
	elif FORCE_TYPE == 'DCTS':
		F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_dcts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
	else:
		sys.exit('ERROR: Invalid forcing option selected.');

if TEST == 'y0':
	y0_min = y[0] + r0;					# We want to keep the forcing at least one gridpoint away from the boundary
	y0_max = y[N-1] - r0;
	y0_set = [];						# Initialise an empty set of forcing latitudes
	y0_index_set = [];
	for j in range(0,N):
		if y0_min <= y[j] <= y0_max:
			y0_set.append(y[j]);		# Build the set of forcing locations, all at least 1 gridpoint away from the boundary.	
			y0_index_set.append(j);
	y0_set = np.array(y0_set);
	y0_index_set = np.array(y0_index_set);
	nn = np.shape(y0_set)[0];
	print(nn);
	a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
			
# Initialise the array which stores the EEF values.
if os.path.isfile(filename + '.npy'):
	EEF_PV = np.load(filename + '.npy');
else:
	if footprintComponents:
		EEF_PV = np.zeros((nn,6,2));
	else:
		EEF_PV = np.zeros((nn,2));

EEF_u = np.zeros((nn,2));
EEF_v = np.zeros((nn,2));

# Now start the loop over each forcing index.
for ii in range(0,nn):
	print(ii);
	
	if EEF_PV[ii,0,0] == 0:

		# If TEST==U0, linear problem has to be redefined each iteration.
		if TEST == 'U0':
			U0_nd = U0_nd_set[ii] * U_ones;	# Redefine U0 in each run.
			a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
				
		# If TEST==y0, matrix only needs to be defined once, but forcing must be defined each iteration.
		if TEST == 'y0':
			y0 = y0_set[ii];				# Redefine y0 and the forcing in each run.
			y0_index = y0_index_set[ii];
			y0_nd = y0 / L;
			# Forcing
			if FORCE_TYPE == 'CTS':
				F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_cts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
			elif FORCE_TYPE == 'DCTS':
				F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_dcts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
			else:
				sys.exit('ERROR: Invalid forcing option selected.');	
	
		# Solver
		if BC == 'NO-SLIP':
			solution = solver.NO_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ro*Ftilde3_nd,N,N2);
		if BC == 'FREE-SLIP':
			solution = solver.FREE_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2);
	
		utilde_nd, vtilde_nd, etatilde_nd = solver.extractSols(solution,N,N2,BC);
		u_nd, v_nd, eta_nd = solver.SPEC_TO_PHYS(utilde_nd,vtilde_nd,etatilde_nd,T_nd,dx_nd,omega_nd,N);
			
		# Take real part.
		u_nd = np.real(u_nd);
		v_nd = np.real(v_nd);
		eta_nd = np.real(eta_nd);
	
		# Normalise all solutions by the (non-dimensional) forcing amplitude. 
		u_nd = u_nd / AmpF_nd;
		v_nd = v_nd / AmpF_nd;
		eta_nd = eta_nd / AmpF_nd;
	
		# In order to calculate the vorticities of the system, we require full (i.e. BG + forced response) u and eta.
		eta_full = np.zeros((N,N,Nt));
		u_full = np.zeros((N,N,Nt));
		for j in range(0,N):
			eta_full[j,:,:] = eta_nd[j,:,:] + H0_nd[j];
			u_full[j,:,:] = u_nd[j,:,:] + U0_nd[j];
	
		# Calculate PV fields and PV fluxes.
		PV_prime, PV_full, PV_BG = PV.potentialVorticity(u_nd,v_nd,eta_nd,u_full,eta_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd);
		uq, Uq, uQ, UQ, vq, vQ = PV.fluxes(u_nd,v_nd,U0_nd,PV_prime,PV_BG,N,Nt);
		P, P_xav = PV.footprint(uq,Uq,uQ,UQ,vq,vQ,x_nd,T_nd,dx_nd,dy_nd,N,Nt);			
		EEF_PV[ii,:] = PV.EEF(P_xav,y_nd,y0_nd,y0_index,dy_nd,omega_nd,N);

		# Calculate momentum fluxes and footprints
		uu, uv, vv = momentum.fluxes(u_nd,v_nd);
		Mu, Mv, Mu_xav, Mv_xav = momentum.footprint(uu,uv,vv,x_nd,T_nd,dx_nd,dy_nd,N,Nt);
		EEF_u[ii,:], EEF_v[ii,:] = momentum.EEF_mom(Mu_xav,Mv_xav,y_nd,y0_nd,y0_index,dy_nd,omega_nd,N);
		
	np.save(filename,EEF_PV);
	np.save(filename_u,EEF_u);
	np.save(filename_v,EEF_v);
	
elapsed = start - time.time();
elapsed = np.ones(1) * elapsed;
print(elapsed);

	
	
	
