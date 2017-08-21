# EQUIV_EDDY_FLUXES.py
#=======================================================
# This is an executable code that solves the 1L shallow water system a number of times, each time storing the equivalent eddy flux.
#=======================================================

import numpy as np
import matplotlib.pyplot as plt

import diagnostics
import PV
import forcing_1L
import solver
#import output

from inputFile_1L import *

#=======================================================

# Can test against U0 or y0, or find the buoyancy vs U0 or y0.
TEST = 'y0';
nn = 5;

if TEST == 'U0':
	U0_set = np.linspace(0.2,0.25,nn);
	U_ones = np.ones(N);
	U0_nd_set = U0_set / U; 
	if FORCE_TYPE == 'CTS':
		F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_cts(x,y,K,y0,r0,N,FORCE,AmpF,g,f,f0,U,L,dx,dy);
	elif FORCE_TYPE == 'DCTS':
		F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_dcts(x,y,K,y0,r0,N,FORCE,AmpF,g,f,f0,U,L,dx,dy);
	else:
		sys.exit('ERROR: Invalid forcing option selected.');

if TEST == 'y0':
	y0_min = y[0] + r0;					# We want to keep the forcing at least one gridpoint away from the boundary
	y0_max = y[N-1] - r0;
	y0_set = [];						# Initialise an empty set of forcing latitudes
	for j in range(0,N):
		if y0_min <= y[j] <= y0_max:
			y0_set.append(y[j]);		# Build the set of forcing locations, all at least 1 gridpoint away from the boundary.	
	y0_set = np.array(y0_set);
	nn = np.shape(y0_set)[0];
	print(nn);
	a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
			
# Initialise the array which stores the EEF values.
if footprintComponents:
	EEF_array = np.zeros((nn,6,2));
else:
	EEF_array = np.zeros((nn,2));
# Now start the loop over each forcing index.
for ii in range(0,nn):
	print(ii);

	# If TEST==U0, linear problem has to be redefined each iteration.
	if TEST == 'U0':
		U0_nd = U0_nd_set[ii] * U_ones;	# Redefine U0 in each run.
		a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
			
	# If TEST==y0, matrix only needs to be defined once, but forcing must be defined each iteration.
	if TEST == 'y0':
		y0 = y0_set[ii];				# Redefine y0 and the forcing in each run.
		y0_nd = y0 / L;
		# Forcing
		if FORCE_TYPE == 'CTS':
			F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_cts(x,y,K,y0,r0,N,FORCE,AmpF,g,f,f0,U,L,dx,dy);
		elif FORCE_TYPE == 'DCTS':
			F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_dcts(x,y,K,y0,r0,N,FORCE,AmpF,g,f,f0,U,L,dx,dy);
		else:
			sys.exit('ERROR: Invalid forcing option selected.');	

	# Solver
	if BC == 'NO-SLIP':
		solution = solver.NO_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N,N2);
	if BC == 'FREE-SLIP':
		solution = solver.FREE_SLIP_SOLVER2(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N,N2);

	utilde_nd, vtilde_nd, etatilde_nd = solver.extractSols(solution,N,N2,BC);
	u_nd, v_nd, eta_nd = solver.SPEC_TO_PHYS(utilde_nd,vtilde_nd,etatilde_nd,T_nd,dx_nd,omega_nd,N);
		
	u_nd = np.real(u_nd);
	v_nd = np.real(v_nd);
	eta_nd = np.real(eta_nd);

	# In order to calculate the vorticities of the system, we require full (i.e. BG + forced response) u and eta.
	eta_full = np.zeros((N,N,Nt));
	u_full = np.zeros((N,N,Nt));
	for j in range(0,N):
		eta_full[j,:,:] = eta_nd[j,:,:] + H0_nd[j];
		u_full[j,:,:] = u_nd[j,:,:] + U0_nd[j];

	# Calculate PV fields.
	PV_prime, PV_full, PV_BG = PV.vort(u_nd,v_nd,eta_nd,u_full,eta_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd);

	# Do footprints
	if footprintComponents:
		P, uq, uQ, Uq, UQ, vq, vQ, P_xav, uq_xav, uQ_xav, Uq_xav, vq_xav, vQ_xav = PV.footprintComponents(u_nd,v_nd,eta_nd,PV_prime,PV_BG,U0_nd,AmpF_nd,x_nd,dx_nd,dy_nd,N,Nt);
		EEF_array[ii,:,:] = PV.EEF_components(P_xav,uq_xav,uQ_xav,Uq_xav,vq_xav,vQ_xav,y_nd,y0_nd,dy_nd,omega_nd,N);
	else: 
		P, P_xav = PV.footprint_1L(u_full,v_nd,eta_full,PV_full,U0_nd,U,Umag,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,AmpF_nd,FORCE,r0,nu,BG,Fpos,ts,period_days,N,Nt,GAUSS);			
		EEF_array[ii,:] = PV.EEF(P_xav,y_nd,y0_nd,dy_nd,omega_nd,N);

if footprintComponents:
	#output.ncSaveEEF_y0_components(EEF_array,y0_set,period_days,nn);
	filename = 'EEF_' + str(int(period_days));
	np.save(filename,EEF_array);
else:
	# not written this nc function yet.
	a=1;




	
	
	
