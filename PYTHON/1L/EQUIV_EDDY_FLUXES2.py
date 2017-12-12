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
import thickness
#import output

from inputFile import *

#=======================================================

start = time.time();

filename = 'EEF_PV';

# Can test against U0 or y0, or find the buoyancy vs U0 or y0.
TEST = 'U0';

if TEST == 'U0':
	nn = 10;
	U0_set = np.linspace(-0.03,-0.02,nn);
	if FORCE_TYPE == 'CTS':
		F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_cts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
	elif FORCE_TYPE == 'DCTS':
		F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_dcts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
	else:
		sys.exit('ERROR: Invalid forcing option selected.');

EEF_PV = np.zeros((nn,2));
P_xav = np.zeros((nn,N))

# Now start the loop over each forcing index.
for ii in range(0,nn):
	print(ii);

	# If TEST==U0, linear problem has to be redefined each iteration.
	if TEST == 'U0':
		# Redefine U0 and H0 in each run.
		for j in range(0,N):
			U0[j] = U0_set[ii];
			H0[j] = - (U0[j] / g) * (f0 * y[j] + beta * y[j]**2 / 2) + Hflat;
		U0_nd = U0 / U;
		H0_nd = H0 / chi; 

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
		solution = solver.FREE_SLIP_SOLVER4(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2);
	
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
	PV_prime1, PV_prime2, PV_prime3 = PV.potentialVorticity_linear(u_nd,v_nd,eta_nd,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd,Ro);
	vq2 = v_nd * PV_prime2;
	vq2 = diagnostics.timeAverage(vq2,T_nd,Nt);
	vq2_y = - diagnostics.diff(vq2,0,0,dy_nd);
	vq2_y = diagnostics.extend(vq2_y);
	P_xav[ii,:] = np.trapz(vq2_y,x_nd,dx_nd,axis=1);
	EEF_PV[ii,:] = PV.EEF(P_xav[ii,:],y_nd,y0_nd,y0_index,dy_nd,N)

#np.save(filename,EEF_PV);

plt.contourf(P_xav);
plt.show();

plt.plot(EEF_PV);
plt.show();

elapsed = time.time() - start;
elapsed = np.ones(1) * elapsed;
print(elapsed);

	
	
	