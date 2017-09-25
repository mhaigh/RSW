# EQUIV_EDDY_FLUXES.py
#=======================================================
# This is an executable code that solves the 1L shallow water system a number of times, each time storing the equivalent eddy flux.
#=======================================================

import os
import sys

import numpy as np
import multiprocessing as mp

import diagnostics
import PV
import forcing_1L
import solver

from inputFile_1L import *

#=======================================================

pe = 1;		# Number of processors

filename = 'EEF__' + str(int(period_days));

y0_min = y[0] + r0;					# We want to keep the forcing at least one gridpoint away from the boundary
y0_max = y[N-1] - r0;
y0_set = [];						# Initialise an empty set of forcing latitudes
for j in range(0,N):
	if y0_min <= y[j] <= y0_max:
		y0_set.append(j);		# Build the set of forcing locations, all at least 1 gridpoint away from the boundary.	
y0_set = np.array(y0_set);
nn = np.shape(y0_set)[0];
shift = y0_set[0]
a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
			
# Now split the input set into pe sample sets.
y0_sets = [];
d = nn / pe;
r = nn % pe;
print(nn,d,r);
for pe_no in range(1,r+1):
	exec('y0_set_' + str(pe_no) + '= y0_set[(pe_no-1)*d:pe_no*d+1]');
	exec('y0_sets.append(y0_set_'+str(pe_no)+')');
for pe_no in range(r+1,pe+1):
	exec('y0_set_' + str(pe_no) + '= y0_set[(pe_no-1)*d:pe_no*d]');
	exec('y0_sets.append(y0_set_'+str(pe_no)+')');
print(y0_sets);

# Initialise the array which stores the EEF values.
if os.path.isfile(filename + '.npy'):
	EEF_array = np.load(filename + '.npy');
else:
	if footprintComponents:
		EEF_array = np.zeros((nn,6,2));
	else:
		EEF_array = np.zeros((nn,2));

def EEF_y0(y0_set,EEF_array):
	from inputFile_1L import *
	
	yn = len(y0_set);
		
	for yi in range(0,yn):	
		ii = y0_set[yi];	
		print(ii);		
		if EEF_array[ii-shift,0,0] == 0:			# Check if any of the array has been updated after initialisation.
	
			y0 = y[ii];						# Redefine y0 and the forcing in each run.
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
		
			# Do footprints
			if footprintComponents:
				P, P_uq, P_uQ, P_Uq, P_vq, P_vQ, P_xav, P_uq_xav, P_uQ_xav, P_Uq_xav, P_vq_xav, P_vQ_xav = PV.footprintComponents(uq,Uq,uQ,vq,vQ,x_nd,T_nd,dx_nd,dy_nd,N,Nt);
				EEF_array[ii-shift,:,:] = PV.EEF_components(P_xav,P_uq_xav,P_uQ_xav,P_Uq_xav,P_vq_xav,P_vQ_xav,y_nd,y0_nd,dy_nd,omega_nd,N);
			else: 
				P, P_xav = PV.footprint_1L(u_full,v_nd,eta_full,PV_full,U0_nd,U,Umag,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,AmpF_nd,FORCE,r0,nu,BG,Fpos,ts,period_days,N,Nt,GAUSS);			
				EEF_array[ii-shift,:] = PV.EEF(P_xav,y_nd,y0_nd,dy_nd,omega_nd,N);

	np.save(filename,EEF_array);

if __name__ == '__main__':
	jobs = [];
	for pi in range(0,pe):
		print(pi);
		p = mp.Process(target=EEF_y0,args=(y0_sets[pi],EEF_array))
		jobs.append(p);
		p.start();

	for p in jobs:
		p.join();


	
	
	
