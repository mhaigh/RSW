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

from inputFile import *

import time

#=======================================================

start = time.time();

pe = 1;		# Number of processors

filename = 'EEF_PV';

y0_min = y[0] + r0;					# We want to keep the forcing at least one gridpoint away from the boundary
y0_max = y[N-1] - r0;
y0_set = [];						# Initialise an empty set of forcing latitudes
for j in range(0,N):
	if y0_min <= y[j] <= y0_max:
		y0_set.append(j);		# Build the set of forcing locations, all at least 1 gridpoint away from the boundary.	
y0_set = np.array(y0_set);
nn = np.shape(y0_set)[0];
a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);

if TEST == 'U0':
	nn = 3;
	U0_set = np.linspace(-0.3,0.5,nn);
	if FORCE_TYPE == 'CTS':
		F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing.forcing_cts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
	elif FORCE_TYPE == 'DCTS':
		F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing.forcing_dcts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
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
	a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
			
# Now split the input set into pe sample sets.
sets = [];
d = nn / pe;
r = nn % pe;
#print(nn,d,r);
for pe_no in range(1,r+1):
	exec('y0_set_' + str(pe_no) + '= y0_set[(pe_no-1)*d:pe_no*d+1]');
	exec('sets.append(y0_set_'+str(pe_no)+')');
for pe_no in range(r+1,pe+1):
	exec('y0_set_' + str(pe_no) + '= y0_set[(pe_no-1)*d:pe_no*d]');
	exec('sets.append(y0_set_'+str(pe_no)+')');
print(sets);

EEF_array = np.zeros((nn,2));
l_array = np.zeros((nn,2));
PV_xav = np.zeros((nn,N));

#=======================================================

# Now define two EEF functions, one for TEST = U0 and one for TEST = y0.

#=======================================================

def EEF_y0(y0_set,pi):
# For every y0 value in y0_set all forcing parameters must be redefined.

	from inputFile_1L import *
	
	yn = len(y0_set);
	
	EEF_array = np.zeros((yn,6,2));
		
	for yi in range(0,yn):					# yi indexes the local EEF_array (i.e. computational domain)
		ii = y0_set[yi];					# ii indexes arrays defined over global domain
		print(ii);		
		if EEF_array[yi,0,0] == 0:			# Check if any of the array has been updated after initialisation.
	
			y0 = y[ii];						# Redefine y0 and the forcing in each run.
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
			P, P_xav[yi,:] = PV.footprint_1L(u_full,v_nd,eta_full,PV_full,U0_nd,U,Umag,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,AmpF_nd,FORCE,r0,nu,BG,Fpos,ts,period_days,N,Nt,GAUSS);			
			EEF_array[yi,:], l_array[yi,:] = PV.EEF(P_xav[yi,:],y_nd,y0_nd,dy_nd,omega_nd,N);

	filename = 'EEF_array_' + str(pi);
	exec('np.save(filename,EEF_array)');

#=======================================================

def EEF_U0(U0_set,pi):
# For every U0 value in U0_set, 

	yn = len(y0_set);
	
	EEF_array = np.zeros((yn,6,2));
		
	for yi in range(0,yn):					# yi indexes the local EEF_array (i.e. computational domain)
		ii = y0_set[yi];					# ii indexes arrays defined over global domain
		print(ii);		

		# Redefine U0 and H0 in each run.
		for j in range(0,N):
			U0[j] = U0_set[ii];
			H0[j] = - (U0[j] / g) * (f0 * y[j] + beta * y[j]**2 / 2) + Hflat;
		U0_nd = U0 / U;
		H0_nd = H0 / chi; 

		a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
		
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
		P, P_xav[yi,:] = PV.footprint_1L(u_full,v_nd,eta_full,PV_full,U0_nd,U,Umag,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,AmpF_nd,FORCE,r0,nu,BG,Fpos,ts,period_days,N,Nt,GAUSS);			
		EEF_array[yi,:], l_array[yi,:] = PV.EEF(P_xav[yi,:],y_nd,y0_nd,dy_nd,omega_nd,N);

	filename = 'EEF_array_' + str(pi);
	exec('np.save(filename,EEF_array)');


if __name__ == '__main__':
	jobs = [];
	for pi in range(0,pe):
		print(pi);
		if TEST == 'y0':
			p = mp.Process(target=EEF_y0,args=(sets[pi],pi));
		elif TEST == 'U0':
			p = mp.Process(target=EEF_U0,args=(sets[pi],pi));
		jobs.append(p);
		p.start();

	for p in jobs:
		p.join();

# Now collect the results by reloading them, and compile into one array
EEF_array = np.zeros((nn,6,2));
yn_count = 0;
for pi in range(0,pe):
	filename = 'EEF_array_' + str(pi) + '.npy';
	exec('EEF_array_tmp = np.load(filename)');
	yn = np.shape(EEF_array_tmp)[0];
	EEF_array[yn_count:yn_count+yn,:,:] = EEF_array_tmp[:,:,:];
	yn_count = yn_count + yn;
	os.remove(filename);
	
np.save('EEF_array',EEF_array);

elapsed = (time.time() - start);
elapsed = np.ones(1) * elapsed;
print(elapsed);
	
	
