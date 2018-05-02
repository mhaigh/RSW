# EQUIV_EDDY_FLUXES.py
#=======================================================
# This is an executable code that solves the 1L shallow water system a number of times, each time storing the equivalent eddy flux.
#=======================================================

import os
import numpy as np
import time

from core import solver, PV, momentum, diagnostics, BG_state

from inputFile import *

#=======================================================

start = time.time();

filename = 'EEF_PV';

# Can test against U0 or y0, or find the buoyancy vs U0 or y0.

#=======================================================

pe = 11		# Number of processors

# Initialise tests

Nu = 50;
test_set = np.linspace(0.015,0.045,Nu) * 3840000.0

y0_min = y[0] + L/3;					# We want to keep the forcing at least one gridpoint away from the boundary
y0_max = y[N-1] - L/3;
y0_set = [];						# Initialise an empty set of forcing latitudes
y0_index_set = [];
for j in range(0,N):
	if y0_min <= y[j] <= y0_max:
		y0_set.append(y[j]);		# Build the set of forcing locations, all at least 1 gridpoint away from the boundary.	
		y0_index_set.append(j);
y0_set = np.array(y0_set);
y0_index_set = np.array(y0_index_set);
nn = np.shape(y0_set)[0];
#print(nn)

# Now split the input set into pe sample sets.
sets = [];
d = Nu // pe + 1	  	# Number of elements in each set (apart from the last one).
r = Nu % d				# Number of elements left over, for the last processor.
#print(nn,d,r);
for pe_no in range(1,pe):
	exec('test_set_' + str(pe_no) + '= test_set[(pe_no-1)*d:pe_no*d]');
	exec('sets.append(test_set_'+str(pe_no)+')');
for pe_no in range(pe,pe+1):
	exec('test_set_' + str(pe_no) + '= test_set[(pe_no-1)*d:pe_no*d+1]');
	exec('sets.append(test_set_'+str(pe_no)+')');
print(sets)

#=======================================================

def EEF_main(test_set,pi)

	NU = len(set_);
	
	EEF_array = np.zeros((NU,nn,2));
	
	# Now start the loop over each forcing index.
	for ui in range(0,NU):

		# Redefine U0 and H0.
		#sigma = test_set[ui]
		Umag = test_set[ui]
		U0, H0 = BG_state.BG_Gaussian(Umag,sigma,JET_POS,Hflat,f0,beta,g,y,L,N)
		U0_nd = U0 / U;
		H0_nd = H0 / chi; 
		a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
	
		for yi in range(0,nn):

			y0 = y0_set[yi];				# Redefine y0 and the forcing in each run.
			y0_index = y0_index_set[yi];
			y0_nd = y0 / L;
			# Forcing
			F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing.forcing_cts2(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,bh,dx_nd,dy_nd)
	
			solution = solver.FREE_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2)
	
			utilde_nd, vtilde_nd, etatilde_nd = solver.extractSols(solution,N,N2,BC)
			u, v, h = solver.SPEC_TO_PHYS(utilde_nd,vtilde_nd,etatilde_nd,T_nd,dx_nd,omega_nd,N)
		
			# Take real part.
			u = np.real(u)
			v = np.real(v)
			h = np.real(h)
	
			# Normalise all solutions by the (non-dimensional) forcing amplitude. 
			u = u / AmpF_nd
			v = v / AmpF_nd
			h = h / AmpF_nd
		
			# In order to calculate the vorticities of the system, we require full (i.e. BG + forced response) u and eta.
			h_full = np.zeros((N,N,Nt));
			u_full = np.zeros((N,N,Nt));
			for j in range(0,N):
				h_full[j,:,:] = h[j,:,:] + H0_nd[j];
				u_full[j,:,:] = u[j,:,:] + U0_nd[j];
	
			# Calculate PV fields and PV fluxes.
			PV_prime, PV_full, PV_BG = PV.potentialVorticity(u,v,h,u_full,h_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd,Ro)
			uq, Uq, uQ, UQ, vq, vQ = PV.fluxes(u,v,U0_nd,PV_prime,PV_BG,N,Nt)
			P, P_xav = PV.footprint(uq,Uq,uQ,UQ,vq,vQ,x_nd,T_nd,dx_nd,dy_nd,N,Nt)			
			EEF_PV[ui,yi,:], l_PV = PV.EEF(P_xav,y_nd,y0_nd,y0_index,dy_nd,N)


if __name__ == '__main__':
	jobs = [];
	for pi in range(0,pe):
		print(pi);
		if TEST == 'y0':
			p = mp.Process(target=EEF_main,args=(sets[pi],pi));
		elif TEST == 'U0':
			p = mp.Process(target=EEF_main,args=(sets[pi],pi));
		jobs.append(p);
		p.start();

	for p in jobs:
		p.join();

# Now collect the results by reloading them, and compile into one array
EEF_array = np.zeros((Nn,));
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
	
elapsed = time.time() - start;
elapsed = np.ones(1) * elapsed;
print(elapsed);

	
	
	
