# EQUIV_EDDY_FLUXES.py
#=======================================================
# This is an executable code that solves the 1L shallow water system a number of times, each time storing the equivalent eddy flux.
#=======================================================

import os
import numpy as np
import multiprocessing as mp
import time
from scipy.ndimage.measurements import center_of_mass

from core import solver, PV, momentum, diagnostics, BG_state

from inputFile import *

#=======================================================

start = time.time();

# Can test against U0 or y0, or find the buoyancy vs U0 or y0.

#=======================================================

pe = 1		# Number of processors

# Initialise tests

Nu = 3;

test_set = np.array((50.,60.,70.))
#test_set = np.array((50.,100.,200.))
#test_set = np.array((60000.,90000.,12000.)) # r0
#test_set = np.linspace(0.015,0.045,Nu) * 3840000.0
#test_set = np.linspace(0.4,1.2,Nu)

y0_min = y[0] + r0# + L/3 					# We want to keep the forcing at least one gridpoint away from the boundary
y0_max = y[N-1] - r0# - L/3
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

def EEF_main(set_,pi):

	NU = len(set_);
	
	EEF_array = np.zeros((NU,nn,2));
	com = np.zeros((NU,nn));
	
	# Now start the loop over each forcing index.
	for ui in range(0,NU):

		# Redefine U0 and H0.
		#sigma = set_[ui]
		#Umag = set_[ui]
		#U0, H0 = BG_state.BG_Gaussian(Umag,sigma,JET_POS,Hflat,f0,beta,g,y,L,N)
		#U0_nd = U0 / U;
		#H0_nd = H0 / chi; 

		# r0
		#r0 = set_[ui]
		#r0_nd = r0 / L

		# period
		period_days = set_[ui]
		period = 3600. * 24. * period_days;
		omega = 1. / (period);          		
		T = np.linspace(0,period,Nt+1);		
		dt = T[1] - T[0];						
		t = T[ts];		

		omega_nd = omega * T_adv;      	
		t_nd = t / T_adv;
		T_nd = T / T_adv;
		dt_nd = dt / T_adv;

		# k
		#nu = set_[ui]
		#Re = L * U / nu	

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

			com[ui,yi] = center_of_mass(np.abs(P_xav))[0]
			i1 = int(com[ui,yi]); i2 = int(i1 + 1); r = com[ui,yi] - i1
	
			# For Gaussian flows, need to calculate EEF about new center of mass.
			# This requires calculation of EEF at two grid points to aid continuity.
			EEF1, l_PV = PV.EEF(P_xav,y_nd,y_nd[i1],i1,dy_nd,N)
			EEF2, l_PV = PV.EEF(P_xav,y_nd,y_nd[i2],i2,dy_nd,N)
			EEF_array[ui,yi,:] = (1 - r) * EEF1 + r * EEF2			

			#EEF_array[ui,yi,:], l_PV = PV.EEF(P_xav,y_nd,y_nd[iii],iii,dy_nd,N)

	filename = 'EEF_array_' + str(pi);
	np.save(filename,EEF_array)
	
	filename_com = 'com' + str(pi)
	np.save(filename_com,com)


if __name__ == '__main__':
	jobs = []
	for pi in range(0,pe):
		p = mp.Process(target=EEF_main,args=(sets[pi],pi))
		jobs.append(p)
		p.start()

	for p in jobs:
		p.join()

# Now collect the results by reloading them, and compile into one array
EEF_array = np.zeros((Nu,nn,2))
com = np.zeros((Nu,nn))
yn_count = 0
for pi in range(0,pe):
	filename = 'EEF_array_' + str(pi) + '.npy'
	filename_com = 'com' + str(pi) + '.npy'

	EEF_array_tmp = np.load(filename)
	com_tmp = np.load(filename_com)	
	
	yn = np.shape(EEF_array_tmp)[0]
	EEF_array[yn_count:yn_count+yn,:,:] = EEF_array_tmp[:,:,:]
	com[yn_count:yn_count+yn,:] = com_tmp[:,:]
	yn_count = yn_count + yn

	os.remove(filename)
	os.remove(filename_com)	
		
np.save('EEF_array',EEF_array)
np.save('com',com)

elapsed = time.time() - start
elapsed = np.ones(1) * elapsed
print(elapsed)

	
	
	
