# EIG.py
#=======================================================

# This code finds the eigenmodes and eigenvalues for the 1-L shallow water system.

#====================================================

import numpy as np
import matplotlib.pyplot as plt

import diagnostics
import eigDiagnostics
import PV
import eigSolver
from output import ncSaveEigenmodes

from inputFile_1L import *

# 1L eigenmode solver
#====================================================
#====================================================

# Numbers to note:
# Dimensional forcing periods and corresponding 0-D omega:
# 50 days --> omega = 0.888...
# 60 days --> omega = 0.740740...
# 70 days --> omega = 0.634920...
 
I = np.complex(0.0,1.0);

# Define the coefficients required by the solver
a1,a2,a3,a4,b1,b4,c1,c2,c3,c4 = eigSolver.EIG_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,gamma_nd,dy_nd,N);
	
k_start = 2;
k_end = k_start + 1;
for ii in range(k_start,k_end):
	# Run the solver for the current k-value.
	k = K_nd[ii];
	if BC == 'NO-SLIP':
		val, u_vec, v_vec, eta_vec = eigSolver.NO_SLIP_EIG(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,N,N2,ii,True);
	if BC == 'FREE-SLIP':
		val, vec = eigSolver.FREE_SLIP_EIG(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,N,N2,ii,False);

	dim = np.size(val);
	# count = number of zero-crossings by the eigenvector 
	# i_count = set of indices ordering modes by their count
	count, i_count = eigDiagnostics.orderEigenmodes(vec,val,N,False);
	
	# Order all relevant arrays
	count = count[i_count];
	vec = vec[:,i_count];
	val = val[i_count];
	
	u_vec, v_vec, eta_vec = eigDiagnostics.vec2vecs(vec,N,dim,BC);

	#====================================================

	# By this point, eigenmodes have been calculated and ordered. The ordering algorithm may make mistakes,
	# however, so the set of modes can be manually updated with the remainder of the code.

	# 1. This first section runs through each mode, allowing to manually update the zero-crossings count.
	# To quit updating the vector's counts, type end.
	# We use the u vector as an example.

	u = np.zeros((N,N),dtype=float);
	count_new = np.array(list(count));
	update_i = [];		# Used to store the set of wi indices that need updating.
	wii = 0;
	while wii < dim:
		print('i_count = ' + str(wii));
		for i in range(0,N):
			for j in range(0,N):
				u[j,i] = np.real(u_vec[j,wii] * np.exp(2 * np.pi * I * (k * x_nd[i])));
		print('count = ' + str(count[wii]));
		plt.subplot(121);
		plt.contourf(u);
		plt.subplot(122);
		plt.plot(np.abs(u_vec[:,wii]),y_nd);
		plt.ylim(-0.5,0.5);
		plt.show();

		# Can use these plots to check any errors, and update count accordingly;
		# comment out if not needed.
		update_count = raw_input('-->');		# The first step updates count, but i_count will no longer match.									
		if update_count != '' and update_count != 'end':	
			print('updated');
			count_new[wii] = int(update_count);	# Stores the new count, to be used to rearrange the vectors.
 			update_i.append(wii);				# Stores the indices of the vectors to be updated.
		elif update_count == 'end':
			wii = wii + dim; 	# End the loop, and don't update any more modes.
		
		wii = wii + 1;

	update_i = np.array(update_i); 	# Convert it into a usable array

	#====================================================

	# 2. According to the updated count, the vectors are reordered.

	i_count_new = np.argsort(count_new);
	count = count_new[i_count_new];
	
	# Update the vectors & eigenvalues	
	vec = vec[:,i_count_new];
	val = val[i_count_new];

	u_vec, v_vec, eta_vec = eigDiagnostics.vec2vecs(vec,N,dim,BC);		

	ncSaveEigenmodes(vec,val,count,y_nd,k,N,dim,BC);
	
	#====================================================
	
	# 3. The new vectors can be checked, again using u_vec as an example.

	u = np.zeros((N,N),dtype=float);
	for wi in range(0,dim):
		for i in range(0,N):
			for j in range(0,N):
				u[j,i] = np.real(u_vec[j,wi] * np.exp(2 * np.pi * I * (k * x_nd[i])));
		print(count_new[wi]);
		u_abs = np.abs
		plt.subplot(121);
		plt.contourf(u);
		plt.subplot(122);
		plt.plot(np.abs(u_vec[:,wi]),y_nd);
		plt.ylim(-0.5,0.5);
		plt.show();

	#====================================================
	#====================================================

	


	# val is a set of dimensionless frequencies; it's easier to work in terms of dimensional periods.
	# Here we order dimensionalise val and define the index arrays that give orderings of val.
	freq = np.real(val);				# Frequencies
	period = 1. / freq;					# Dimensionless period
	growth = np.imag(val);				# Growth rates
	period_days = T_adv / (freq * 24. * 3600.);		# Periods corresponding to eigenvalues (days).
	i_freq = np.argsort(np.abs(freq));						# Indices of frequencies in ascending order.
	i_growth = np.argsort(growth);					# Indices of growth rates in ascending order.

	# We want to look at a eigenfunctions with eigenvalues within some frequency/period range, given by pmin, pmax.
	pmin = - 100;		# Min and max period in days
	pmax = - 40;
	# Some empty lists, to be added to and converted to numpy arrays.
	freq_set = [];
	freq_index = [];
	u_set = [];
	v_set = [];
	eta_set = [];
	for wi in range(0,dim):
		#plt.plot(u_vec[:,wi]);
		#plt.show();
		if pmin <= period_days[wi] <= pmax:
		#if pmin <= period_days[wi] <= pmax:
			freq_index.append(wi);
			freq_set.append(freq[wi]);
			u_set.append(u_vec[:,wi]);
			v_set.append(v_vec[:,wi]);
			eta_set.append(eta_vec[:,wi]);
	u_set = np.array(u_set);
	v_set = np.array(v_set);
	eta_set = np.array(eta_set);
	# Note that using np.array reverses the indexing convention.

	# Now define u, v and eta in (x,y)-space.
	Nf = len(freq_set);
	u = np.zeros((N,N,Nf));
	v = np.zeros((N,N,Nf));
	eta = np.zeros((N,N,Nf));
	for wi in range(0,Nf):
		for i in range(0,N):
			for j in range(0,N):
				u[j,i,wi] = np.real(u_set[wi,j] * np.exp(2 * np.pi * I * (k * x_nd[i])));
				v[j,i,wi] = np.real(v_set[wi,j] * np.exp(2 * np.pi * I * (k * x_nd[i])));
				eta[j,i,wi] = np.real(eta_set[wi,j] * np.exp(2 * np.pi * I * (k * x_nd[i])));

		u_full = np.zeros((N,N));
		eta_full = np.zeros((N,N));
		for i in range(0,N):
			u_full[:,i] = u[:,i,wi] + U0_nd[:];
			eta_full[:,i] = eta[:,i,wi] + H0_nd[:];

		PV_full, PV_prime = eigDiagnostics.PV(u[:,:,wi],v[:,:,wi],eta[:,:,wi],u_full,eta_full,H0_nd,U0_nd,f_nd,dx_nd,dy_nd,N);
		print(period_days[freq_index[wi]]);
		P, P_xav = eigDiagnostics.footprint(u[:,:,wi],v[:,:,wi],PV_full,x_nd,dx_nd,dy_nd,N);
		
		eigDiagnostics.eigPlots(u[:,:,wi],v[:,:,wi],eta[:,:,wi],x_nd,y_nd);
		#diagnostics.pvPlots(PV_full,PV_prime,P,x_nd,y_nd);

		plt.plot(y_nd,u_set[wi,:]);
		plt.show();

#====================================================

# Selects all the frequencies inside some desired range, e.g. unstable modes.
#====================================================

#scatter_set = [];
for jj in range(0,dim):
	if -0.015 <= np.imag(val[jj]):
		scatter_set.append([k,24*3600*val[jj]/T_adv]);

scatter_set = np.array(scatter_set);	
print(scatter_set);	
	
plt.scatter(scatter_set[:,0],np.real(scatter_set[:,1]));
plt.xlabel('Wavenumber',fontsize=18);
plt.ylabel('Frequency (day-1)',fontsize=18);	
plt.show();
		











	
