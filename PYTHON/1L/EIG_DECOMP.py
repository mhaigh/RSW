# EIG_DECOMP.py
#=======================================================

# This code decomposes a 1L SW solution, as produced by RSW_1L.py into into eigenmodes components,
# producing a set of weights corresponding the to set of eigenmodes found by runnning modules located in eigSolver.py
# First it defines the solution, which can either be defined from a file (FILE) or by running RSW_1L.py (NEW).

#====================================================

import sys

import numpy as np
#import matplotlib.pyplot as plt
import itertools as it

from eig import eigSolver, eigDiagnostics
from core import diagnostics, solver, forcing, energy, PV
from output import output, output_read

from inputFile import *

#====================================================

def EIG_DECOMP_main(U0_nd,H0_nd,dim):

	# Dimensions
	Nm = dim;		# How many modes to use in the decomposition at each wavenumber (dim is maximum).
	Nk = N;

	# The 1L SW solution
	#====================================================

	I = np.complex(0.0,1.0);

	SOL = 'NEW';

	# Define the solution in (k,y)-space - can be from FILE or a NEW run.
	if SOL == 'FILE':
		solution = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/REF/solution_NORTH256.npy');
	if SOL == 'NEW':
		# Forcing
		#if FORCE_TYPE == 'CTS':
		#	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing.forcing_cts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
		#elif FORCE_TYPE == 'DCTS':
		#	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing.forcing_dcts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd)
		#elif FORCE_TYPE == 'DELTA':
		#	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing.forcing_delta(AmpF_nd,y0_index,dx_nd,N);
		# Coefficients
		a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
		# Solver
		if BC == 'NO-SLIP':
			solution = solver.NO_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2);
		elif BC == 'FREE-SLIP':
			solution = solver.FREE_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2);
		else:
			sys.exit('ERROR: choose valid BC');

	solution = solution / AmpF_nd

	print('solved');

	#====================================================
	
	# Eigenmode analysis
	#====================================================
	#====================================================
	
	# Initisialastion steps
	#====================================================

	VEC = 'FILE';		# From FILE, requires pre-saved vectors which take up lots of memory.
	LOOP = 'FULL';		# FULL, PART

	
	if LOOP == 'FULL':
		loop = range(0,N);
		Nk_neg = 9;
		Nk_pos = 8;
	elif LOOP == 'PART':
		Nk_neg = 6; Nk_pos = 6;		# How many positive/negative wavenumbers to perform this decomposition at.
		loop = it.chain(range(0,Nk_pos+1),range(N-Nk_neg,N));
		Nk = Nk_neg + Nk_pos + 1;
	else:
		sys.exit('ERROR: LOOP must be FULL or FILE.');


	theta = np.zeros((Nm,Nk),dtype=complex); 	# Initialise the set of weights; these will be complex.
	proj = np.zeros((dim,N),dtype=complex);		# The projection. Sums the Nm most dominant modes, each of length dim, for Nk i-values.
	dom_index = np.zeros((Nm,Nk),dtype=int);	# To store the indices of the Nm-most dominant modes.

	var = np.zeros((Nk));
	mean = np.zeros((Nk));

	scatter_k = np.zeros(Nm * Nk);		# An empty array for saving k-values, for use in the scatter plot of dominant modes.
	scatter_l = np.zeros(Nm * Nk);		# An empty array for storing the count, psuedo-wavenumber l.
	scatter_p = np.zeros(Nm * Nk);		# An empty array for storing periods of the dominant wavenumbers.
	theta_abs = np.zeros((Nm,Nk));		# For storing the absolute value of each weight.
	theta_abs_tot = np.zeros(Nk);		# For storing sum of absolute values of each set of decomposition weights.
	cx = np.zeros((Nm,Nk));				# For storing the zonal phase speed of each mode,
	cy = np.zeros((Nm,Nk));				# and the meridional phase speed.
	p = np.zeros((Nk,2));					# For storing weighted phase speeds at each wavenumber. (i=0 => y)

	# Analysis
	#====================================================

	# Loop over desired wavenumbers (for tests, this may not be the full range of wavenumbers)
	# ii indexes arrays storing information at ALL wavenumbers k
	# i indexes arrays storing information ONLY at wavenumbers used in the decomposition.
	for ii in loop:	 	
	
		print(' ');
		print('ii = ' + str(ii));
		k = K_nd[ii];
		print('k = ' + str(k));
		i = k % Nk;
		print('i = ' + str(i));
		i = int(i+0.01);
	
		# Eigenmodes, eigenvalues and count.
		#====================================================

		# This section returns three arrays: 1. val, 2. vec, 3. count
		# 1.) val = val[0:dim] stores the eigenvalues/frequencies.
		# 2.) vec = vec[]
		# 3.) count = count[]

		# Run the solver for the current k-value.
		if VEC == 'NEW':	# Solve the eigenmode problem anew.
			a1,a2,a3,a4,b1,b4,c1,c2,c3,c4 = eigSolver.EIG_COEFFICIENTS2(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,gamma_nd,dy_nd,N);
			if BC == 'NO-SLIP':
				val, u_vec, v_vec, h_vec = eigSolver.NO_SLIP_EIG(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,N,N2,ii,True);
			if BC == 'FREE-SLIP':
				val, vec = eigSolver.FREE_SLIP_EIG(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,N,N2,ii,False);

			# Order modes by meridional pseudo wavenumber (count).
			count, i_count = eigDiagnostics.orderEigenmodes(vec,val,x_nd,k,T_nd[ts],N,dim,BC)
			#count, i_count = eigDiagnostics.orderEigenmodes2(vec,val,N,False)
			count = count[i_count]
			vec = vec[:,i_count]
			val = val[i_count]

			# Each eigenmode is currently a unit vector, but we normalise so that each mode contains unit energy.
			#==
	
			# Extract the three components.
			u_vec, v_vec, h_vec = eigDiagnostics.vec2vecs(vec,N,dim,BC);	
		
			# Calculate the contained in each component.
			E = np.zeros(dim);	
			for wi in range(0,dim):
				EE = energy.E_anomaly_EIG(u_vec[:,wi],v_vec[:,wi],h_vec[:,wi],H0_nd,U0_nd,Ro,y_nd,dy_nd);
				# Normalise each vector by the square root of the energy.
				u_vec[:,wi], v_vec[:,wi], h_vec[:,wi] = u_vec[:,wi] / np.sqrt(EE), v_vec[:,wi] / np.sqrt(EE), h_vec[:,wi] / np.sqrt(EE);

			# Rebuild the vector. This should have unit energy perturbation. 
			# (There are more direct ways of executing this normalisation, but this method is the safest.)
			vec = eigDiagnostics.vecs2vec(u_vec,v_vec,h_vec,N,dim,BC);
		
			# Comment out this line, depending on which EIG_COEFFICIENTS function is being called.
			#val = val / (2. * np.pi * I * Ro);

		elif VEC == 'FILE':	# Load eigenmodes and eigenvalues from file.
	
			# Low-res
			path = '/home/mike/Documents/GulfStream/RSW/DATA/1L/EIG/128/west/';
			#path = '/home/mike/Documents/GulfStream/RSW/DATA/1L/EIG/128/nu='+str(int(nu))+'/';
			ncFile = path + 'RSW1L_Eigenmodes_k' + str(int(k)) + '_N129.nc';
		
			# High-res
			#path = '/media/mike/Seagate Expansion Drive/Documents/GulfStream/RSW/DATA/1L/EIG/256/16/'
			#path =  '/home/mike/Documents/GulfStream/RSW/DATA/1L/EIG/256/west/';
			#ncFile = path + 'RSW1L_Eigenmodes_k' + str(int(k)) + '_N257.nc';
			print('Reading from ' + ncFile + '...');
			val, vec, count = output_read.ncReadEigenmodes(ncFile);
		else:
			sys.exit('VEC must be FILE or NEW');
	
		# Expresses the eigenvalues (frequencies) in terms of periods (units days).
		freq = np.real(val);
		period_days = T_adv / (freq * 24. * 3600.);

		dim = np.size(val);


		#====================================================
	
		# Now we have the solution and the eigenmodes.
		# The decomposition follows the following steps:
		# 1. Define the solution to be decomposed as Phi.
		# 2. Decompose Phi into vec using a linear solver; theta_tmp stores the weights.
		# 3. Arrange the weights in descending order according to their complex amplitude.
		# 4. Sum the Nm-most dominant weights.

		Phi = solution[:,ii];		# 1. Assign the solution corresponding to wavenumber k=K_nd[ii].

		theta_tmp = np.linalg.solve(vec,Phi); 				# 2.
		theta_abs_tmp = np.abs(theta_tmp);
		dom_index_tmp = np.argsort(-theta_abs_tmp);			# 3. The indices of the modes, ordered by 'dominance'.
		theta_abs_tot[i] = sum(theta_abs_tmp[dom_index_tmp[0:dim]]);	# Change dim to Nm if necessary

		# Now loop over each mode (at wavenumber k)
		for mi in range(0,Nm):
			#print(dom_index_tmp[mi]);
			#print('count = ' + str(count[dom_index_tmp[mi]]));
			#print(np.abs(theta_tmp[dom_index_tmp[mi]]));

			dom_index[mi,i] = dom_index_tmp[mi];
			theta[mi,i] = theta_tmp[dom_index_tmp[mi]];
			# All weights are now ordered in terms of their absolute value.

			# Absolute value of each mode
			theta_abs[mi,i] = np.abs(theta[mi,i]);		

			# Zonal & meridional phase speed of each mode
			cx[mi,i] = freq[dom_index[mi,i]] / k;	
			if count[dom_index[mi,i]] != 0:
				cy[mi,i] = freq[dom_index[mi,i]] / count[dom_index[mi,i]];
			
			if mi < 95:
				# The projection.
				proj[:,ii] = proj[:,ii] + theta_tmp[dom_index_tmp[mi]] * vec[:,dom_index_tmp[mi]];	# 4.
		
			# Scatter plot arrays.
			scatter_k[i*Nm+mi] = k;	
			scatter_l[i*Nm+mi] = count[dom_index[mi,i]];
			scatter_p[i*Nm+mi] = period_days[dom_index[mi,i]];
			#plt.plot(vec[0:N,dom_index_tmp[mi]],y_nd);
			#plt.ylim(-0.5,0.5);
			#plt.show();
		#plt.plot(theta_abs[:,i]);
		#plt.show();
	
		# Statistics: mean, variance, weighted average 
		# Should normalise so that all have the same mean (i.e. mean = 1/Nm);
		mean[i] = theta_abs_tot[i] / Nm;
		var[i] = sum((theta_abs[:,i] - mean[i])**2) / (Nm * mean[i]**2);
		p[i,1] = sum((theta_abs[:,i] * cx[:,i])) / theta_abs_tot[i];
		p[i,0] = sum((theta_abs[:,i] * cy[:,i])) / theta_abs_tot[i];
		print(cy[10,i]);

	return theta, mean, var, p, proj, solution, scatter_k, scatter_l, scatter_p

# End function
#====================================================

# Preamble before looping through function.
#====================================================

# Define BG flow set. 
nn = 1;
U0_set = np.linspace(0.04,0.04,nn);

# Define forcing. This is constant throughout if only BG flow varies.
if FORCE_TYPE == 'CTS':
	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing.forcing_cts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd)
elif FORCE_TYPE == 'CTS2':
	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing.forcing_cts2(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,bh,dx_nd,dy_nd);
elif FORCE_TYPE == 'DCTS':
	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing.forcing_dcts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
else:
	sys.exit('ERROR: Invalid forcing option selected.');

Nm = dim;
Nk = N;

# Initialise output arrays
#====================================================

solution = np.zeros((Nm,Nk,nn),dtype=complex);
theta = np.zeros((Nm,Nk,nn),dtype=complex); 	# Initialise the set of weights; these will be complex.
proj = np.zeros((dim,N,nn),dtype=complex);		# The projection. Sums the Nm most dominant modes, each of length dim, for Nk i-values.
dom_index = np.zeros((Nm,Nk,nn),dtype=int);	# To store the indices of the Nm-most dominant modes.

var = np.zeros((Nk,nn));
mean = np.zeros((Nk,nn));

scatter_k = np.zeros((Nm*Nk,nn));		# An empty array for saving k-values, for use in the scatter plot of dominant modes.
scatter_l = np.zeros((Nm*Nk,nn));		# An empty array for storing the count, psuedo-wavenumber l.
scatter_p = np.zeros((Nm*Nk,nn));		# An empty array for storing periods of the dominant wavenumbers.
theta_abs = np.zeros((Nm,Nk,nn));		# For storing the absolute value of each weight.
theta_abs_tot = np.zeros((Nk,nn));		# For storing sum of absolute values of each set of decomposition weights.
c = np.zeros((Nm,Nk,nn));				# For storing the zonal phase speed of each mode.
p = np.zeros((Nk,2,nn));					# For storing weighted phase speed at each wavenumber.

# Main function
#====================================================

if __name__ == '__main__':

	for ui in range(0,nn):

		# Before executing the main function, (re)define BG state.
		for j in range(0,N):
			U0[j] = U0_set[ui];
			H0[j] = - (U0[j] / g) * (f0 * y[j] + beta * y[j]**2 / 2) + Hflat;
			U0_nd = U0 / U;
			H0_nd = H0 / chi;

		# Decompose solution into eigenmodes.
		theta[:,:,ui], mean[:,ui], var[:,ui], p[:,:,ui], proj[:,:,ui], solution[:,:,ui], scatter_k[:,ui], scatter_l[:,ui], scatter_p[:,ui] = EIG_DECOMP_main(U0_nd,H0_nd,dim);


theta_abs = np.absolute(theta[:,:,0])
plt.plot(theta_abs[:,10])
plt.plot(theta_abs[:,40])
plt.plot(theta_abs[:,80])
plt.plot(theta_abs[:,120])
plt.show()

#sys.exit()
# Save data

np.save('theta',theta);
np.save('mean',mean);
np.save('var',var);
np.save('p',p);

#sys.exit();

# Plotting
#====================================================

if False:
	#plt.plot(np.fft.fftshift(K_nd),np.fft.fftshift(p[:,1,:]*U),label='zonal');
	plt.plot(np.fft.fftshift(K_nd),np.fft.fftshift(p[:,0,:]*U),label='merid');
	plt.xlabel('k');
	plt.ylabel('Weighted phase speed');
	plt.legend();
	plt.grid();
	plt.show();


	plt.plot(np.fft.fftshift(K_nd),np.fft.fftshift(mean[:,0]));
	plt.ylabel('MEAN ABS WEIGHT');
	plt.xlabel('WAVENUMBER');
	plt.show();

	plt.plot(np.fft.fftshift(K_nd),np.fft.fftshift(var[:,0]));
	plt.ylabel('VARIANCE');
	plt.xlabel('WAVENUMBER');
	plt.show();

#sys.exit();

#====================================================

# Intialise the projections
utilde_proj = np.zeros((N,N),dtype=complex);
vtilde_proj = np.zeros((N,N),dtype=complex);
etatilde_proj = np.zeros((N,N),dtype=complex);

proj = np.squeeze(proj,axis=2)
# Now look at the solution given by the dominant modes
if BC == 'FREE-SLIP':
	for j in range(0,N):
		utilde_proj[j,:] = proj[j,:];
		etatilde_proj[j,:] = proj[N+N2+j,:];
	for j in range(0,N2):
		vtilde_proj[j+1,:] = proj[N+j,:];

if BC == 'NO-SLIP':
	for j in range(0,N2):
		utilde_proj[j+1,:] = proj[j,:];
		vtilde_proj[j+1,:] = proj[N2+j,:];
	for j in range(0,N):		
		etatilde_proj[j,:] = proj[2*N2+j,:];

u_proj, v_proj, eta_proj = solver.SPEC_TO_PHYS(utilde_proj,vtilde_proj,etatilde_proj,T_nd,dx_nd,omega_nd,N);

u_proj = np.real(u_proj);
v_proj = np.real(v_proj);
eta_proj = np.real(eta_proj);

u_full = np.zeros((N,N,Nt));
h_full = np.zeros((N,N,Nt));
for i in range(0,N):
	for ti in range(0,Nt):
		u_full[:,i,ti] = u_proj[:,i,ti] + U0_nd[:];
		h_full[:,i,ti] = eta_proj[:,i,ti] + H0_nd[:];


utilde_nd, vtilde_nd, etatilde_nd = solver.extractSols(np.squeeze(solution,axis=2),N,N2,BC);
u, v, h = solver.SPEC_TO_PHYS(utilde_nd,vtilde_nd,etatilde_nd,T_nd,dx_nd,omega_nd,N);

u = np.real(u);
v = np.real(v);
h = np.real(h);

error_norm = np.sum(u[:,:,ts]**2)
error = np.sum((u[:,:,ts] - u_proj[:,:,ts])**2) / error_norm
print(error)

#====================================================

eigDiagnostics.eigPlots(u_proj[:,:,ts],v_proj[:,:,ts],eta_proj[:,:,ts],u[:,:,ts],v[:,:,ts],h[:,:,ts],x_nd,y_nd,x_grid,y_grid,True);

#====================================================

# PV

PV_prime, PV_full, PV_BG = PV.potentialVorticity(u_proj,v_proj,eta_proj,u_full,h_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd,Ro);
#PV_prime1, PV_prime2, PV_prime3 = PV.potentialVorticity_linear(u,v,h,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd,Ro);
uq, Uq, uQ, UQ, vq, vQ = PV.fluxes(u,v,U0_nd,PV_prime,PV_BG,N,Nt);
P, P_xav = PV.footprint(uq,Uq,uQ,UQ,vq,vQ,x_nd,T_nd,dx_nd,dy_nd,N,Nt);			

plt.plot(P_xav)
plt.show()

# We know that the selection of modes is dominated by factors other than the forcing freqeuncy.
sys.exit();

eigDiagnostics.scatterWeight(scatter_k,scatter_l,theta,theta_abs_tot,dom_index,Nm,Nk_neg,Nk_pos,Fpos);	
eigDiagnostics.scatterPeriod(scatter_k,scatter_l,scatter_p,dom_index,Nm,Nk_neg,Nk_pos,Fpos);	

#PV_full, PV_prime = eigDiagnostics.PV(u_proj,v_proj,eta_proj,u_full,h_full,H0_nd,U0_nd,f_nd,dx_nd,dy_nd,N);
#P, P_xav = eigDiagnostics.footprint(u_proj,v_proj,PV_full,x_nd,dx_nd,dy_nd,N);

#diagnostics.pvPlots(PV_full,PV_prime,P,x_nd,y_nd);

#plt.plot(y_nd,P_xav);
#plt.show();

#====================================================

	
	
	
