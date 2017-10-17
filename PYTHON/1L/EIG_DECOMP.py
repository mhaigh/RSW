# EIG_DECOMP.py
#=======================================================

# This code decomposes a 1L SW solution, as produced by RSW_1L.py into into eigenmodes components,
# producing a set of weights corresponding the to set of eigenmodes found by runnning modules located in eigSolver.py
# First it defines the solution, which can either be defined from a file (FILE) or by running RSW_1L.py (NEW).

#====================================================

import sys

import numpy as np
import matplotlib.pyplot as plt
import itertools as it

import eigSolver
import eigDiagnostics
import diagnostics
import solver
import output
import output_read

from inputFile_ref import *

#====================================================

# The 1L SW solution
#====================================================

I = np.complex(0.0,1.0);

SOL = 'NEW';

# Define the solution in (k,y)-space - can be from FILE or a NEW run.
if SOL == 'FILE':
	solution = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/REF/solution_NORTH256.npy');
if SOL == 'NEW':
	# Import the relevant modules.
	import forcing_1L
		# Forcing
	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_cts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
	# Coefficients
	a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
	# Solver
	if BC == 'NO-SLIP':
		solution = solver.NO_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N,N2);
	elif BC == 'FREE-SLIP':
		uBC, etaBC = solver.BC_COEFFICIENTS(Ro,Re,f_nd,H0_nd,dy_nd,N);
		solution = solver.FREE_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2);
	else:
		sys.exit('ERROR: choose valid BC');

# Split the solution into its components
utilde_nd = np.zeros((N,N),dtype=complex);
vtilde_nd = np.zeros((N,N),dtype=complex);
etatilde_nd = np.zeros((N,N),dtype=complex);
if BC == 'FREE-SLIP':
	for j in range(0,N):
		utilde_nd[j,:] = solution[j,:];
		etatilde_nd[j,:] = solution[N+N2+j,:];
	for j in range(0,N2):
		vtilde_nd[j+1,:] = solution[N+j,:];
if BC == 'NO-SLIP':
	for j in range(0,N2):
		utilde_nd[j+1,:] = solution[j,:];
		vtilde_nd[j+1,:] = solution[N2+j,:];
	for j in range(0,N):		
		etatilde_nd[j,:] = solution[2*N2+j,:];
u_nd, v_nd, eta_nd = solver.SPEC_TO_PHYS(utilde_nd,vtilde_nd,etatilde_nd,T_nd,dx_nd,omega_nd,N);

u_nd = np.real(u_nd);
v_nd = np.real(v_nd);
eta_nd = np.real(eta_nd);

print('solved');

#====================================================

# The eigmodes, eigenvalues, and count (all ordered)
#====================================================

VEC = 'FILE';		# From FILE, requires pre-saved vectors which take up lots of memory.
LOOP = 'FULL';		# FULL, PART

Nm = 8;						# How many modes to use in the decomposition at each wavenumber (dim is maximum).
if LOOP == 'FULL':
	loop = range(0,N);
	Nk = N;
	# 
	Nk_neg = 6;
	Nk_pos = 6;
elif LOOP == 'PART':
	Nk_neg = 6; Nk_pos = 6;		# How many positive/negative wavenumbers to perform this decomposition at.
	loop = it.chain(range(0,Nk_pos+1),range(N-Nk_neg,N));
	Nk = Nk_neg + Nk_pos + 1;
else:
	sys.exit('ERROR: LOOP must be FULL or FILE.');

theta = np.zeros((Nm,Nk),dtype=complex); 	# Initialise the set of weights; these will be complex.
proj = np.zeros((dim,N),dtype=complex);		# The projection. Sums the Nm most dominant modes, each of length dim, for Nk i-values.
dom_index = np.zeros((Nm,Nk),dtype=int);	# To store the indices of the Nm-most dominant modes.

scatter_k = np.zeros(Nm * Nk);		# An empty array for saving k-values, for use in the scatter plot of dominant modes.
scatter_l = np.zeros(Nm * Nk);		# An empty array for storing the count, psuedo-wavenumber l.
scatter_p = np.zeros(Nm * Nk);	# An empty array for storing periods of the dominant wavenumbers.
theta_abs_tot = np.zeros(Nk);		# For storing sum of absolute values of each set of decomposition weights.

# Loop over desired wavenumbers (for tests, this may not be the full range of wavenumbers)
# ii indexes arrays storing information at ALL wavenumbers k
# i indexes arrays storing information ONLY at wavenumbers used in the decomposition.
for ii in loop:	 	

	print(' ');
	print('ii = ' + str(ii));
	k = K_nd[ii];
	print('k = ' + str(k));
	i = k % Nk
	print('i = ' + str(i));
	i = int(i+0.01);

	# Run the solver for the current k-value.
	if VEC == 'NEW':	# Solve the eigenmode problem anew.
		a1,a2,a3,a4,b4,c1,c2,c3,c4 = eigSolver.EIG_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,gamma_nd,dy_nd,N);
		if BC == 'NO-SLIP':
			val, vec = eigSolver.NO_SLIP_EIG(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,N,N2,i,False);
			count = np.zeros(dim);
		if BC == 'FREE-SLIP':
			val, vec = eigSolver.FREE_SLIP_EIG(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,N,N2,ii,False);
			count = np.zeros(dim);
	elif VEC == 'FILE':	# Load eigenmodes and eigenvalues from file.
		path = '/home/mike/Documents/GulfStream/RSW/DATA/1L/EIG/128/nu='+str(int(nu))+'/';
		ncFile = path + 'RSW1L_Eigenmodes_k' + str(int(k)) + '_N129.nc';
		print('Reading from ' + ncFile + '...');
		val, vec, count = output_read.ncReadEigenmodes(ncFile);
	else:
		sys.exit('VEC must be FILE or NEW');
	
	val = val / (-2 * np.pi * I * Ro);
	# Expresses the eigenvalues (frequencies) in terms of periods.
	freq = np.real(val);
	period_days = T_adv / (freq * 24. * 3600.);
	print(np.transpose(period_days[dom_index[0:Nm,0]]));

	# This section returns three arrays: 1. val, 2. vec, 3. count
	# 1.) val = val[0:dim] stores the eigenvalues/frequencies.
	# 2.) vec = vec[]
	# 3.) count = count[]
	
	# Now we have the solution and the eigenmodes.
	# The decomposition follows the following steps:
	# 1. Define the solution to be decomposed as Phi.
	# 2. Decompose Phi into vec using a linear solver; theta_tmp stores the weights.
	# 3. Arrange the weights in descending order according to their complex amplitude.
	# 4. Sum the Nm-most dominant weights.

	Phi = solution[:,ii];		# 1. Assign the solution corresponding to wavenumber k=K_nd[i].

	theta_tmp = np.linalg.solve(vec,Phi); 				# 2.
	theta_abs_tmp = np.abs(theta_tmp);
	dom_index_tmp = np.argsort(-theta_abs_tmp);			# 3. The indices of the modes, ordered by 'dominance'.
	theta_abs_tot[i] = sum(theta_abs_tmp[dom_index_tmp[0:dim]]);	# Change dim to Nm if necessary
	for mi in range(0,Nm):
		#print(dom_index_tmp[mi]);
		#print('count = ' + str(count[dom_index_tmp[mi]]));
		#print(np.abs(theta_tmp[dom_index_tmp[mi]]));
		dom_index[mi,i] = dom_index_tmp[mi];
		theta[mi,i] = theta_tmp[dom_index_tmp[mi]];
		proj[:,ii] = proj[:,ii] + theta_tmp[dom_index_tmp[mi]] * vec[:,dom_index_tmp[mi]];	# 4.
		scatter_k[i*Nm+mi] = k;	
		scatter_l[i*Nm+mi] = count[dom_index[mi,i]];
		scatter_p[i*Nm+mi] = period_days[dom_index[mi,i]];
		#plt.plot(vec[0:N,dom_index_tmp[mi]],y_nd);
		#plt.ylim(-0.5,0.5);
		#plt.show();

	#plt.subplot(121);
	#plt.plot(np.real(proj[0:N,i]),y_nd);
	#plt.plot(np.real(Phi[0:N]),y_nd);
	#plt.ylim(-0.5,0.5);
	#plt.subplot(122);
	#plt.plot(vec[0:N,dom_index_tmp[0:Nm]],y_nd);
	#plt.ylim(-0.5,0.5);
	#plt.show();

#====================================================

# Intialise the projections
utilde_proj = np.zeros((N,N),dtype=complex);
vtilde_proj = np.zeros((N,N),dtype=complex);
etatilde_proj = np.zeros((N,N),dtype=complex);

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

u_proj = np.real(u_proj[:,:,ts]);
v_proj = np.real(v_proj[:,:,ts]);
eta_proj = np.real(eta_proj[:,:,ts]);

u_full = np.zeros((N,N));
eta_full = np.zeros((N,N));
for i in range(0,N):
	u_full[:,i] = u_proj[:,i] + U0_nd[:];
	eta_full[:,i] = eta_proj[:,i] + H0_nd[:];

#====================================================

eigDiagnostics.eigPlots(u_proj,v_proj,eta_proj,u_nd[:,:,ts],v_nd[:,:,ts],eta_nd[:,:,ts],x_nd,y_nd,x_grid,y_grid,True);

#====================================================

eigDiagnostics.scatterWeight(scatter_k,scatter_l,theta,theta_abs_tot,dom_index,Nm,Nk_neg,Nk_pos,Fpos);	
eigDiagnostics.scatterPeriod(scatter_k,scatter_l,scatter_p,dom_index,Nm,Nk_neg,Nk_pos,Fpos);	

#PV_full, PV_prime = eigDiagnostics.PV(u_proj,v_proj,eta_proj,u_full,eta_full,H0_nd,U0_nd,f_nd,dx_nd,dy_nd,N);
#P, P_xav = eigDiagnostics.footprint(u_proj,v_proj,PV_full,x_nd,dx_nd,dy_nd,N);

#diagnostics.pvPlots(PV_full,PV_prime,P,x_nd,y_nd);

#plt.plot(y_nd,P_xav);
#plt.show();

#====================================================

	
	
	
