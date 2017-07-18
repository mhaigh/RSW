# EIG_DECOMP.py
#=======================================================

# This code decomposes a 1L SW solution, as produced by RSW_1L.py into into eigenmodes components,
# producing a set of weights corresponding the to set of eigenmodes found by runnning modules located in eigSolver.py
# First it defines the solution, which can either be defined from a file (FILE) or by running RSW_1L.py (NEW).

#====================================================

import numpy as np
import eigSolver
import matplotlib.pyplot as plt
import eigDiagnostics
import diagnostics
import solver

from inputFile_1L import *

#====================================================

# The 1L SW solution
#====================================================

SOL = 'NEW';

# Define the solution in (k,y)-space - can be from FILE or a NEW run.
if SOL == 'FILE':
	solution = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/REF/solution_NORTH256.npy');
if SOL == 'NEW':
	# Import the relevant modules.
	import forcing_1L
		# Forcing
	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.Forcing(x,y,K,y0,r0,N,FORCE,AmpF,g,f,f0,U,L,dx,dy);
	# Coefficients
	a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
	# Solver
	if BC == 'NO-SLIP':
		solution = solver.NO_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N,N2);
	elif BC == 'FREE-SLIP':
		uBC, etaBC = solver.BC_COEFFICIENTS(Ro,Re,f_nd,H0_nd,dy_nd,N);
		solution = solver.FREE_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,uBC,etaBC,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N,N2);
	else:
		print('ERROR: choose BCs');

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

print('solved');

#====================================================

# Now we have the solution. The remainder of the code has multiple functions, which can be selected by defining OPT.
# OPT = 1 --> 
# OPT = 2 -->
#====================================================

VEC = 'NEW';	# From FILE, requires pre-saved vectors which take up lots of memory.


Nm = 3;		# How many modes to use in the decomposition at each wavenumber
Nk = Nm;	# How many positive/negative wavenumbers to perform this decomposition at,
			# totataling 

theta = np.zeros((Nm,N),dtype=complex); 		# Initialise the set of weights; these will be complex.
proj = np.zeros((dim,N),dtype=complex);			# The projection. Sums the Nm most dominant modes, each of length dim, for N i-values.
val = np.zeros((Nm,N),dtype=complex);			# Eigenvalues, will only need to save Nm of them.
vec = np.zeros((Nm,dim,N),dtype=complex);		# The eigenvectors

# Define the coefficients that the solver requires
a1,a2,a3,a4,b1,b4,c1,c2,c3,c4 = eigSolver.EIG_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,gamma_nd,dy_nd,N);
uBC, etaBC = eigSolver.BC_COEFFICIENTS(Ro,Re,f_nd,H0_nd,dy_nd,N);

# Loop over all wavenumbers
for i in range(0,Nk+1):
	print(i);

	Phi = solution[:,i];		# Assign the solution corresponding to wavenumber k=K_nd[i].
	
	theta_tmp, val_tmp, vec_tmp = eigSolver.eigDecomp(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,N,N2,i,BC,VEC,Phi);

	dom_index = np.argsort(-(np.abs(theta_tmp))**2);	# The indices of the modes, ordered by 'dominance'.
	for mi in range(0,Nm):
		val[mi,i] = val_tmp[dom_index[mi]];
		theta[mi,i] = theta_tmp[dom_index[mi]];
		vec[mi,:,i] = vec_tmp[:,dom_index[mi]];
		proj[:,i] = proj[:,i] + theta[mi,i] * vec[mi,:,i];

for i in range(N-Nk-1,N):
	print(i);

	Phi = solution[:,i];		# Assign the solution corresponding to wavenumber k=K_nd[i].
	
	theta_tmp, val_tmp, vec_tmp = eigSolver.eigDecomp(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,N,N2,i,BC,VEC,Phi);

	dom_index = np.argsort(-(np.abs(theta_tmp))**2);	# The indices of the modes, ordered by 'dominance'.
	for mi in range(0,Nm):
		val[mi,i] = val_tmp[dom_index[mi]];
		theta[mi,i] = theta_tmp[dom_index[mi]];
		vec[mi,:,i] = vec_tmp[:,dom_index[mi]];
		proj[:,i] = proj[:,i] + theta[mi,i] * vec[mi,:,i];

	print(T_adv/(np.real(val[0,i]*24*3600)));
	#plt.plot(vec[0,0:N,i]);
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

u_proj = u_proj[:,:,ts];
v_proj = v_proj[:,:,ts];
eta_proj = eta_proj[:,:,ts];

u_full = np.zeros((N,N));
eta_full = np.zeros((N,N));
for i in range(0,N):
	u_full[:,i] = u_proj[:,i] + U0_nd[:];
	eta_full[:,i] = eta_proj[:,i] + H0_nd[:];

#====================================================

eigDiagnostics.eigPlots(u_proj,v_proj,eta_proj,u_nd[:,:,ts],v_nd[:,:,ts],eta_nd[:,:,ts],x_nd,y_nd,True);

#====================================================

#PV_full, PV_prime = eigDiagnostics.PV(u_proj,v_proj,eta_proj,u_full,eta_full,H0_nd,U0_nd,f_nd,dx_nd,dy_nd,N);
#P, P_xav = eigDiagnostics.footprint(u_proj,v_proj,PV_full,x_nd,dx_nd,dy_nd,N);

#diagnostics.pvPlots(PV_full,PV_prime,P,x_nd,y_nd);

#plt.plot(y_nd,P_xav);
#plt.show();

#====================================================

	
	
	
