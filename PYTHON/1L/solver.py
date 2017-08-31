# solver.py
#=======================================================

# Contains a set of functions to be called by RSW_1L_visc.py. 
# They constitute the solvers (for two choices of BCs) as well as a function that defines all parameters of the linear solver.

#=======================================================

import numpy as np
from diagnostics import diff

#=======================================================

# NO_SLIP_SOLVER_COEFFICIENTS
#=======================================================
def SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N):
# Here we also divide coeficients of derivatives by the relevant space-step: 
# 2*dy_nd for first-order derivatives, dy_nd**2 for second-order derivatives.

	I = np.complex(0,1);		# Define I = sqrt(-1)

	if Re == None:
		Ro_Re = 0;
	else:
		Ro_Re = Ro / Re;

	# Coefficients with no k- or y-dependence
	a2 = - Ro_Re / (dy_nd**2);	# Note: b3=a2, so instead of defining b3, we just use a2 in its place, saving time.
	b4 = 1. / (2. * dy_nd);

	# Coefficients with k-dependence only
	a4 = np.zeros(N,dtype=complex);		
	for i in range(0,N):
		a4[i] = 2. * np.pi * I * K_nd[i];

	# Coefficients with y-dependence only

	a3 = Ro * diff(U0_nd,2,0,dy_nd) - f_nd;		# For uniform U0_nd, the first term is zero	

	c2 = diff(H0_nd,2,0,dy_nd);					# For zero BG flow, H0_nd=Hflat=const, i.e. c2=0
	
	c3 = H0_nd / (2. * dy_nd);

	# a3, c2 and c3 have length N, so when used in defining u and v equations in the matrix below,
	# we must add 1 to each index in order to skip out the dead gridpoint.

	# Note that b1=f_nd so will not be defined, but we use f_nd directly.

	# Coefficients dependent on k and y
	delta = np.zeros((N,N));				# Appears in many of the coeffs, essentially a space-saver.
	a1 = np.zeros((N,N),dtype=complex);		# Note: b2=a1
	c1 = np.zeros((N,N),dtype=complex);
	for i in range(0,N):
		for j in range(0,N):
			delta[j,i] = 2. * np.pi * (omega_nd + U0_nd[j] * K_nd[i]);
			a1[j,i] = I * delta[j,i] * Ro + 4. * np.pi**2 * K_nd[i]**2 * Ro_Re + gamma_nd;
			c1[j,i] = 2. * np.pi * I * K_nd[i] * H0_nd[j];

	c4 = I * delta;

	return a1,a2,a3,a4,b4,c1,c2,c3,c4;

#=======================================================

# BC_COEFFICIENTS
#=======================================================
def BC_COEFFICIENTS(Ro,Re,f_nd,H0_nd,dy_nd,N):
# Some extra coefficients that are required by the free-slip solver in order to impose
# the du/dy=0 boundary conditions. These terms are not needed in the no-slip solver.

	if Re == None:
		Ro_Re = 0;
	else:
		Ro_Re = Ro / Re;

	uBC = 2 * Ro_Re / dy_nd**2;		# The extra term required for the u equation BCs.
	
	etaBC = np.zeros(2);
	if Re != None:
		etaBC[0] = f_nd[0] * H0_nd[0] * dy_nd * Re / (2 * Ro);		# The extra term required for the eta equation BCs.
		etaBC[1] = f_nd[N-1] * H0_nd[N-1] * dy_nd * Re / (2 * Ro);	

	return uBC, etaBC;

#=======================================================

# NO_SLIP_SOLVER 
#=======================================================
def NO_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N,N2):
# Called by RSW_1L.py if BC = 'NO-SLIP'.

	dim = 2 * N2 + N;
	#print(dim);

	A = np.zeros((dim,dim),dtype=complex);	# For the no-slip, no-normal flow BC.
	# eta has N gridpoints in y, whereas u and v have N2=N-2, after removing the two 'dead' gridpoints.
	# We primarily consider forcing away from the boundaries so that it's okay applying no-slip BCs.

	# Initialise the forcing.
	F = np.zeros((dim),dtype=complex);

	solution = np.zeros((dim,N),dtype=complex);

	for i in range(0,N):
		#print(i);

		# First the boundary terms. Some of these could be defined in the upcoming loop,
		# but for simplicity we define all boundary terms here.
	
		# u equation BCs
		# South
		A[0,0] = a1[1,i] - 2. * a2;			# u[1]
		A[0,1] = a2;						# u[2]
		A[0,N2] = a3[1];					# v[1]
		A[0,2*N2+1] = a4[i];				# eta[1] (eta[0] lies at i=2*N2 in A)
		# North
		A[N2-1,N2-1] = a1[N2,i] - 2. * a2;	# u[N2]=u[N-2]
		A[N2-1,N2-2] = a2;					# u[N2-1]=u[N-3]
		A[N2-1,2*N2-1] = a3[N2];			# v[N2]=v[N-2]
		A[N2-1,3*N2] = a4[i];				# eta[N2] (eta[N2] lies at i=3*N2 in A)

		# v equation BCs
		# South
		A[N2,0] = f_nd[1];					# u[1]
		A[N2,N2] = a1[1,i] - 2. * a2		# v[1]
		A[N2,N2+1] = a2;					# v[2]
		A[N2,2*N2] = - b4;					# eta[0] 
		A[N2,2*N2+2] = b4;					# eta[2]
		# North
		A[2*N2-1,N2-1] = f_nd[N2];					# u[N-2] 
		A[2*N2-1,2*N2-1] = a1[N2,i] - 2. * a2;		# v[N-2]
		A[2*N2-1,2*N2-2] = a2;						# v[N-3]
		A[2*N2-1,2*N2+N-1] = b4;					# eta[N-1]
		A[2*N2-1,2*N2+N-3] = - b4;					# eta[N-3]

		# eta equation BCs (Here we have to define BCs at j=2*N2,2*N2+1,3*N2+1,3*N2)
		# South
		A[2*N2,N2] = 2. * c3[0]					# v[1] (factor of 2 because we use one-sided FD for v_y at the boundaries)
		A[2*N2,2*N2] = c4[0,i];					# eta[0]
		A[2*N2+1,0] = c1[1,i];					# u[1]
		A[2*N2+1,N2] = c2[1];					# v[1]
		A[2*N2+1,N2+1] = c3[1];					# v[2] (factor of 1 here, back to using centered FD)
		A[2*N2+1,2*N2+1] = c4[1,i];				# eta[1]
		# North
		A[2*N2+N-1,2*N2-1] = - 2. * c3[N2+1];	# v[N-2] (factor of 2 because we use one-sided FD for v_y at the boundaries)
		A[2*N2+N-1,3*N2+1] = c4[N-1,i];			# eta[N-1]
		A[3*N2,N2-1] = c1[N2,i];  				# u[N-2]
		A[3*N2,2*N2-1] = c2[N2];				# v[N-2]
		A[3*N2,2*N2-2] = - c3[N2];				# v[N-3] (factor of 1 here, back to using centered FD)
		A[3*N2,3*N2] = c4[N2,i];				# eta[N-2]

		# Now the inner values - two loops required: one for u,v, and one for eta
		for j in range(1,N2-1):
			# u equation
			A[j,j] = a1[j+1,i] - 2. * a2;			# u[j+1]
			A[j,j-1] = a2;							# u[j]
			A[j,j+1] = a2;							# u[j+2]
			A[j,N2+j] = a3[j+1];					# v[j]
			A[j,2*N2+j+1] = a4[i];					# eta[j] (must add 1 to each eta index to skip out eta[0])
			
			# v equation
			A[N2+j,j] = f_nd[j+1];					# u[j+1]
			A[N2+j,N2+j] = a1[j+1,i] - 2. * a2;		# v[j+1]
			A[N2+j,N2+j-1] = a2;					# v[j]
			A[N2+j,N2+j+1] = a2;					# v[j+2]
			A[N2+j,2*N2+j] = - b4;					# eta[j]
			A[N2+j,2*N2+j+2] = b4;					# eta[j+2]

		for j in range(2,N2):
			# eta equation
			A[2*N2+j,j-1] = c1[j,i];		# u[j] (the J=j-1 index of A(,J) corresponds to u[j], u without the deadpoints)
			A[2*N2+j,N2+j-1] = c2[j];		# v[j]
			A[2*N2+j,N2+j] = c3[j];			# v[j+1]
			A[2*N2+j,N2+j-2] = - c3[j];		# v[j-1]
			A[2*N2+j,2*N2+j] = c4[j,i];		# eta[j]
	
		for j in range(0,N2):	
			F[j] = Ftilde1_nd[j+1,i];		# Have to add 1 to y-index to skip the first 'dead' gridpoint.
			F[N2+j] = Ftilde2_nd[j+1,i];
		for j in range(0,N):
			F[2*N2+j] = Ftilde3_nd[j,i];	

		solution[:,i] = np.linalg.solve(A,F);

		#import matplotlib.pyplot as plt
		#plt.plot(solution[0:N,i]);
		#plt.show();

	return solution;

#=======================================================

# FREE_SLIP_SOLVER
#=======================================================
def FREE_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,uBC,etaBC,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N,N2):
# Called by RSW_1L.py if BC = 'FREE-SLIP'.
# Currently uses the ghost point idea for dealing with boundary terms in the u eq, goes to 1-sided fd for v in eta eq

	dim = N2 + 2 * N;
	#print(dim);

	A = np.zeros((dim,dim),dtype=complex);	# For the free-slip, no-normal flow BC.

	# Initialise the forcing.
	F = np.zeros((dim),dtype=complex);

	# Initialise the solution.
	solution = np.zeros((dim,N),dtype=complex);

	for i in range(0,N):
		#print(i);
				
		# First the boundary terms. Some of these could be defined in the upcoming loop,
		# but for simplicity we define all boundary terms here.
	
		# u equation BCs
		# South
		A[0,0] = a1[0,i] - 2 * a2;		# u[0]		# See documentation for reasoning behind BCs here
		A[0,1] = a2;					# u[1]
		A[0,N2+N] = a4[i];				# eta[0] 
		# North
		A[N-1,N-1] = a1[N-1,i] - 2 * a2;	# u[N-1]
		A[N-1,N-2] = a2;					# u[N-2]
		A[N-1,2*N+N2-1] = a4[i];			# eta[N-1]

		# v equation BCs
		# South
		A[N,1] = f_nd[1];					# u[1]
		A[N,N] = a1[1,i] - 2. * a2;			# v[1]
		A[N,N+1] = a2;						# v[2]
		A[N,N+N2] = - b4;					# eta[0] 	
		A[N,N+N2+2] = b4;					# eta[2]
		# North
		A[N+N2-1,N2] = f_nd[N2];					# u[N-2] 
		A[N+N2-1,N+N2-1] = a1[N2,i] - 2. * a2;		# v[N-2]
		A[N+N2-1,N+N2-2] = a2;						# v[N-3]
		A[N+N2-1,2*N+N2-1] = b4;					# eta[N-1]
		A[N+N2-1,2*N+N2-3] = - b4;					# eta[N-3]

		# eta equation BCs (Here we have to define BCs at j=N+N2,N+N2+1,3*N-3,3*N-4)
		# South
		A[N+N2,0] = c1[0,i];				# u[0]		# Again, see documentation for BC reasoning
		A[N+N2,N] = 2. * c3[0];				# v[1] (factor of 2 because we use one-sided FD for v_y at the boundaries)
		A[N+N2,N+N2] = c4[0,i];				# eta[0]
		A[N+N2+1,1] = c1[1,i];				# u[1]
		A[N+N2+1,N] = c2[1];				# v[1]
		A[N+N2+1,N+1] = c3[1];				# v[2] (factor of 1 here, back to using centered FD)
		A[N+N2+1,N+N2+1] = c4[1,i];			# eta[1]
		# North
		A[2*N+N2-1,N-1] = c1[N-1,i];			# u[N-1]
		A[2*N+N2-1,N+N2-1] = - 2. * c3[N-1];	# v[N-2] (factor of 2 because we use one-sided FD for v_y at the boundaries)
		A[2*N+N2-1,2*N+N2-1] = c4[N-1,i];		# eta[N-1]
		A[N+2*N2,N2] = c1[N2,i];  				# u[N-2]
		A[N+2*N2,N+N2-1] = c2[N2];				# v[N-2]
		A[N+2*N2,N+N2-2] = - c3[N2];			# v[N-3] (factor of 1 here, back to using centered FD)
		A[N+2*N2,N+2*N2] = c4[N2,i];			# eta[N-2]

		# Now the inner values - two loops required: one for u,v, and one for eta
		for j in range(1,N-1):
			# u equation
			A[j,j] = a1[j,i] - 2. * a2;			# u[j+1]
			A[j,j-1] = a2;						# u[j]
			A[j,j+1] = a2;						# u[j+2]
			A[j,N+j-1] = a3[j];					# v[j]
			A[j,N+N2+j] = a4[i];				# eta[j] 
			
		for j in range(1,N2-1):
			# v equation
			A[N+j,j+1] = f_nd[j+1];					# u[j+1]
			A[N+j,N+j] = a1[j+1,i] - 2. * a2;		# v[j+1]
			A[N+j,N+j-1] = a2;						# v[j-1]
			A[N+j,N+j+1] = a2;						# v[j+1]
			A[N+j,N+N2+j] = - b4;					# eta[j-1]
			A[N+j,N+N2+j+2] = b4;					# eta[j+1]

		for j in range(2,N-2):
			# eta equation
			A[N2+N+j,j] = c1[j,i];			# u[j] 
			A[N2+N+j,N+j-1] = c2[j];		# v[j]
			A[N2+N+j,N+j] = c3[j];			# v[j+1]
			A[N2+N+j,N+j-2] = - c3[j];		# v[j-1]
			A[N2+N+j,N+N2+j] = c4[j,i];		# eta[j]
	
		# Now define the forcing from the forcing input file.
		for j in range(0,N2):
			F[N+j] = Ftilde2_nd[j+1,i];
		for j in range(0,N):
			F[j] = Ftilde1_nd[j,i];
			F[N2+N+j] = Ftilde3_nd[j,i];
			
		solution[:,i] = np.linalg.solve(A,F);
		
		#import matplotlib.pyplot as plt
		#plt.plot(solution[2*N:3*N,i]);
		#plt.show();
	
	return solution;

#=======================================================

# FREE_SLIP_SOLVER
#=======================================================
def FREE_SLIP_SOLVER2(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N,N2):
# Called by RSW_1L.py if BC = 'FREE-SLIP'.

	dim = N2 + 2 * N;
	#print(dim);
	
	A = np.zeros((dim,dim),dtype=complex);	# For the free-slip, no-normal flow BC.

	# Initialise the forcing.
	F = np.zeros((dim),dtype=complex);

	# Initialise the solution.
	solution = np.zeros((dim,N),dtype=complex);

	for i in range(0,N):
		#print(i);
				
		# First the boundary terms. Some of these could be defined in the upcoming loop,
		# but for simplicity we define all boundary terms here.
	
		# u equation BCs
		# South
		A[0,0] = a1[0,i]- 2 * a2;		# u[0]		# See documentation for reasoning behind BCs here
		A[0,1] = a2;					# u[1]
		#A[0,2]=a2;
		A[0,N2+N] = a4[i];				# eta[0] 
		# North
		A[N-1,N-1] = a1[N-1,i]- 2 * a2;	# u[N-1]
		A[N-1,N-2] = a2;					# u[N-2]
		#A[N-1,N-3]=a2;
		A[N-1,2*N+N2-1] = a4[i];			# eta[N-1]

		# v equation BCs
		# South
		A[N,1] = f_nd[1];					# u[1]
		A[N,N] = a1[1,i] - 2. * a2;			# v[1]
		A[N,N+1] = a2;						# v[2]
		A[N,N+N2] = - b4;					# eta[0] 	
		A[N,N+N2+2] = b4;					# eta[2]
		# North
		A[N+N2-1,N2] = f_nd[N2];					# u[N-2] 
		A[N+N2-1,N+N2-1] = a1[N2,i] - 2. * a2;		# v[N-2]
		A[N+N2-1,N+N2-2] = a2;						# v[N-3]
		A[N+N2-1,2*N+N2-1] = b4;					# eta[N-1]
		A[N+N2-1,2*N+N2-3] = - b4;					# eta[N-3]

		# eta equation BCs (Here we have to define BCs at j=N+N2,N+N2+1,3*N-3,3*N-4)
		# South
		A[N+N2,0] = c1[0,i];				# u[0]		# Again, see documentation for BC reasoning
		A[N+N2,N] = c3[0];				# v[1] (factor of 2 because we use one-sided FD for v_y at the boundaries)
		A[N+N2,N+N2] = c4[0,i];				# eta[0]
		A[N+N2+1,1] = c1[1,i];				# u[1]
		A[N+N2+1,N] = c2[1];				# v[1]
		A[N+N2+1,N+1] = c3[1];				# v[2] (factor of 1 here, back to using centered FD)
		A[N+N2+1,N+N2+1] = c4[1,i];			# eta[1]
		# North
		A[2*N+N2-1,N-1] = c1[N-1,i];			# u[N-1]
		A[2*N+N2-1,N+N2-1] = - c3[N-1];	# v[N-2] (factor of 2 because we use one-sided FD for v_y at the boundaries)
		A[2*N+N2-1,2*N+N2-1] = c4[N-1,i];		# eta[N-1]
		A[N+2*N2,N2] = c1[N2,i];  				# u[N-2]
		A[N+2*N2,N+N2-1] = c2[N2];				# v[N-2]
		A[N+2*N2,N+N2-2] = - c3[N2];			# v[N-3] (factor of 1 here, back to using centered FD)
		A[N+2*N2,N+2*N2] = c4[N2,i];			# eta[N-2]

		# Now the inner values - two loops required: one for u,v, and one for eta
		for j in range(1,N-1):
			# u equation
			A[j,j] = a1[j,i] - 2. * a2;			# u[j+1]
			A[j,j-1] = a2;						# u[j]
			A[j,j+1] = a2;						# u[j+2]
			A[j,N+j-1] = a3[j];					# v[j]
			A[j,N+N2+j] = a4[i];				# eta[j] 
			
		for j in range(1,N2-1):
			# v equation
			A[N+j,j+1] = f_nd[j+1];					# u[j+1]
			A[N+j,N+j] = a1[j+1,i] - 2. * a2;		# v[j+1]
			A[N+j,N+j-1] = a2;						# v[j-1]
			A[N+j,N+j+1] = a2;						# v[j+1]
			A[N+j,N+N2+j] = - b4;					# eta[j-1]
			A[N+j,N+N2+j+2] = b4;					# eta[j+1]

		for j in range(2,N-2):
			# eta equation
			A[N2+N+j,j] = c1[j,i];			# u[j] 
			A[N2+N+j,N+j-1] = c2[j];		# v[j]
			A[N2+N+j,N+j] = c3[j];			# v[j+1]
			A[N2+N+j,N+j-2] = - c3[j];		# v[j-1]
			A[N2+N+j,N+N2+j] = c4[j,i];		# eta[j]
	
		# Now define the forcing from the forcing input file.
		for j in range(0,N2):
			F[N+j] = Ftilde2_nd[j+1,i];
		for j in range(0,N):
			F[j] = Ftilde1_nd[j,i];
			F[N2+N+j] = Ftilde3_nd[j,i];
			
		solution[:,i] = np.linalg.solve(A,F);
		
		#import matplotlib.pyplot as plt
		#plt.plot(solution[2*N:3*N,i]);
		#plt.show();
	
	return solution;

# FREE_SLIP_SOLVER
#=======================================================
def FREE_SLIP_SOLVER3(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,uBC,etaBC,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N,N2):
# Called by RSW_1L.py if BC = 'FREE-SLIP'.

	dim = N2 + 2 * N;
	#print(dim);

	A = np.zeros((dim,dim),dtype=complex);	# For the free-slip, no-normal flow BC.

	# Initialise the forcing.
	F = np.zeros((dim),dtype=complex);

	# Initialise the solution.
	solution = np.zeros((dim,N),dtype=complex);

	for i in range(0,N):
		#print(i);
				
		# First the boundary terms. Some of these could be defined in the upcoming loop,
		# but for simplicity we define all boundary terms here.
	
		# u equation BCs
		# South
		A[0,0] = a1[0,i] + uBC;			# u[0]		# See documentation for reasoning behind BCs here
		A[0,1] = - uBC;					# u[1]
		A[0,N2+N] = a4[i];				# eta[0] 
		# North
		A[N-1,N-1] = a1[N-1,i] + uBC;	# u[N-1]
		A[N-1,N-2] = - uBC;				# u[N-2]
		A[N-1,2*N+N2-1] = a4[i];		# eta[N-1]

		# v equation BCs
		# South
		A[N,1] = f_nd[1];					# u[1]
		A[N,N] = a1[1,i] - 2. * a2;			# v[1]
		A[N,N+1] = a2;						# v[2]
		A[N,N+N2] = - b4;					# eta[0] 	
		A[N,N+N2+2] = b4;					# eta[2]
		# North
		A[N+N2-1,N2] = f_nd[N2];					# u[N-2] 
		A[N+N2-1,N+N2-1] = a1[N2,i] - 2. * a2;		# v[N-2]
		A[N+N2-1,N+N2-2] = a2;						# v[N-3]
		A[N+N2-1,2*N+N2-1] = b4;					# eta[N-1]
		A[N+N2-1,2*N+N2-3] = - b4;					# eta[N-3]

		# eta equation BCs (Here we have to define BCs at j=N+N2,N+N2+1,3*N-3,3*N-4)
		# South
		A[N+N2,0] = c1[0,i] + etaBC[0];		# u[0]		# Again, see documentation for BC reasoning
		A[N+N2,N] = 2. * c3[0];				# v[1] (factor of 2 because we use one-sided FD for v_y at the boundaries)
		A[N+N2,N+N2] = c4[0,i];				# eta[0]
		A[N+N2+1,1] = c1[1,i];				# u[1]
		A[N+N2+1,N] = c2[1];				# v[1]
		A[N+N2+1,N+1] = c3[1];				# v[2] (factor of 1 here, back to using centered FD)
		A[N+N2+1,N+N2+1] = c4[1,i];			# eta[1]
		# North
		A[2*N+N2-1,N-1] = c1[N-1,i] + etaBC[1];	# u[N-1]
		A[2*N+N2-1,N+N2-1] = - 2. * c3[N-1];	# v[N-2] (factor of 2 because we use one-sided FD for v_y at the boundaries)
		A[2*N+N2-1,2*N+N2-1] = c4[N-1,i];		# eta[N-1]
		A[N+2*N2,N2] = c1[N2,i];  				# u[N-2]
		A[N+2*N2,N+N2-1] = c2[N2];				# v[N-2]
		A[N+2*N2,N+N2-2] = - c3[N2];			# v[N-3] (factor of 1 here, back to using centered FD)
		A[N+2*N2,N+2*N2] = c4[N2,i];			# eta[N-2]

		# Now the inner values - two loops required: one for u,v, and one for eta
		for j in range(1,N-1):
			# u equation
			A[j,j] = a1[j,i] - 2. * a2;			# u[j+1]
			A[j,j-1] = a2;						# u[j]
			A[j,j+1] = a2;						# u[j+2]
			A[j,N+j-1] = a3[j];					# v[j]
			A[j,N+N2+j] = a4[i];				# eta[j] 
			
		for j in range(1,N2-1):
			# v equation
			A[N+j,j+1] = f_nd[j+1];					# u[j+1]
			A[N+j,N+j] = a1[j+1,i] - 2. * a2;		# v[j+1]
			A[N+j,N+j-1] = a2;						# v[j-1]
			A[N+j,N+j+1] = a2;						# v[j+1]
			A[N+j,N+N2+j] = - b4;					# eta[j-1]
			A[N+j,N+N2+j+2] = b4;					# eta[j+1]

		for j in range(2,N-2):
			# eta equation
			A[N2+N+j,j] = c1[j,i];			# u[j] 
			A[N2+N+j,N+j-1] = c2[j];		# v[j]
			A[N2+N+j,N+j] = c3[j];			# v[j+1]
			A[N2+N+j,N+j-2] = - c3[j];		# v[j-1]
			A[N2+N+j,N+N2+j] = c4[j,i];		# eta[j]
	
		# Now define the forcing from the forcing input file.
		for j in range(0,N2):
			F[N+j] = Ftilde2_nd[j+1,i];
		for j in range(0,N):
			F[j] = Ftilde1_nd[j,i];
			F[N2+N+j] = Ftilde3_nd[j,i];
			
		solution[:,i] = np.linalg.solve(A,F);

		#import matplotlib.pyplot as plt
		#plt.plot(solution[0:N,i]);
		#plt.show();
	
	return solution;

#=======================================================

# extractSols
#=======================================================	
def extractSols(solution,N,N2,BC):

	# Intialise the solutions
	utilde_nd = np.zeros((N,N),dtype=complex);
	vtilde_nd = np.zeros((N,N),dtype=complex);
	etatilde_nd = np.zeros((N,N),dtype=complex);
    
	if BC == 'NO-SLIP':
		for j in range(0,N2):
			utilde_nd[j+1,:] = solution[j,:];
			vtilde_nd[j+1,:] = solution[N2+j,:];
		for j in range(0,N):		
			etatilde_nd[j,:] = solution[2*N2+j,:];
		
	elif BC == 'FREE-SLIP':
		for j in range(0,N):
			utilde_nd[j,:] = solution[j,:];
			etatilde_nd[j,:] = solution[N+N2+j,:];
		for j in range(0,N2):
			vtilde_nd[j+1,:] = solution[N+j,:];
		
	else:
		print('ERROR');

	return utilde_nd, vtilde_nd, etatilde_nd;

#=======================================================

# SPEC_TO_PHYS
def SPEC_TO_PHYS(utilde,vtilde,etatilde,T,dx,omega,N):
# Function takes the spectral-physical solutions produced by the solver and returns the time-dependent solutions in physical space.
	
	I = np.complex(0,1);

	Nt = np.size(T) - 1;

	u = np.zeros((N,N,Nt),dtype=complex);
	v = np.zeros((N,N,Nt),dtype=complex);
	eta = np.zeros((N,N,Nt),dtype=complex);	

	for ti in range(0,Nt):
		# Calculate the solutions in physical space at some instant t. Solutions are divided (later) by extra factor AmpF, so that they are normalised by the forcing amplitude.
		u[:,:,ti] = np.exp(2.*np.pi*I*omega*T[ti])*np.fft.ifft(utilde,axis=1) / dx; 
		v[:,:,ti] = np.exp(2.*np.pi*I*omega*T[ti])*np.fft.ifft(vtilde,axis=1) / dx;		
		eta[:,:,ti] = np.exp(2.*np.pi*I*omega*T[ti])*np.fft.ifft(etatilde,axis=1) / dx;

	return u, v, eta;







