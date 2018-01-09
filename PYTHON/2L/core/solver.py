# solver.py
#=======================================================

# Contains a set of functions to be called by RSW_2L_visc.py. 
# They constitute the solvers (for two choices of BCs) as well as a function that defines all parameters of the linear solver.

#=======================================================

import numpy as np
from diagnostics import diff

#=======================================================

# SOLVER_COEFFICIENTS
#=======================================================
def SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U1_nd,U2_nd,H1_nd,H2_nd,rho1_nd,rho2_nd,omega_nd,gamma_nd,dy_nd,N):
# Here we also divide coeficients of derivatives by the relevant space-step: 
# 2*dy_nd for first-order derivatives, dy_nd**2 for second-order derivatives.
# Note: to speed up the algorithm, some coefficients aren't unnecessarily defined as they have a duplicate featuring in another equation.


	I = np.complex(0,1);

	# 1 - Define all k- and y-dependent coefficients
	delta1 = np.zeros((N,N));
	delta2 = np.zeros((N,N));
	a1 = np.zeros((N,N),dtype=complex);
	c1 = np.zeros((N,N),dtype=complex);
	c4 = np.zeros((N,N),dtype=complex);
	c5 = np.zeros((N,N),dtype=complex);
	d1 = np.zeros((N,N),dtype=complex);
	f1 = np.zeros((N,N),dtype=complex);
	f4 = np.zeros((N,N),dtype=complex);
	for j in range(0,N):
		for i in range(0,N):
			delta1[j,i] = 2. * np.pi * (omega_nd + U1_nd[j] * K_nd[i]);
			delta2[j,i] = 2. * np.pi * (omega_nd + U2_nd[j] * K_nd[i]);
			a1[j,i] = I * delta1[j,i] * Ro + 4 * np.pi**2 * K_nd[i]**2 * Ro / Re;
			c1[j,i] = 2. * np.pi * I * H1_nd[j] * K_nd[i]
			c4[j,i] = I * delta1[j,i];
			c5[j,i] = - I * delta1[j,i];
			d1[j,i] = I * delta2[j,i] * Ro + 4 * np.pi**2 * K_nd[i]**2 * Ro / Re + gamma_nd;
			f1[j,i] = 2. * np.pi * I * K_nd[i] * H2_nd[j];
			f4[j,i] = I * delta2[j,i];

	# 2 - Define all y-dependent coefficients
	a3 = Ro * diff(U1_nd,2,0,dy_nd) - f_nd;
	b1 = f_nd;
	d3 = Ro * diff(U2_nd,2,0,dy_nd) - f_nd;
	f2 = diff(H2_nd,2,0,dy_nd);
	f3 = H2_nd / (2 * dy_nd);

	# 3 - Define all k-dependent coefficients
	a4 = 2. * np.pi * I * K_nd;
	d4 = 2. * np.pi * I * K_nd * rho1_nd;
	d5 = 2. * np.pi * I * K_nd * rho2_nd;

	# 4 - Define all constant coefficients
	a2 = - Ro / (Re * dy_nd**2);
	b4 = 1. / (2 * dy_nd); 
	c2 = diff(H1_nd,2,0,dy_nd);
	c3 = H1_nd / (2 * dy_nd);
	e4 = rho1_nd / (2 * dy_nd);
	e5 = rho2_nd / (2 * dy_nd);

	# Summary of coefficients not explicitly defined for sake of algorithm speed:
	# b2 = a1, b3 = a2, d2 = a2, e1 = f_nd, e2 = d1, e3 = a2.
	
	return a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,c5,d1,d3,d4,d5,e4,e5,f1,f2,f3,f4;

#=======================================================

# NO_SLIP_SOLVER 
#=======================================================
def NO_SLIP_SOLVER(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,c5,d1,d3,d4,d5,e4,e5,f1,f2,f3,f4,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,Ftilde4_nd,Ftilde5_nd,Ftilde6_nd,N,N2):
# Called by RSW_2L.py if BC = 'NO-SLIP'.

	dim = 6 * N - 8; 	 # u and v have N-2 gridpoints, eta has N gridpoints
	#print(dim);

	A = np.zeros((dim,dim),dtype=complex);	# For the free-slip, no-normal flow BC.
	# eta has N gridpoints in y, whereas u and v have N2=N-2, after removing the two 'dead' gridpoints.
	# We primarily consider forcing away from the boundaries so that it's okay applying no-slip BCs.
	# We define a so that the 6 equations are ordered as follows: u1, u2, v1, v2, eta0, eta1.

	solution = np.zeros((dim,N),dtype=complex);	

	for i in range(0,N):
		#print(i);
			
		# Boundary Conditions

		# u1 equation
		# South
		A[0,0] = a1[1,i] - 2 * a2;		# u1[1]
		A[0,1] = a2;					# u1[2]
		A[0,2*N2] = a3[1];				# v1[1]
		A[0,4*N2+1] = a4[i];			# eta0[1]
		# North
		A[N2-1,N2-1] = a1[N2,i] - 2 * a2;	# u1[N-2]
		A[N2-1,N2-2] = a2;					# u1[N-3]
		A[N2-1,3*N2-1] = a3[N2];			# v1[N-3]
		A[N2-1,4*N2+N-2] = a4[i];			# eta0[N-2]

		# u2 equation
		# South
		A[N2,N2] = d1[1,i] - 2 * a2;	# u2[1]
		A[N2,N2+1] = a2;				# u2[2]
		A[N2,3*N2] = d3[1];				# v2[1]
		A[N2,4*N2+1] = d4[i];			# eta0[1]
		A[N2,4*N2+N+1] = d5[i];			# eta1[1]
		# North
		A[2*N2-1,2*N2-1] = d1[N2,i] - 2 * a2;	# u2[N-2]
		A[2*N2-1,2*N2-2] = a2;					# u2[N-2]
		A[2*N2-1,4*N2-1] = d3[N2];				# v2[N-2]
		A[2*N2-1,4*N2+N-2] = d4[i];				# eta0[N-2]
		A[2*N2-1,4*N2+2*N-2] = d5[i];			# eta1[N-2]
		
		# v1 equation
		# South
		A[2*N2,0] = b1[1];					# u1[1]
		A[2*N2,2*N2] = a1[1,i] - 2 * a2;	# v1[1]
		A[2*N2,2*N2+1] = a2;				# v1[2]
		A[2*N2,4*N2+2] = b4;				# eta0[2]
		A[2*N2,4*N2] = - b4;				# eta0[0]
		# North
		A[3*N2-1,N2-1] = b1[N2];				# u1[N-2]
		A[3*N2-1,3*N2-1] = a1[N2,i] - 2 * a2;	# v1[N-2]
		A[3*N2-1,3*N2-2] = a2;					# v1[N-3]
		A[3*N2-1,4*N2+N-1] = b4;				# eta0[N-1]
		A[3*N2-1,4*N2+N-3] = - b4;				# eta0[N-3]

		# v2 equation
		# South
		A[3*N2,N2] = b1[1];					# u2[1]
		A[3*N2,3*N2] = d1[1,i] - 2 * a2;	# v2[1]
		A[3*N2,3*N2+1] = a2;				# v2[2]
		A[3*N2,4*N2+2] = e4;				# eta0[2]
		A[3*N2,4*N2] = - e4;				# eta0[0]
		A[3*N2,4*N2+N+2] = e5;				# eta1[2]
		A[3*N2,4*N2+N] = - e5;				# eta1[0]
		# North
		A[4*N2-1,2*N2-1] = b1[N2];				# u2[N-2]
		A[4*N2-1,4*N2-1] = d1[N2,i] - 2 * a2;	# v2[N-2]
		A[4*N2-1,4*N2-2] = a2;					# v2[N-3]
		A[4*N2-1,4*N2+N-1] = e4;				# eta0[N-1]
		A[4*N2-1,4*N2+N-3] = - e4;				# eta0[N-3]
		A[4*N2-1,4*N2+2*N-1] = e5;				# eta1[N-1]
		A[4*N2-1,4*N2+2*N-3] = - e5;			# eta1[N-3]
	
		# h1 equation
		# South
		A[4*N2,2*N2] = 2 * c3[0];		# v1[1] (one-sided FD)
		A[4*N2,4*N2] = c4[0,i];			# eta0[0]
		A[4*N2,4*N2+N] = c5[0,i];		# eta1[0]
		A[4*N2+1,0] = c1[1,i];			# u1[1]
		A[4*N2+1,2*N2] = c2[1];			# v1[1]
		A[4*N2+1,2*N2+1] = c3[1];		# v1[2]
		A[4*N2+1,4*N2+1] = c4[1,i];		# eta0[1]
		A[4*N2+1,4*N2+N+1] = c5[1,i];	# eta1[1]
		# North
		A[4*N2+N-1,3*N2-1] = - 2 * c3[N-1];		# v1[N-2]
		A[4*N2+N-1,4*N2+N-1] = c4[N-1,i];		# eta0[N-1]
		A[4*N2+N-1,4*N2+2*N-1] = c5[N-1,i];		# eta1[N-1]
		A[4*N2+N-2,N2-1] = c1[N2,i];			# u1[N-2]
		A[4*N2+N-2,3*N2-1] = c2[N2];			# v1[N-2]
		A[4*N2+N-2,3*N2-2] = - c3[N2];			# v1[N-3]
		A[4*N2+N-2,4*N2+N-2] = c4[N2,i];		# eta0[N-2]
		A[4*N2+N-2,4*N2+2*N-2] = c5[N2,i];		# eta1[N-2]

		# h2 equation
		# South
		A[4*N2+N,3*N2] = 2 * f3[0];			# v2[1] (one-sided FD)
		A[4*N2+N,4*N2+N] = f4[0,i];			# eta1[0]
		A[4*N2+N+1,N2] = f1[1,i];			# u2[1]
		A[4*N2+N+1,3*N2] = f2[1];			# v2[1]
		A[4*N2+N+1,3*N2+1] = f3[1];			# v2[2]
		A[4*N2+N+1,4*N2+N+1] = f4[1,i];		# eta1[1]
		# North
		A[4*N2+2*N-1,4*N2-1] = - 2 * f3[N-1];	# v2[N-1]
		A[4*N2+2*N-1,4*N2+2*N-1] = f4[N-1,i];	# eta1[N-1]
		A[4*N2+2*N-2,2*N2-1] = f1[N2,i];		# u2[N-2]
		A[4*N2+2*N-2,4*N2-1] = f2[N2];			# v2[N-2]
		A[4*N2+2*N-2,4*N2-2] = f3[N2];			# v2[N-3]
		A[4*N2+2*N-2,4*N2+2*N-2] = f4[N2,i];	# eta1[N-2]		

		# Inner domain values

		# A loop for the u and v equations
		for j in range(1,N-3):
			# u1 equation
			A[j,j] = a1[j+1,i] - 2 * a2;	# u1[j+1]
			A[j,j+1] = a2;					# u1[j+2]
			A[j,j-1] = a2;					# u1[j]
			A[j,2*N2+j] = a3[j+1];			# v1[j+1] 
			A[j,4*N2+j+1] = a4[i];			# eta0[j+1]

			# u2 equation
			A[N2+j,N2+j] = d1[j+1,i] - 2 * a2;	# u2[j+1]
			A[N2+j,N2+j+1] = a2;				# u2[j+2]
			A[N2+j,N2+j-1] = a2;				# u2[j]
			A[N2+j,3*N2+j] = d3[j+1];			# v2[j+1]
			A[N2+j,4*N2+j+1] = d4[i];			# eta0[j+1]
			A[N2+j,4*N2+N+j+1] = d5[i];			# eta1[j+1]
	
			# v1 equation
			A[2*N2+j,j] =	b1[j+1];				# u1[j+1]
			A[2*N2+j,2*N2+j] = a1[j+1,i] - 2 * a2;	# v1[j+1]
			A[2*N2+j,2*N2+j+1] = a2;				# v1[j+2]
			A[2*N2+j,2*N2+j-1] = a2;				# v1[j]
			A[2*N2+j,4*N2+j+2] = b4;				# eta0[j+2]
			A[2*N2+j,4*N2+j] = - b4;				# et0[j] 

			# v2 equation
			A[3*N2+j,N2+j] = b1[j+1];				# u2[j+1]
			A[3*N2+j,3*N2+j] = d1[j+1,i] - 2 * a2;	# v2[j+1]
			A[3*N2+j,3*N2+j+1] = a2;				# v2[j+2]
			A[3*N2+j,3*N2+j-1] = a2;				# v2[j]
			A[3*N2+j,4*N2+j+2] = e4;				# eta0[j+2]
			A[3*N2+j,4*N2+j] = - e4;				# eta0[j]
			A[3*N2+j,4*N2+N+j+2] = e5;				# eta1[j+2]
			A[3*N2+j,4*N2+N+j] = - e5;				# eta1[j]			
			
		# A loop for the h equations

		for j in range(2,N-2):
			# h1 equation
			A[4*N2+j,j-1] = c1[j,i];			# u1[j]
			A[4*N2+j,2*N2+j-1] = c2[j];			# v1[j]
			A[4*N2+j,2*N2+j] = c3[j];			# v1[j+1]
			A[4*N2+j,2*N2+j-2] = - c3[j];		# v1[j-1]
			A[4*N2+j,4*N2+j] = c4[j,i];			# eta0[j]
			A[4*N2+j,4*N2+N+j] = c5[j,i];		# eta1[j] 

			# h2 equation
			A[4*N2+N+j,N2+j-1] = f1[j,i];		# u2[j]
			A[4*N2+N+j,3*N2+j-1] = f2[j];		# v2[j]
			A[4*N2+N+j,3*N2+j] = f3[j];			# v2[j+1]
			A[4*N2+N+j,3*N2+j-2] = - f3[j];		# v2[j-1]
			A[4*N2+N+j,4*N2+N+j] = f4[j,i];		# eta1[j]
			
		# Forcing
		F = np.zeros(dim,dtype=complex);
		# Now assign the values to F
		for j in range(0,N-2):
			F[j] = Ftilde1_nd[j+1,i];			# Forcing the u1 equation
			F[2*N2+j] = Ftilde2_nd[j+1,i];		# Forcing the v1 equation
		for j in range(0,N):	
			F[4*N2+j] = Ftilde3_nd[j,i];	# Forcing the h1 equation
			F[4*N2+N+j] = Ftilde6_nd[j,i];	# Forcing the h2 equation

		solution[:,i] = np.linalg.solve(A,F);

	u1tilde_nd = np.zeros((N,N),dtype=complex);
	u2tilde_nd = np.zeros((N,N),dtype=complex);
	v1tilde_nd = np.zeros((N,N),dtype=complex);
	v2tilde_nd = np.zeros((N,N),dtype=complex);
	eta0tilde_nd = np.zeros((N,N),dtype=complex);
	eta1tilde_nd = np.zeros((N,N),dtype=complex);
	for j in range(0,N2):
		u1tilde_nd[j+1,:] = solution[j,:];
		u2tilde_nd[j+1,:] = solution[N2+j,:];
		v1tilde_nd[j+1,:] = solution[2*N2+j,:];
		v2tilde_nd[j+1,:] = solution[3*N2+j,:];
	for j in range(0,N):
		eta0tilde_nd[j,:] = solution[4*N2+j,:];
		eta1tilde_nd[j,:] = solution[4*N2+N+j,:];	

	return u1tilde_nd, u2tilde_nd, v1tilde_nd, v2tilde_nd, eta0tilde_nd, eta1tilde_nd;

#=======================================================

# FREE_SLIP_SOLVER
#=======================================================
def FREE_SLIP_SOLVER(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,c5,d1,d3,d4,d5,e4,e5,f1,f2,f3,f4,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,Ftilde4_nd,Ftilde5_nd,Ftilde6_nd,N,N2):
# Called by RSW_2L.py if BC = 'FREE-SLIP'.
# The boundary conditions are: u_y = v = v_yy = 0.

	dim = 6 * N - 4; 	 # u and eta have N gridpoints, v have N-2 gridpoints
	#print(dim);

	A = np.zeros((dim,dim),dtype=complex);	# For the free-slip, no-normal flow BC.
	# u and eta have N gridpoints in y, whereas v has N2=N-2, after removing the two 'dead' gridpoints.
	# We define a so that the 6 equations are ordered as follows: u1, u2, v1, v2, eta0, eta1.

	solution = np.zeros((dim,N),dtype=complex);		

	for i in range(0,N):
		#print(i);

		# Boundary conditions

		# u1 equation
		# South
		A[0,0] = a1[0,i] - 2. * a2;	# u1[0]
		A[0,1] = 2. * a2;			# u1
		A[0,2*N+2*N2] = a4[i];		# eta0[0]
		# North
		A[N-1,N-1] = a1[N-1,i] - 2. * a2;	# u1[N-1]
		A[N-1,N-2] = 2. * a2;				# u1[N-2]
		A[N-1,3*N+2*N2-1] = a4[i];			# eta0[N-1]

		# u2 equation
		# South
		A[N,N] = d1[0,i] - 2. * a2;	# u2[0]
		A[N,N+1] = 2. * a2;			# u2[1]
		A[N,2*N+2*N2] = d4[i];		# eta0[0]
		A[N,3*N+2*N2] = d5[i];		# eta1[0]
		# North
		A[2*N-1,2*N-1] = d1[N-1,i] - 2. * a2;	# u2[N-1]
		A[2*N-1,2*N-2] = 2. * a2;				# u2[N-2]
		A[2*N-1,3*N+2*N2-1] = d4[i];			# eta0[N-1]
		A[2*N-1,4*N+2*N2-1] = d5[i];			# eta1[N-1]
		
		# edited up to here.....
		# v1 equation
		# South
		A[2*N,1] = f_nd[1];				# u1[1]
		A[2*N,2*N] = a1[1,i] - 2. * a2;	# v1[1]
		A[2*N,2*N+1] = a2;				# v1[2]
		A[2*N,2*N+2*N2] = - b4;			# eta0[0]
		A[2*N,2*N+2*N2+2] = b4;			# eta0[2]
		# North
		A[2*N+N2-1,N-2] = f_nd[N2];					# u1[N-2]
		A[2*N+N2-1,2*N+N2-1] = a1[N-2,i] - 2. * a2;	# v1[N-2]
		A[2*N+N2-1,2*N+N2-2] = a2;					# v1[N-3]
		A[2*N+N2-1,3*N+2*N2-1] = b4;				# eta0[N-1]
		A[2*N+N2-1,3*N+2*N2-3] = -b4;				# eta0[N-3]

		# v2 equation
		# South
		A[2*N+N2,N+1] = f_nd[1];				# u2[1]
		A[2*N+N2,2*N+N2] = d1[1,i] - 2. * a2;	# v2[1]
		A[2*N+N2,2*N+N2+1] = a2;				# v2[2]
		A[2*N+N2,2*N+2*N2+2] = e4;				# eta0[2]
		A[2*N+N2,2*N+2*N2] = - e4;				# eta0[0]
		A[2*N+N2,3*N+2*N2+2] = e5;				# eta1[2]
		A[2*N+N2,3*N+2*N2] = - e5;				# eta1[0]
		# North
		A[2*N+2*N2-1,2*N-2] = f_nd[N2];					# u2[N-2]
		A[2*N+2*N2-1,2*N+2*N2-1] = d1[N2,i] - 2. * a2;	# v2[N-2]
		A[2*N+2*N2-1,2*N+2*N2-2] = a2;					# v2[N-3]
		A[2*N+2*N2-1,3*N+2*N2-1] = e4;					# eta0[N-1]
		A[2*N+2*N2-1,3*N+2*N2-3] = - e4;				# eta0[N-3]
		A[2*N+2*N2-1,4*N+2*N2-1] = e5;					# eta1[N-1]
		A[2*N+2*N2-1,4*N+2*N2-3] = - e5;				# eta1[N-3]

		# h1 = eta0 - eta1 equation
		# South
		A[2*N+2*N2,0] = c1[0,i];				# u1[0]
		A[2*N+2*N2,2*N] = 2. * c3[0];			# v1[1] (one-sided FD approx.)
		A[2*N+2*N2,2*N+2*N2] = c4[0,i];			# eta0[0]
		A[2*N+2*N2,3*N+2*N2] = c5[0,i];			# eta1[0]
		A[2*N+2*N2+1,1] = c1[0,i];				# u1[1]
		A[2*N+2*N2+1,2*N] = c2[1];				# v1[1]
		A[2*N+2*N2+1,2*N+1] = c3[1];			# v1[2]
		A[2*N+2*N2+1,2*N+2*N2+1] = c4[1,i];		# eta0[1]
		A[2*N+2*N2+1,3*N+2*N2+1] = c5[1,i];		# eta1[1]
		# North
		A[3*N+2*N2-1,N-1] = c1[N-1,i];				# u1[N-1]
		A[3*N+2*N2-1,2*N+N2-1] = - 2. * c3[N-1];	# v1[N-2] (one-sided FD approx.)
		A[3*N+2*N2-1,3*N+2*N2-1] = c4[N-1,i];		# eta0[N-1]
		A[3*N+2*N2-1,4*N+2*N2-1] = c5[N-1,i];		# eta1[N-1]
		A[3*N+2*N2-2,N-2] = c1[N2,i];				# u1[N-2]
		A[3*N+2*N2-2,2*N+N2-1] = c2[N2];			# v1[N-2]
		A[3*N+2*N2-2,2*N+N2-2] = - c3[N2];			# v1[N-3]	
		A[3*N+2*N2-2,3*N+2*N2-2] = c4[N2,i];		# eta0[N-2]
		A[3*N+2*N2-2,4*N+2*N2-2] = c5[N2,i];		# eta1[N-2]

		# h2 = eta1 equation
		# South
		A[3*N+2*N2,N] = f1[0,i];				# u2[0]
		A[3*N+2*N2,2*N+N2] = 2. * f3[0];		# v2[0] (one-sided FD approx.)
		A[3*N+2*N2,3*N+2*N2] = f4[0,i];			# eta1[0]
		A[3*N+2*N2+1,N+1] = f1[1,i];			# u2[1]
		A[3*N+2*N2+1,2*N+N2] = f2[1];			# v2[1]
		A[3*N+2*N2+1,2*N+N2] = f3[1];			# v2[2]
		A[3*N+2*N2+1,3*N+2*N2+1] = f4[1,i];		# eta1[1]
		# North
		A[4*N+2*N2-1,N-1] = f1[N-1,i];					# u2[N-1]
		A[4*N+2*N2-1,2*N+2*N2-1] = - 2. * f3[N-1];		# v2[N-1] (one-sided FD approx.)
		A[4*N+2*N2-1,4*N+2*N2-1] = f4[N-1,i];			# eta1[N-1]
		A[4*N+2*N2-2,N-2] = f1[N2,i];					# u2[N-2]
		A[4*N+2*N2-2,2*N+2*N2-1] = f2[N2];				# v2[N-2]
		A[4*N+2*N2-2,2*N+2*N2-2] = - f3[N2];			# v2[N-3]
		A[4*N+2*N2-2,4*N+2*N2-2] = f4[N2,i];			# eta1[N-2]

		# Inner domain values	
	
		# A loop for remaining values of the u equations.
		for j in range(1,N-1):
			# u1 equation
			A[j,j] = a1[j,i] - 2. * a2;		# u1[j]
			A[j,j+1] = a2;					# u1[j+1]
			A[j,j-1] = a2;					# u1[j-1]
			A[j,2*N+j-1] = a3[j];			# v1[j]
			A[j,2*N+2*N2+j] = a4[i];		# eta0[j]
			# u2 equation
			A[N+j,N+j] = d1[j,i] - 2. * a2;		# u2[j]
			A[N+j,N+j+1] = a2;					# u2[j+1]
			A[N+j,N+j-1] = a2;					# u2[j-1]
			A[N+j,2*N+N2+j-1] = d3[j];			# v2[j]
			A[N+j,2*N+2*N2+j] = d4[i];			# eta0[j]
			A[N+j,3*N+2*N2+j] = d5[i];			# eta1[j]	

		# A loop for the remaining values of the v equations.
		for j in range(1,N-3):
			# v1 equation
			A[2*N+j,j+1] = f_nd[j+1];				# u1[j+1]
			A[2*N+j,2*N+j] = a1[j+1,i] - 2. * a2;	# v1[j+1]
			A[2*N+j,2*N+j+1] = a2;					# v1[j+2]
			A[2*N+j,2*N+j-1] = a2;					# v1[j]
			A[2*N+j,2*N+2*N2+j+2] = b4;				# eta0[j+2]
			A[2*N+j,2*N+2*N2+j] = - b4;				# eta0[j]
			# v2 equation
			A[2*N+N2+j,N+j+1] = f_nd[j+1];				# u2[j+1]
			A[2*N+N2+j,2*N+N2+j] = d1[j+1,i] - 2. * a2;	# v2[j+1]
			A[2*N+N2+j,2*N+N2+j+1] = a2;				# v2[j+2]
			A[2*N+N2+j,2*N+N2+j-1] = a2;				# v2[j]
			A[2*N+N2+j,2*N+2*N2+j+2] = e4;				# eta0[j+2]
			A[2*N+N2+j,2*N+2*N2+j] = - e4;				# eta0[j]
			A[2*N+N2+j,3*N+2*N2+j+2] = e5;				# eta1[j+2]
			A[2*N+N2+j,3*N+2*N2+j] = - e5;				# eta1[j]

		# A loop for the remaining values of the h/eta equations.
		for j in range(2,N-2):
			# h1 equation
			A[2*N+2*N2+j,j] = c1[j,i];				# u1[j]
			A[2*N+2*N2+j,2*N+j-1] = c2[j];			# v1[j]
			A[2*N+2*N2+j,2*N+j] = c3[j];			# v1[j+1]
			A[2*N+2*N2+j,2*N+j-2] = - c3[j];		# v1[j-1]
			A[2*N+2*N2+j,2*N+2*N2+j] = c4[j,i];		# eta0[j]
			A[2*N+2*N2+j,3*N+2*N2+j] = c5[j,i];		# eta1[j]
			# h2 equation
			A[3*N+2*N2+j,N+j] = f1[j,i];				# u2[j]
			A[3*N+2*N2+j,2*N+N2+j-1] = f2[j];			# v2[j]
			A[3*N+2*N2+j,2*N+N2+j] = f3[j];				# v2[j+1]
			A[3*N+2*N2+j,2*N+N2+j-2] = - f3[j];			# v2[j-1]
			A[3*N+2*N2+j,3*N+2*N2+j] = f4[j,i];			# eta1[j]
		
		# Now assign the values to F
		for j in range(0,N):	
			F[j] = Ftilde1_nd[j,i];				# Forcing the u1 equation
			F[2*N+2*N2+j] = Ftilde3_nd[j,i];	# Forcing the h1 equation
			F[3*N+2*N2+j] = Ftilde6_nd[j,i];	# Forcing the h2 equation
		for j in range(0,N-2):
			F[2*N+j] = Ftilde2_nd[j+1,i];		# Forcing the v1 equation
		

		solution[:,i] = np.linalg.solve(A,F);

	u1tilde_nd = np.zeros((N,N),dtype=complex);
	u2tilde_nd = np.zeros((N,N),dtype=complex);
	v1tilde_nd = np.zeros((N,N),dtype=complex);
	v2tilde_nd = np.zeros((N,N),dtype=complex);
	eta0tilde_nd = np.zeros((N,N),dtype=complex);
	eta1tilde_nd = np.zeros((N,N),dtype=complex);
	for j in range(0,N):
		u1tilde_nd[j,:] = solution[j,:];
		u2tilde_nd[j,:] = solution[N+j,:];
		eta0tilde_nd[j,:] = solution[2*N+2*N2+j];
		eta1tilde_nd[j,:] = solution[3*N+2*N2+j];		
	for j in range(0,N2):
		v1tilde_nd[j+1,:] = solution[2*N+j,:];
		v2tilde_nd[j+1,:] = solution[2*N+N2+j,:];

	return u1tilde_nd, u2tilde_nd, v1tilde_nd, v2tilde_nd, eta0tilde_nd, eta1tilde_nd;

#=======================================================


# FREE_SLIP_SOLVER2
#=======================================================
def FREE_SLIP_SOLVER2(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,c5,d1,d3,d4,d5,e4,e5,f1,f2,f3,f4,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,Ftilde4_nd,Ftilde5_nd,Ftilde6_nd,N,N2):
# Called by RSW_2L.py if BC = 'FREE-SLIP'.

	dim = 6 * N - 4; 	 # u and eta have N gridpoints, v have N-2 gridpoints
	#print(dim);

	A = np.zeros((dim,dim),dtype=complex);	# For the free-slip, no-normal flow BC.
	# u and eta have N gridpoints in y, whereas v has N2=N-2, after removing the two 'dead' gridpoints.
	# We define a so that the 6 equations are ordered as follows: u1, u2, v1, v2, eta0, eta1.

	# Initialise the forcing.
	F = np.zeros((dim),dtype=complex);
	# Now define the forcing - these will be used later in the loop over i-values.
	F1_nd, F2_nd, F3_nd, F4_nd, F5_nd, F6_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd, Ftilde4_nd, Ftilde5_nd, Ftilde6_nd = forcing_2L.Forcing(x,y,K,y0,r0,N,FORCE1,AmpF,g,f,f0,U,L,rho1_nd,rho2_nd,dx,dy);

	solution = np.zeros((dim,N),dtype=complex);		

	for i in range(0,N):
		#print(i);

		# Boundary conditions

		# u1 equation
		# South
		A[0,0] = a1[0,i] + a2; 		# u1[0]
		A[0,1] = - 2. * a2;			# u1[1]
		A[0,2] = a2; 				# u1[2]
		A[0,2*N+2*N2] = a4[i];		# eta0[0]
		# North
		A[N-1,N-1] = a1[N-1,i] + a2;	# u1[N-1]
		A[N-1,N-2] = - 2. * a2;			# u1[N-2]
		A[N-1,N-3] = a2;				# u1[N-3]
		A[N-1,3*N+2*N2-1]				# eta0[N-1]

		# u2 equation
		# South
		A[N,N] = d1[0,i] + a2;		# u2[0]
		A[N,N+1] = - 2. * a2;		# u2[1]
		A[N,N+2] = a2;				# u2[2]
		A[N,2*N+2*N2] = d4[i];		# eta0[0]
		A[N,3*N+2*N2] = d5[i];		# eta1[0]
		# North
		A[2*N-1,2*N-1] = d1[N-1,i] + a2;	# u2[N-1]
		A[2*N-1,2*N-2] = - 2. * a2;			# u2[N-2]
		A[2*N-1,2*N-2] = a2;				# u2[N-3]
		A[2*N-1,3*N+2*N2-1] = d4[i];		# eta0[N-1]
		A[2*N-1,4*N+2*N2-1] = d5[i];		# eta1[N-1]
		
		# v1 equation
		# South
		A[2*N,1] = f_nd[1];				# u1[1]
		A[2*N,2*N] = a1[1,i] - 2. * a2;	# v1[1]
		A[2*N,2*N+1] = a2;				# v1[2]
		A[2*N,2*N+2*N2] = - b4;			# eta0[0]
		A[2*N,2*N+2*N2+2] = b4;			# eta0[2]
		# North
		A[2*N+N2-1,N-2] = f_nd[N2];					# u1[N-2]
		A[2*N+N2-1,2*N+N2-1] = a1[N-2,i] - 2. * a2;	# v1[N-2]
		A[2*N+N2-1,2*N+N2-2] = a2;					# v1[N-3]
		A[2*N+N2-1,3*N+2*N2-1] = b4;				# eta0[N-1]
		A[2*N+N2-1,3*N+2*N2-3] = -b4;				# eta0[N-3]

		# v2 equation
		# South
		A[2*N+N2,N+1] = f_nd[1];				# u2[1]
		A[2*N+N2,2*N+N2] = d1[1,i] - 2. * a2;	# v2[1]
		A[2*N+N2,2*N+N2+1] = a2;				# v2[2]
		A[2*N+N2,2*N+2*N2+2] = e4;				# eta0[2]
		A[2*N+N2,2*N+2*N2] = - e4;				# eta0[0]
		A[2*N+N2,3*N+2*N2+2] = e5;				# eta1[2]
		A[2*N+N2,3*N+2*N2] = - e5;				# eta1[0]
		# North
		A[2*N+2*N2-1,2*N-2] = f_nd[N2];					# u2[N-2]
		A[2*N+2*N2-1,2*N+2*N2-1] = d1[N2,i] - 2. * a2;	# v2[N-2]
		A[2*N+2*N2-1,2*N+2*N2-2] = a2;					# v2[N-3]
		A[2*N+2*N2-1,3*N+2*N2-1] = e4;					# eta0[N-1]
		A[2*N+2*N2-1,3*N+2*N2-3] = - e4;				# eta0[N-3]
		A[2*N+2*N2-1,4*N+2*N2-1] = e5;					# eta1[N-1]
		A[2*N+2*N2-1,4*N+2*N2-3] = - e5;				# eta1[N-3]

		# h1 = eta0 - eta1 equation
		# South
		A[2*N+2*N2,0] = c1[0,i];				# u1[0]
		A[2*N+2*N2,2*N] = 2. * c3[0];			# v1[1] (one-sided FD approx.)
		A[2*N+2*N2,2*N+2*N2] = c4[0,i];			# eta0[0]
		A[2*N+2*N2,3*N+2*N2] = c5[0,i];			# eta1[0]
		A[2*N+2*N2+1,1] = c1[0,i];				# u1[1]
		A[2*N+2*N2+1,2*N] = c2[1];				# v1[1]
		A[2*N+2*N2+1,2*N+1] = c3[1];			# v1[2]
		A[2*N+2*N2+1,2*N+2*N2+1] = c4[1,i];		# eta0[1]
		A[2*N+2*N2+1,3*N+2*N2+1] = c5[1,i];		# eta1[1]
		# North
		A[3*N+2*N2-1,N-1] = c1[N-1,i];				# u1[N-1]
		A[3*N+2*N2-1,2*N+N2-1] = - 2. * c3[N-1];	# v1[N-2] (one-sided FD approx.)
		A[3*N+2*N2-1,3*N+2*N2-1] = c4[N-1,i];		# eta0[N-1]
		A[3*N+2*N2-1,4*N+2*N2-1] = c5[N-1,i];		# eta1[N-1]
		A[3*N+2*N2-2,N-2] = c1[N2,i];				# u1[N-2]
		A[3*N+2*N2-2,2*N+N2-1] = c2[N2];			# v1[N-2]
		A[3*N+2*N2-2,2*N+N2-2] = - c3[N2];			# v1[N-3]	
		A[3*N+2*N2-2,3*N+2*N2-2] = c4[N2,i];		# eta0[N-2]
		A[3*N+2*N2-2,4*N+2*N2-2] = c5[N2,i];		# eta1[N-2]

		# h2 = eta1 equation
		# South
		A[3*N+2*N2,N] = f1[0,i];				# u2[0]
		A[3*N+2*N2,2*N+N2] = 2. * f3[0];		# v2[0] (one-sided FD approx.)
		A[3*N+2*N2,3*N+2*N2] = f4[0,i];			# eta1[0]
		A[3*N+2*N2+1,N+1] = f1[1,i];			# u2[1]
		A[3*N+2*N2+1,2*N+N2] = f2[1];			# v2[1]
		A[3*N+2*N2+1,2*N+N2] = f3[1];			# v2[2]
		A[3*N+2*N2+1,3*N+2*N2+1] = f4[1,i];		# eta1[1]
		# North
		A[4*N+2*N2-1,N-1] = f1[N-1,i];					# u2[N-1]
		A[4*N+2*N2-1,2*N+2*N2-1] = - 2. * f3[N-1];		# v2[N-1] (one-sided FD approx.)
		A[4*N+2*N2-1,4*N+2*N2-1] = f4[N-1,i];			# eta1[N-1]
		A[4*N+2*N2-2,N-2] = f1[N2,i];					# u2[N-2]
		A[4*N+2*N2-2,2*N+2*N2-1] = f2[N2];				# v2[N-2]
		A[4*N+2*N2-2,2*N+2*N2-2] = - f3[N2];			# v2[N-3]
		A[4*N+2*N2-2,4*N+2*N2-2] = f4[N2,i];			# eta1[N-2]

		# Inner domain values	
	
		# A loop for remaining values of the u equations.
		for j in range(1,N-1):
			# u1 equation
			A[j,j] = a1[j,i] - 2. * a2;		# u1[j]
			A[j,j+1] = a2;					# u1[j+1]
			A[j,j-1] = a2;					# u1[j-1]
			A[j,2*N+j-1] = a3[j];			# v1[j]
			A[j,2*N+2*N2+j] = a4[i];		# eta0[j]
			# u2 equation
			A[N+j,N+j] = d1[j,i] - 2. * a2;		# u2[j]
			A[N+j,N+j+1] = a2;					# u2[j+1]
			A[N+j,N+j-1] = a2;					# u2[j-1]
			A[N+j,2*N+N2+j-1] = d3[j];			# v2[j]
			A[N+j,2*N+2*N2+j] = d4[i];			# eta0[j]
			A[N+j,3*N+2*N2+j] = d5[i];			# eta1[j]	

		# A loop for the remaining values of the v equations.
		for j in range(1,N-3):
			# v1 equation
			A[2*N+j,j+1] = f_nd[j+1];				# u1[j+1]
			A[2*N+j,2*N+j] = a1[j+1,i] - 2. * a2;	# v1[j+1]
			A[2*N+j,2*N+j+1] = a2;					# v1[j+2]
			A[2*N+j,2*N+j-1] = a2;					# v1[j]
			A[2*N+j,2*N+2*N2+j+2] = b4;				# eta0[j+2]
			A[2*N+j,2*N+2*N2+j] = - b4;				# eta0[j]
			# v2 equation
			A[2*N+N2+j,N+j+1] = f_nd[j+1];				# u2[j+1]
			A[2*N+N2+j,2*N+N2+j] = d1[j+1,i] - 2. * a2;	# v2[j+1]
			A[2*N+N2+j,2*N+N2+j+1] = a2;				# v2[j+2]
			A[2*N+N2+j,2*N+N2+j-1] = a2;				# v2[j]
			A[2*N+N2+j,2*N+2*N2+j+2] = e4;				# eta0[j+2]
			A[2*N+N2+j,2*N+2*N2+j] = - e4;				# eta0[j]
			A[2*N+N2+j,3*N+2*N2+j+2] = e5;				# eta1[j+2]
			A[2*N+N2+j,3*N+2*N2+j] = - e5;				# eta1[j]

		# A loop for the remaining values of the h/eta equations.
		for j in range(2,N-2):
			# h1 equation
			A[2*N+2*N2+j,j] = c1[j,i];				# u1[j]
			A[2*N+2*N2+j,2*N+j-1] = c2[j];			# v1[j]
			A[2*N+2*N2+j,2*N+j] = c3[j];			# v1[j+1]
			A[2*N+2*N2+j,2*N+j-2] = - c3[j];		# v1[j-1]
			A[2*N+2*N2+j,2*N+2*N2+j] = c4[j,i];		# eta0[j]
			A[2*N+2*N2+j,3*N+2*N2+j] = c5[j,i];		# eta1[j]
			# h2 equation
			A[3*N+2*N2+j,N+j] = f1[j,i];				# u2[j]
			A[3*N+2*N2+j,2*N+N2+j-1] = f2[j];			# v2[j]
			A[3*N+2*N2+j,2*N+N2+j] = f3[j];				# v2[j+1]
			A[3*N+2*N2+j,2*N+N2+j-2] = - f3[j];			# v2[j-1]
			A[3*N+2*N2+j,3*N+2*N2+j] = f4[j,i];			# eta1[j]
		
		# Now assign the values to F
		for j in range(0,N):	
			F[j] = Ftilde1_nd[j,i];				# Forcing the u1 equation
			F[2*N+2*N2+j] = Ftilde3_nd[j,i];	# Forcing the h1 equation
			F[3*N+2*N2+j] = Ftilde6_nd[j,i];	# Forcing the h2 equation
		for j in range(0,N-2):
			F[2*N+j] = Ftilde2_nd[j+1,i];		# Forcing the v1 equation
		

		solution[:,i] = np.linalg.solve(A,F);

	u1tilde_nd = np.zeros((N,N),dtype=complex);
	u2tilde_nd = np.zeros((N,N),dtype=complex);
	v1tilde_nd = np.zeros((N,N),dtype=complex);
	v2tilde_nd = np.zeros((N,N),dtype=complex);
	eta0tilde_nd = np.zeros((N,N),dtype=complex);
	eta1tilde_nd = np.zeros((N,N),dtype=complex);
	for j in range(0,N):
		u1tilde_nd[j,:] = solution[j,:];
		u2tilde_nd[j,:] = solution[N+j,:];
		eta0tilde_nd[j,:] = solution[2*N+2*N2+j];
		eta1tilde_nd[j,:] = solution[3*N+2*N2+j];		
	for j in range(0,N2):
		v1tilde_nd[j+1,:] = solution[2*N+j,:];
		v2tilde_nd[j+1,:] = solution[2*N+N2+j,:];

	return u1tilde_nd, u2tilde_nd, v1tilde_nd, v2tilde_nd, eta0tilde_nd, eta1tilde_nd;

#=======================================================

# SPEC_TO_PHYS
#=======================================================
def SPEC_TO_PHYS(u1tilde_nd,u2tilde_nd,v1tilde_nd,v2tilde_nd,eta0tilde_nd,eta1tilde_nd,T_nd,Nt,dx_nd,omega_nd,N):
# Function takes the spectral-physical solutions produced by the solver and returns the time-dependent solutions in physical space.
	
	I = np.complex(0,1);	
	
	u1_nd = np.zeros((N,N,Nt),dtype=complex);
	u2_nd = np.zeros((N,N,Nt),dtype=complex);
	v1_nd = np.zeros((N,N,Nt),dtype=complex);
	v2_nd = np.zeros((N,N,Nt),dtype=complex);
	eta0_nd = np.zeros((N,N,Nt),dtype=complex);
	eta1_nd = np.zeros((N,N,Nt),dtype=complex);
	
	for ti in range(0,Nt):
		# Calculate the solutions in physical space at some instant t. Solutions are divided (later) by extra factor AmpF, so that they are normalised by the forcing amplitude.
		u1_nd[:,:,ti] = np.exp(2.*np.pi*I*omega_nd*T_nd[ti])*np.fft.ifft(u1tilde_nd,axis=1) / dx_nd;
		u2_nd[:,:,ti] = np.exp(2.*np.pi*I*omega_nd*T_nd[ti])*np.fft.ifft(u2tilde_nd,axis=1) / dx_nd; 
		v1_nd[:,:,ti] = np.exp(2.*np.pi*I*omega_nd*T_nd[ti])*np.fft.ifft(v1tilde_nd,axis=1) / dx_nd;
		v2_nd[:,:,ti] = np.exp(2.*np.pi*I*omega_nd*T_nd[ti])*np.fft.ifft(v2tilde_nd,axis=1) / dx_nd; 		
		eta0_nd[:,:,ti] = np.exp(2.*np.pi*I*omega_nd*T_nd[ti])*np.fft.ifft(eta0tilde_nd,axis=1) / dx_nd;
		eta1_nd[:,:,ti] = np.exp(2.*np.pi*I*omega_nd*T_nd[ti])*np.fft.ifft(eta1tilde_nd,axis=1) / dx_nd;

	return u1_nd, u2_nd, v1_nd, v2_nd, eta0_nd, eta1_nd;







