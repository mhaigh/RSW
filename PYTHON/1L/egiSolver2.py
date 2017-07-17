def EIG_COEFFICIENTS(Ro,Re,k,f_nd,U0_nd,H0_nd,gamma_nd,dy_nd,N):
# Here we also divide coeficients of derivatives by the relevant space-step: 
# 2*dy_nd for first-order derivatives, dy_nd**2 for second-order derivatives.

	I = np.complex(0.0,1.0);		# Define I = sqrt(-1)

	# Coefficients with no k- or y-dependence
	a2 = I / (2. * np.pi * Re * dy_nd**2);	# Note: b3=a2, so instead of defining b3, we just use a2 in its place, saving time.
	b4 = - I / (4. * np.pi * Ro * dy_nd);

	# Coefficients with k-dependence only
	a4 = k / Ro;

	# Coefficients with y-dependence only

	a3 = I * (f_nd / Ro - diff(U0_nd,2,0,dy_nd)) / (2. * np.pi);		# For uniform U0_nd, the first term is zero	

	b1 = - I * f_nd / (2. * np.pi * Ro);
	
	c2 = - I * diff(H0_nd,2,0,dy_nd) / (2. * np.pi);					# For zero BG flow, H0_nd=Hflat=const, i.e. c2=0
	
	c3 = - I * H0_nd / (4. * np.pi * dy_nd);

	# a3, c2 and c3 have length N, so when used in defining u and v equations in the matrix below,
	# we must add 1 to each index in order to skip out the dead gridpoint.

	
	# Coefficients dependent on k and y
	a1 = np.zeros(N,dtype=complex);		# Note: b2=a1
	c1 = np.zeros(N,dtype=complex);
	c4 = np.zeros(N,dtype=complex);
	for j in range(0,N):
		a1[j] = U0_nd[j] * k - 2. * np.pi * I * k**2 / Re - I * gamma_nd / (2. * np.pi * Ro);
		c1[j] = k * H0_nd[j];
		c4[j] = U0_nd[j] * k;	

	return a1,a2,a3,a4,b1,b4,c1,c2,c3,c4;
 	

 	# Coefficients dependent on k and y

 	a1 = np.zeros(N,dtype=complex);		# Note: b2=a1

 	c1 = np.zeros(N,dtype=complex);

 	c4 = np.zeros(N,dtype=complex);

 	for j in range(0,N):

 		a1[j] = U0_nd[j] * k - 2. * np.pi * I * k**2 / Re - I * gamma_nd / (2. * np.pi * Ro);

 		c1[j] = k * H0_nd[j];

 		c4[j] = U0_nd[j] * k;	

#====

	return a1,a2,a3,a4,b1,b4,c1,c2,c3,c4;

def NO_SLIP_EIG(a1,a2,a3,a4,b1,b4,c1,c2,c3,c4,Ro,N,N2):
# Called by EIG.py if BC = 'NO-SLIP'.
# Also impose that eta_y = 0 on the boundaries, allowing to write 1st-order derivatives of v at the boundaries in a one-sided FD approx.

	dim = 2 * N2 + N;
	#print(dim);

	A = np.zeros((dim,dim),dtype=complex);	
				
	# First the boundary terms. Some of these could be defined in the upcoming loop,
	# but for simplicity we define all boundary terms here.

	# u equation BCs
	# South
	A[0,0] = a1[1] - 2. * a2;		# u[1]
	A[0,1] = a2;					# u[2]
	A[0,N2] = a3[1];				# v[1]
	A[0,2*N2+1] = a4;				# eta[1] 
	# North
	A[N2-1,N2-1] = a1[N2] - 2. * a2;	# u[N-2]
	A[N2-1,N2-2] = a2;					# u[N-3]
	A[N2-1,2*N2-1] = a3[N2]				# v[N-2] 
	A[N2-1,2*N2+N-2] = a4;				# eta[N-2]

	# v equation BCs
	# South
	A[N2,0] = b1[1];					# u[1]
	A[N2,N2] = a1[1] - 2. * a2;			# v[1]
	A[N2,N2+1] = a2;					# v[2]
	A[N2,2*N2] = - b4;					# eta[0]	
	A[N2,2*N2+2] = b4;					# eta[2]
	# North
	A[2*N2-1,N2-1] = b1[N2];					# u[N-2] 
	A[2*N2-1,2*N2-1] = a1[N2] - 2. * a2;		# v[N-2]
	A[2*N2-1,2*N2-2] = a2;						# v[N-3]
	A[2*N2-1,2*N2+N-1] = b4;					# eta[N-1]
	A[2*N2-1,2*N2+N-3] = - b4;					# eta[N-3]

	# eta equation BCs (Here we have to define BCs at j=2*N2,2*N2+1,2*N2+N-1,2*N2+N-2)
	# South
	A[2*N2,N2] = 2. * c3[0];			# v[1] (factor of 2 because we use one-sided FD for v_y at the boundaries)
	A[2*N2,2*N2] = c4[0];				# eta[0]
	A[2*N2+1,0] = c1[1];				# u[1]
	A[2*N2+1,N2] = c2[1];				# v[1]
	A[2*N2+1,N2+1] = c3[1];				# v[2]
	A[2*N2+1,2*N2+1] = c4[1];			# eta[1]
	# North
	A[2*N2+N-1,2*N2-1] = - 2. * c3[N-1];	# v[N-2]
	A[2*N2+N-1,2*N2+N-1] = c4[N-1];			# eta[N-1]
	A[2*N2+N-2,N2-1] = c1[N2]; 				# u[N-2]
	A[2*N2+N-2,2*N2-1] = c2[N2];			# v[N-2]
	A[2*N2+N-2,2*N2-2] = - c3[N2];			# v[N-3]
	A[2*N2+N-2,2*N2+N-2] = c4[N2];			# eta[N-2]

	# Now the inner values - two loops required: one for u,v, and one for eta
	for j in range(1,N2-1):
		# u equation
		A[j,j] = a1[j+1] - 2. * a2;			# u[j+1]
		A[j,j-1] = a2;						# u[j]
		A[j,j+1] = a2;						# u[j+2]
		A[j,N2+j] = a3[j+1];				# v[j+1]
		A[j,2*N2+j+1] = a4;					# eta[j+1] 
			
	for j in range(1,N2-1):
		# v equation
		A[N2+j,j] = b1[j+1];					# u[j+1]
		A[N2+j,N2+j] = a1[j+1] - 2. * a2;		# v[j+1]
		A[N2+j,N2+j-1] = a2;					# v[j]
		A[N2+j,N2+j+1] = a2;					# v[j+2]
		A[N2+j,2*N2+j] = - b4;					# eta[j]
		A[N2+j,2*N2+j+2] = b4;					# eta[j+2]

	for j in range(2,N-2):
		# eta equation
		A[2*N2+j,j-1] = c1[j];			# u[j] 
		A[2*N2+j,N2+j-1] = c2[j];		# v[j]
		A[2*N2+j,N2+j] = c3[j];			# v[j+1]
		A[2*N2+j,N2+j-2] = - c3[j];		# v[j-1]
		A[2*N2+j,2*N2+j] = c4[j];		# eta[j]
	
	val,vec = np.linalg.eig(A);

	# Keep the set of eigenvalues as is; separate out the eigenvectors to be returned.
	u_vec = np.zeros((N,dim),dtype=complex);
	u_vec[1:N-1,:] = vec[0:N2,:];
	v_vec = np.zeros((N,dim),dtype=complex);
	v_vec[1:N-1,:] = vec[N:N+N2,:];
	eta_vec = vec[2*N2:2*N2+N,:];	

	return val, u_vec, v_vec, eta_vec;
