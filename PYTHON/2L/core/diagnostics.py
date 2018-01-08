# diagnostics.py
# File containing functions to be called by the master script RSW_visc_1L.py.
# Diagnostic functions include plotting tools, differentiation operators and error calculators.
#====================================================

import numpy as np
import matplotlib.pyplot as plt

from math import fmod

#====================================================

# solutionPlots
# Solution plots
def solutionPlots(x_nd,y_nd,u1_nd,v1_nd,eta0_nd,u2_nd,v2_nd,eta1_nd,ts,FORCE1,BG1,Fpos,N):

	u1_nd = extend(u1_nd);
	u2_nd = extend(u2_nd);
	v1_nd = extend(v1_nd);
	v2_nd = extend(v2_nd);
	eta0_nd = extend(eta0_nd);
	eta1_nd = extend(eta1_nd);
	
	plt.figure(1,figsize=[15,7]);
	plt.subplot(231)
	plt.contourf(x_nd,y_nd,u1_nd[:,:,ts])
	plt.text(0.2,0.4,'u1',fontsize=18)
	plt.colorbar()
	plt.subplot(232)
	plt.contourf(x_nd,y_nd,v1_nd[:,:,ts])
	plt.text(0.2,0.4,'v1',fontsize=18)
	plt.colorbar()
	plt.subplot(233)
	plt.contourf(x_nd,y_nd,eta0_nd[:,:,ts])
	plt.text(0.2,0.4,'eta0',fontsize=18)
	plt.colorbar()
	plt.subplot(234)
	plt.contourf(x_nd,y_nd,u2_nd[:,:,ts])
	plt.text(0.2,0.4,'u2',fontsize=18)
	plt.colorbar()
	plt.subplot(235)
	plt.contourf(x_nd,y_nd,v2_nd[:,:,ts])
	plt.text(0.2,0.4,'v2',fontsize=18)
	plt.colorbar()
	plt.subplot(236)
	plt.contourf(x_nd,y_nd,eta1_nd[:,:,ts])
	plt.text(0.2,0.4,'eta1',fontsize=18)
	plt.colorbar()
	plt.tight_layout()
	plt.show()
	

	#plt.savefig('/home/mike/Documents/GulfStream/Code/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/' + str(Fpos) + '_'  + str(N) + '.png');
	plt.show();

#====================================================

# isolutionPlots
# Plots of phase and amplitude 
def isolutionPlots(x_nd,y_nd,u_nd,v_nd,eta_nd,PV_prime,ts,FORCE,BG,Fpos,N):

	u_nd = extend(u_nd);
	v_nd = extend(v_nd);
	eta_nd = extend(eta_nd);
	
	plt.figure(3);

	plt.subplot(231);
	plt.contourf(np.angle(u_nd[:,:,ts]));
	plt.title('u');
	plt.ylabel('PHASE');
	plt.subplot(234);
	plt.contourf(np.absolute(u_nd[:,:,ts]));
	plt.ylabel('AMP');
	plt.subplot(232);
	plt.contourf(np.angle(v_nd[:,:,ts]));
	plt.title('v');
	plt.subplot(235);
	plt.contourf(np.absolute(v_nd[:,:,ts]));
	plt.subplot(233);
	plt.contourf(np.angle(eta_nd[:,:,ts]));
	plt.title('eta');
	plt.subplot(236);
	plt.contourf(np.absolute(eta_nd[:,:,ts]));
	plt.tight_layout()

	plt.show()

#====================================================

# bgPlots
# Background state plots
def bgPlots(y_nd,H0_nd,U0_nd,PV_BG):

	plt.figure(2);
	plt.subplot(131);
	plt.plot(H0_nd,y_nd);
	plt.yticks((-1./2,0,1./2));
	plt.title('BG SSH');
	plt.subplot(132);
	plt.plot(U0_nd,y_nd);
	plt.yticks((-1./2,0,1./2));
	plt.title('BG flow U0');
	plt.subplot(133);
	plt.plot(PV_BG,y_nd);
	plt.yticks((-1./2,0,1./2));
	plt.title('BG PV');
	plt.show()

#====================================================

# forcingPlots
# Forcing plots 
def forcingPlots(x_nd,y_nd,F1_nd,F2_nd,F3_nd,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N):
# Function that plots the forcing, and its Fourier representation.

	plt.figure(1);

	plt.subplot(331);
	plt.contourf(x_nd,y_nd,F1_nd);
	plt.xticks((-1./2,0,1./2));
	plt.yticks((-1./2,0,1./2));
	plt.colorbar();
	plt.subplot(332);
	plt.contourf(x_nd,y_nd,F2_nd);
	plt.xticks((-1./2,0,1./2));
	plt.yticks((-1./2,0,1./2));
	plt.colorbar();
	plt.subplot(333);
	plt.contourf(x_nd,y_nd,F3_nd);
	plt.xticks((-1./2,0,1./2));
	plt.yticks((-1./2,0,1./2));
	plt.colorbar()

	plt.subplot(334);
	plt.contourf(np.real(Ftilde1_nd));
	plt.colorbar()
	plt.subplot(335);
	plt.contourf(np.real(Ftilde2_nd));
	plt.colorbar()
	plt.subplot(336);
	plt.contourf(np.real(Ftilde3_nd));
	plt.colorbar()

	plt.subplot(337);
	plt.contourf(np.imag(Ftilde1_nd));
	plt.colorbar()
	plt.subplot(338);
	plt.contourf(np.imag(Ftilde2_nd));
	plt.colorbar()
	plt.subplot(339);
	plt.contourf(np.imag(Ftilde3_nd));
	plt.colorbar()

	plt.show();

#====================================================

# diff
# Function for differentiating a vector
def diff(f,d,p,delta):

# f[y,x] is the function to be differentiated.
# d is the direction in which the differentiation is to be taken:
# d=0 for differentiation over the first index, d=1 for the second.
# d=2 for a 1D vector
# p is a periodic switch:
# p=1 calculates periodic derivative.
	
	if d != 2:
		dimx = np.shape(f)[1];		# Finds the number of gridpoints in the x and y directions
		dimy = np.shape(f)[0];
		df = np.zeros((dimy,dimx),dtype=f.dtype);
	else:
		dimy = np.shape(f)[0];
		df = np.zeros(dimy,dtype=f.dtype);	
	
	if p == 0:
		# Solid boundary derivative.
		# Note multiplication by 2 of boundary terms are to
		# invert the division by 2 at the end of the module.
		if d == 0:
		
			df[1:dimy-1,:] = f[2:dimy,:] - f[0:dimy-2,:];
		
			df[0,:] = 2 * (f[1,:] - f[0,:]);
			df[dimy-1,:] = 2 * (f[dimy-1,:] - f[dimy-2,:]);
	
		elif d == 1:

			df[:,1:dimx-1] = f[:,2:dimx] - f[:,0:dimx-2];

			df[:,0] = 2 * (f[:,1] - f[:,0]);
			df[:,dimx-1] = 2 * (f[:,0] - f[:,dimx-2]);
		
		elif d == 2:

			df[1:dimy-1] = f[2:dimy] - f[0:dimy-2];

			df[0] = 2 * (f[1] - f[0]);
			df[dimy-1] = 2 * (f[dimy-1] - f[dimy-2]);

		else:
			print('error')

	elif p == 1:
		# Periodic option

		if d == 0:

			df[1:dimy-1,:] = f[2:dimy,:] - f[0:dimy-2,:];	

			df[0,:] = f[1,:] - f[dimy-1,:];
			df[dimy-1,:] = f[0,:] - f[dimy-2,:];
	
		elif d == 1:

			df[:,1:dimx-1] = f[:,2:dimx]-f[:,0:dimx-2];

			df[:,0] = f[:,1] - f[:,dimx-1];
			df[:,dimx-1] = f[:,0] - f[:,dimx-2];

		elif d == 2:

			df[1:dimy-1]=f[2:dimy] - f[0:dimy-2];

			df[0] = f[1] - f[0];
			df[dimy-1] = f[0] - f[dimy-2];
		
		else:
			print('error')

	else:
		print('error')

	df = 0.5 * df / delta;

	return df

#====================================================

# ddt
# A time-derivative function.
def ddt(f,delta):
# Only takes two inputs, the function itself and the size of the timestep as defined in input_File_1L.
# f should be a 3-D vector, where the 3rd index is the time index.
	dimx, dimy, dimt = np.shape(f);
	df = np.zeros((dimy,dimx,dimt),dtype=f.dtype);
	
	df[:,:,1:dimt-1] = f[:,:,2:dimt] - f[:,:,0:dimt-2];	# centered fd

	df[:,:,0] = f[:,:,1] - f[:,:,dimt-1];				# The two boundary terms
	df[:,:,dimt-1] = f[:,:,0] - f[:,:,dimt-2];

	df = 0.5 * df / delta;
	
	return df

#====================================================

# error
def error(u_nd,v_nd,eta_nd,dx_nd,dy_nd,dt_nd,U0_nd,H0_nd,Ro,gamma_nd,Re,f_nd,F1_nd,F2_nd,F3_nd,T_nd,ts,omega_nd,N):
# This function calculates the error of the 1L SW solutions, and is to be used in the main code RSW_visc_1L.py.
# It takes the half-spectral, half-physical solutions (u_tilde etc) and performs the inverse Fourier transform, but retains both the real and imaginary parts of the solution.
# Both parts of the solution are required as the solutions' time-dependence are complex.

	I = np.complex(0,1);
	
	# Now we calculate all the relevant x and y derivatives
	u_y = diff(u_nd[:,:,ts],0,0,dy_nd);
	u_yy = diff(u_y[:,:],0,0,dy_nd);
	u_x = diff(u_nd[:,:,ts],1,1,dx_nd);
	u_xx = diff(u_x[:,:],1,1,dx_nd);

	v_y = diff(v_nd[:,:,ts],0,0,dy_nd);
	v_yy = diff(v_y[:,:],0,0,dy_nd);
	v_x = diff(v_nd[:,:,ts],1,1,dx_nd);
	v_xx = diff(v_x[:,:],1,1,dx_nd);

	eta_x = diff(eta_nd[:,:,ts],1,1,dx_nd);
	eta_y = diff(eta_nd[:,:,ts],0,0,dy_nd);

	U0_y = diff(U0_nd,2,0,dy_nd);
	H0_y = diff(H0_nd,2,0,dy_nd);

	# t derivatives
	u_t = ddt(u_nd,dt_nd);
	v_t = ddt(v_nd,dt_nd);
	eta_t = ddt(eta_nd,dt_nd);

	error1 = np.zeros((N,N),dtype=complex);
	error2 = np.zeros((N,N),dtype=complex);
	error3 = np.zeros((N,N),dtype=complex);
	for j in range(0,N):
		for i in range(0,N):
			error1[j,i] = Ro * (u_t[j,i,ts] + U0_nd[j] * u_x[j,i]) + gamma_nd * u_nd[j,i,ts] - Ro * (u_xx[j,i] + u_yy[j,i]) / Re + (Ro * U0_y[j] - f_nd[j]) * v_nd[j,i,ts] + eta_x[j,i] - 1 * F1_nd[j,i] * np.exp(I * omega_nd * T_nd[ts]);
			error2[j,i] = Ro * (v_t[j,i,ts] + U0_nd[j] * v_x[j,i]) + gamma_nd * v_nd[j,i,ts] - Ro * (v_xx[j,i] + v_yy[j,i]) / Re + f_nd[j] * u_nd[j,i,ts] + eta_y[j,i] - F2_nd[j,i] * np.exp(I * omega_nd * T_nd[ts]);
			error3[j,i] = eta_t[j,i,ts] + U0_nd[j] * eta_x[j,i] + H0_nd[j] * u_x[j,i] + H0_nd[j] * v_y[j,i] + H0_y[j] * v_nd[j,i,ts] - F3_nd[j,i] * np.exp(I * omega_nd * T_nd[ts]);
	
	plt.subplot(131)
	plt.contourf(error1)
	plt.colorbar()
	plt.subplot(132)
	plt.contourf(error2)
	plt.colorbar()
	plt.subplot(133)
	plt.contourf(error3)
	plt.colorbar()
	plt.show()

	error1 = np.real(error1);

	zero = np.zeros((N,N+1));
	e1 = np.sqrt((error1**2).mean())
	e2 = np.sqrt((error2**2).mean())
	e3 = np.sqrt((error3**2).mean())

	return e1, e2, e3

#====================================================

# extend
# A function used to replace the extra x-gridpoint on a solution.
def extend(f):

	dimx = np.shape(f)[1];
	dimy = np.shape(f)[0];
	if f.size != dimx * dimy:
		dimt = np.shape(f)[2];

		f_new = np.zeros((dimy,dimx+1,dimt),dtype=f.dtype);
		for i in range(0,dimx):
			f_new[:,i,:] = f[:,i,:];
	
		f_new[:,dimx,:] = f[:,0,:];
	
	else:
		f_new = np.zeros((dimy,dimx+1),dtype=f.dtype);
		for i in range(0,dimx):
			f_new[:,i] = f[:,i];
	
		f_new[:,dimx] = f[:,0];

	return f_new	 


#====================================================

# permute
# A function that permutes a solution in x, to 'move' the location of the forcing.
def permute(u,N,Nt,x0):

	u_perm = np.zeros((N,N,Nt));
	xi = - N / 4;
	for i in range(0,N):
		ii = fmod(i+xi,N);
		u_perm[:,ii,:] = u[:,i,:]
	
	return u_perm




