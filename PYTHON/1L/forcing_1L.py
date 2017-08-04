# forcing_1L
#=======================================================

import numpy as np
import matplotlib.pyplot as plt
from diagnostics import diff, extend
# Forcing 
#=======================================================

# This function is used to define the three forcing terms on the 1-L SW equations.
# F1 is the forcing in physical space applied to the zonal momentum equation.
# F2 is the forcing in physical space applied to the meridional momentum equation.
# F3 is the forcing in physical space applied to the continuity equation.
# These forcings should be defined in normalised form (i.e. sum(Fi=1)), so that they are
# dimensionless and their amplitudes and dimensions are stored in the alphai coefficients.
# For geostrophically balanced forcing, the alphai are functionally related.

#=======================================================

def Forcing(x,y,K,y0,r0,N,FORCE,AmpF,g,f,f0,U,L,dx,dy):
	
	I = np.complex(0,1);
	F1 = np.zeros((N,N+1));
	F2 = np.zeros((N,N+1));
	F3 = np.zeros((N,N+1));
	if FORCE == 'BALANCED':
		count = 0;
		mass = 0;
		for i in range(0,N+1):
			for j in range(0,N):
				r = np.sqrt(x[i]**2 + (y[j]-y0)**2);
				if r < r0:
					count = count + 1;
					if r == 0:
						F1[j,i] = 0;
						F2[j,i] = 0;						
					else:	
						F1[j,i] = AmpF * np.pi * g * (y[j]-y0) / (2 * r0 * f[j] * r) * np.sin((np.pi / 2) * r / r0);
						F2[j,i] = - AmpF * np.pi * g * x[i] / (2 * r0 * f[j] * r) * np.sin((np.pi / 2) * r / r0);
					F3[j,i] = AmpF * np.cos((np.pi / 2) * r / r0);
					mass = mass + F3[j,i];
		mass = mass / (N*(N+1) - count);
		for i in range(0,N+1):
			for j in range(0,N):
				F3[j,i] = F3[j,i] - mass;
		#F3x = diff(F3,1,1,dx);
		#F3y = diff(F3,0,0,dy);
		#for j in range(0,N):
		#	F1[j,:] = - g * F3y[j,:] / f[j];
		#	F2[j,:] = g * F3x[j,:] / f[j];
			

	if FORCE == 'BUOYANCY':
		count = 0;
		mass = 0;
		for i in range(0,N):
			for j in range(0,N):
				r = np.sqrt(x[i]**2 + (y[j]-y0)**2);
				if r<r0:
					count = count + 1;
					F3[j,i] = AmpF * np.cos((np.pi / 2) * r / r0);
					mass = mass + F3[j,i];
		mass = mass / (N**2 - count);
		for i in range(0,N):
			for j in range(0,N):
				r = np.sqrt(x[i]**2 + (y[j]-y0)**2);
				if r >= r0:
					F3[j,i] = - mass;

	PLOT = False;
	if PLOT:
		aa = 1./24;
		Fmax = np.max(np.max(F3,axis=0));
		Flim = np.max(abs(F3/(1.05*Fmax)));
		#F3[0,0] = - Fmax;
		x_nd = x / L;
		y_nd = y / L;
		plt.contourf(x_nd,y_nd,F3/(1.1*Fmax));
		plt.xlabel('x',fontsize=22);
		plt.ylabel('y',fontsize=22);
		plt.text(y_nd[N-130],x_nd[N-130],'F3',color='k',fontsize=22);
		plt.arrow(-aa,2*aa+0.25,2*aa,0,head_width=0.7e5/L, head_length=0.7e5/L,color='k');
		plt.arrow(2*aa,aa+.25,0,-2*aa,head_width=0.7e5/L, head_length=0.7e5/L,color='k');
		plt.arrow(aa,-2*aa+.25,-2*aa,0,head_width=0.7e5/L, head_length=0.7e5/L,color='k');
		plt.arrow(-2*aa,-aa+.25,0,2*aa,head_width=0.7e5/L, head_length=0.7e5/L,color='k');
		plt.xticks((-1./2,0,1./2),['-0.5','0','0.5'],fontsize=14);
		plt.yticks((-1./2,0,1./2),['-0.5','0','0.5'],fontsize=14);
		plt.clim(-Flim,Flim);
		plt.colorbar();
		plt.tight_layout();
		plt.show();

		

	if FORCE == 'VORTICITY':
		for i in range(0,N):
			for j in range(0,N):
				r = np.sqrt(x[i]**2 + (y[j]-y0)**2);
				if r<r0:
					F1[j,i] = AmpF * np.pi * g * (y[j]-y0) / (2 * r0 * f[j] * r) * np.sin((np.pi / 2) * r / r0);
					F2[j,i] = - AmpF * np.pi * g * x[i] / (2 * r0 * f[j] * r) * np.sin((np.pi / 2) * r / r0);
	
	
	# Lastly, Fourier transform the three forcings in the x-direction
		
	Ftilde1 = dx * np.fft.hfft(F1,N,axis=1);	# Multiply by dx_nd as FFT differs by this factor compared to FT.
	Ftilde3 = dx * np.fft.hfft(F3,N,axis=1); 
	Ftilde2 = np.zeros((N,N),dtype=complex);
	for j in range(0,N):
		for i in range(0,N):
			Ftilde2[j,i] = 2 * np.pi * g * I * K[i] * Ftilde3[j,i] / f[j];

	# Nondimensionalise forcing terms
	#=======================================================

	F1_nd = F1 / (f0 * U);
	F2_nd = F2 / (f0 * U);
	F3_nd = F3 * g / (f0 * U**2); 

	Ftilde1_nd = Ftilde1 / (f0 * U * L);
	Ftilde2_nd = Ftilde2 / (f0 * U * L);
	Ftilde3_nd = Ftilde3 * g / (f0 * U**2 * L);

	return F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd;

#=======================================================

# forcingInv
def forcingInv(Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,x_nd,y_nd,dx_nd,N):
# A function that calculates the inverse of the forcing the check that the original forcing is found.

	F1i = np.real(np.fft.ihfft(Ftilde1_nd,axis=1)) / dx_nd;
	F2 = np.fft.ifft(Ftilde2_nd,axis=1) / dx_nd;
	F3i = np.real(np.fft.ihfft(Ftilde3_nd,axis=1)) / dx_nd;

	F1 = np.zeros((N,N));
	F3 = np.zeros((N,N));
	for i in range(0,N/2):
		F1[:,i] = F1i[:,i];
		F1[:,N-1-i] = F1i[:,i]; 
		F3[:,i] = F3i[:,i];
		F3[:,N-1-i] = F3i[:,i];
	
	F1 = extend(F1);
	F2 = extend(F2);
	F3 = extend(F3);

	plt.figure(2);

	plt.subplot(331);
	plt.contourf(x_nd,y_nd,F1);
	plt.colorbar();
	plt.subplot(332);
	plt.contourf(x_nd,y_nd,F2);
	plt.colorbar();
	plt.subplot(333);
	plt.contourf(x_nd,y_nd,F3);
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

#=======================================================

# forcingDiff
def forcingDiff(Ftilde_nd,y_nd,dy_nd,N,i):
# Plots the y-derivatives of the 1D forcing at a given wavenumber.
# A function used to test the diff algortithm, improve it, and resolve the error issue.

	Ftilde_y = diff(Ftilde_nd[:,i],2,0,dy_nd);
	Ftilde_yy = diff(Ftilde_y,2,0,dy_nd);

	plt.subplot(131);
	plt.plot(Ftilde_nd[:,i],y_nd);
	plt.subplot(132);
	plt.plot(Ftilde_y,y_nd);
	plt.subplot(133);
	plt.plot(Ftilde_yy,y_nd);
	plt.show();

	import sys
	sys.exit();







