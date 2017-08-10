# Other
#========================================

import numpy as np
import matplotlib.pyplot as plt

import diagnostics

from inputFile_1L import *

#========================================
option = 4;

if option == 0:
	cos = np.cos(2*np.pi*x_nd/(x_nd[N]-x_nd[0]));
	cos_y = diagnostics.diff(cos,2,1,dy_nd) / (2*np.pi);
	cos_yy = diagnostics.diff(cos_y,2,1,dy_nd)  / (2*np.pi);

	dim = len(cos);
	print(cos[0],cos[dim-1]);
	plt.subplot(131);
	plt.plot(cos);
	plt.grid();
	plt.subplot(132);
	plt.plot(cos_y);
	plt.grid();
	plt.subplot(133);
	plt.plot(cos_yy);
	plt.grid();
	plt.show();

# Plots a selection of Gaussian profiles
if option == 1:
	if BG == 'UNIFORM':
		a = 2;

	if BG == 'GAUSSIAN':

		U01 = np.zeros(N);
		U02 = np.zeros(N);
		U03 = np.zeros(N);
		U04 = np.zeros(N);
		U05 = np.zeros(N);
		U06 = np.zeros(N);

		H01 = np.zeros(N);
		H02 = np.zeros(N);
		H03 = np.zeros(N);
		H04 = np.zeros(N);
		H05 = np.zeros(N);
		H06 = np.zeros(N);

		# U01 - reference Gaussian jet BG flow
		Umag = 0.2;
		sigma = 0.3 * Ly;			# Increasing sigma decreases the sharpness of the jet
		l = Ly / 2;
		a = Umag / (np.exp(l**2 / sigma**2) - 1);	# Maximum BG flow velocity Umag
		for j in range(0,N):
			U01[j] = a * np.exp((l**2 - y[j]**2) / sigma**2) - a;		# -a ensures U0 is zero on the boundaries
			H01[j] = - a * (np.sqrt(np.pi) * f0 * sigma * np.exp(l**2 / sigma**2) * erf(y[j] / sigma) / 2 
					- beta * sigma**2 * np.exp((l**2 - y[j]**2) / sigma**2) / 2
					- f0 * y[j] - beta * y[j]**2 / 2) / g + Hflat; #erf(0);
	
		# U02 - 1.5 * strenth BG flow
		Umag = 0.3;
		sigma = 0.3 * Ly;			
		l = Ly / 2;
		a = Umag / (np.exp(l**2 / sigma**2) - 1);	
		for j in range(0,N):
			U02[j] = a * np.exp((l**2 - y[j]**2) / sigma**2) - a;		
			H02[j] = - a * (np.sqrt(np.pi) * f0 * sigma * np.exp(l**2 / sigma**2) * erf(y[j] / sigma) / 2 
					- beta * sigma**2 * np.exp((l**2 - y[j]**2) / sigma**2) / 2
					- f0 * y[j] - beta * y[j]**2 / 2) / g + Hflat;
	
		# U03 - half strength BG flow
		Umag = 0.1;
		sigma = 0.3 * Ly;			
		l = Ly / 2;
		a = Umag / (np.exp(l**2 / sigma**2) - 1);	
		for j in range(0,N):
			U03[j] = a * np.exp((l**2 - y[j]**2) / sigma**2) - a;		
			H03[j] = - a * (np.sqrt(np.pi) * f0 * sigma * np.exp(l**2 / sigma**2) * erf(y[j] / sigma) / 2 
					- beta * sigma**2 * np.exp((l**2 - y[j]**2) / sigma**2) / 2
					- f0 * y[j] - beta * y[j]**2 / 2) / g + Hflat;
	
	
		# U04 - 'sharp' jet
		Umag = 0.2;
		sigma = 0.2 * Ly;			
		l = Ly / 2;
		a = Umag / (np.exp(l**2 / sigma**2) - 1);	
		for j in range(0,N):
			U04[j] = a * np.exp((l**2 - y[j]**2) / sigma**2) - a;		
			H04[j] = - a * (np.sqrt(np.pi) * f0 * sigma * np.exp(l**2 / sigma**2) * erf(y[j] / sigma) / 2 
					- beta * sigma**2 * np.exp((l**2 - y[j]**2) / sigma**2) / 2
					- f0 * y[j] - beta * y[j]**2 / 2) / g + Hflat;

		# U04 - 'sharper' jet
		Umag = 0.15;
		sigma = 0.2 * Ly;			
		l = Ly / 2;
		a = Umag / (np.exp(l**2 / sigma**2) - 1);	
		for j in range(0,N):
			U05[j] = a * np.exp((l**2 - y[j]**2) / sigma**2) - a;		
			H05[j] = - a * (np.sqrt(np.pi) * f0 * sigma * np.exp(l**2 / sigma**2) * erf(y[j] / sigma) / 2 
					- beta * sigma**2 * np.exp((l**2 - y[j]**2) / sigma**2) / 2
					- f0 * y[j] - beta * y[j]**2 / 2) / g + Hflat;
	
		# U05 - 'wide' jet
		Umag = 0.2;
		sigma = 0.4 * Ly;			
		l = Ly / 2;
		a = Umag / (np.exp(l**2 / sigma**2) - 1);	
		for j in range(0,N):
			U06[j] = a * np.exp((l**2 - y[j]**2) / sigma**2) - a;		
			H06[j] = - a * (np.sqrt(np.pi) * f0 * sigma * np.exp(l**2 / sigma**2) * erf(y[j] / sigma) / 2 
					- beta * sigma**2 * np.exp((l**2 - y[j]**2) / sigma**2) / 2
					- f0 * y[j] - beta * y[j]**2 / 2) / g + Hflat;
	
		U01 = U01 / U;
		U02 = U02 / U;
		U03 = U03 / U;
		U04 = U04 / U;
		U05 = U05 / U;
		U06 = U06 / U;
	
		H01 = H01 / H;
		H02 = H02 / H;
		H03 = H03 / H;
		H04 = H04 / H;
		H05 = H05 / H;
		H06 = H06 / H;
	
		plt.figure(1);
		plt.subplot(121);
		plt.plot(U01,y_nd,'b-',label='REFRENCE',linewidth=2);
		plt.plot(U02,y_nd,'r-',label='STRONG',linewidth=2);
		plt.plot(U03,y_nd,'g-',label='WEAK',linewidth=2);
		plt.plot(U04,y_nd,'k-',label='SHARP',linewidth=2);
		plt.plot(U05,y_nd,'m-',label='SHARPER',linewidth=2);
		plt.plot(U06,y_nd,'c-',label='WIDE',linewidth=2);
		plt.title('BG FLOW U0');
		plt.yticks((-1./2,0,1./2));
		plt.legend();
		plt.subplot(122);
		plt.plot(H01,y_nd,'b-',label='REFRENCE',linewidth=2);
		plt.plot(H02,y_nd,'r-',label='STRONG',linewidth=2);
		plt.plot(H03,y_nd,'g-',label='WEAK',linewidth=2);
		plt.plot(H04,y_nd,'k-',label='SHARP',linewidth=2);
		plt.plot(H05,y_nd,'m-',label='SHARPER',linewidth=2);
		plt.plot(H06,y_nd,'c-',label='WIDE',linewidth=2);
		plt.title('BG SSH H0')
		plt.yticks((-1./2,0,1./2));
		plt.legend();
		plt.tight_layout();
		plt.show();

# Plots the Rossby def. rad. against y
elif option == 2:
	Ld = np.zeros(N);
	for j in range(0,N):
		Ld[j] = np.sqrt(g * r0) / f[j];
	plt.plot(Ld,y)
	plt.show()

elif option == 3:
	E50 = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/EEFs/EEF_om50_y0.npy');
	E60 = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/EEFs/EEF_om60_y0.npy');
	E70 = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/EEFs/EEF_om70_y0.npy');
	NU = np.shape(E50)[0];
	#U0 = np.linspace(-0.3,0.5,NU);
	plt.plot(E50,'b-',label='50 days',linewidth=2);
	plt.plot(E60,'r-',label='60 days',linewidth=2);
	plt.plot(E70,'g-',label='70 days',linewidth=2);
	#plt.axhline(0,color='k',ls='--');
	#plt.axvline(0,color='k',ls='--');
	#plt.xlim(U0[0],U0[NU-1]);
	plt.title('Equivalent Eddy FLuxes',fontsize=18);
	plt.ylabel('EEF',fontsize=18);
	#plt.xlabel('U0',fontsize=18);
	plt.legend();
	plt.tight_layout();
	plt.show()

elif option == 4:
	
	N_set = [32,64,128,256,512,1024];
	p_set = [5,6,7,8,9,10];

	error1 = [1.16981294267e-19,1.45815391948e-19,1.46022581344e-19,1.10596154905e-19,8.14040437205e-20];
	error2 = [6.70739732654e-20,7.4529325698e-20,5.18529854606e-20,4.16850431505e-20,4.2138040004e-20];
	error3 = [1.16981294267e-19,1.45815391948e-19,1.46022581344e-19,1.10596154905e-19,8.14040437205e-20];

	plt.plot(N_set,error1);
	plt.show();




