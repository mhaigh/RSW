# plotEEFs
#========================================

import numpy as np
import matplotlib.pyplot as plt

from inputFile_1L import *

#========================================

# Options are: 1 for y0 test; 2 for U0 test.  
option = 2;

if option == 1:
	# Must make sure input file has BG = GAUSSIAN/UNIFORM appropriately. 
	# Comment out U0 plotting terms if necessary.
	E = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/EEFs/EEF_y0_GAUSSIAN.npy');
	E50 = E[:,0];
	E60 = E[:,1];
	E70 = E[:,2];

	NU = np.shape(E50)[0];
	N_unforced = int(N - NU); 					# The number of y gridpoints for which the test was not run
	
	y_forced = np.zeros(NU);
	U0_forced = np.zeros(NU);
	Q_forced = np.zeros(NU);
	for j in range(0,NU):
		y_forced[j] = y_nd[j+N_unforced/2];
		U0_forced[j] = U0_nd[j+N_unforced/2];
		Q_forced[j] = Q[j+N_unforced/2];
	ylin = np.zeros(NU);
	for j in range(0,NU/2):
		ylin[j] = 0.06 + - 5 * y_forced[j];
		ylin[j+NU/2] = 1 + 0.06 + 2 * y_forced[j];
 
	plt.plot(y_forced,E50,'b-',label='50 days',linewidth=1.3);
	plt.plot(y_forced,E60,'r-',label='60 days',linewidth=1.3);
	plt.plot(y_forced,E70,'g-',label='70 days',linewidth=1.3);
	#plt.plot(y_forced,ylin);
	plt.plot(y_forced,2.*U0_forced,'k--',label='BG FLOW',linewidth=1.3);
	#plt.plot(1000000*Q_forced,y_forced,'k--',label='BG',linewidth=2);
	plt.xticks((-1./2,0,1./2),['-1/2','0','1/2']);
	plt.xlim(-0.5,0.5);
	#plt.ylim(-0.2,1.6);
	#plt.title('Equivalent Eddy Fluxes - Gaussian BG flow',fontsize=18);
	#plt.title('Equivalent Eddy Fluxes - Uniform BG flow = 15 cm/s',fontsize=18);
	plt.xlabel('FORCING LATITUDE',fontsize=22);
	#plt.ylabel('EEF',fontsize=18);
	plt.legend(prop={'size':20});
	plt.tight_layout();
	plt.show();
	
	a = 40;
	plt.plot(E50[NU/2-a:NU/2+a],y_forced[NU/2-a:NU/2+a],'b-',label='50 days',linewidth=1);
	plt.plot(E60[NU/2-a:NU/2+a],y_forced[NU/2-a:NU/2+a],'r-',label='60 days',linewidth=1);
	plt.plot(E70[NU/2-a:NU/2+a],y_forced[NU/2-a:NU/2+a],'g-',label='70 days',linewidth=1);
	#plt.title('Equivalent Eddy FLuxes',fontsize=18);
	plt.tight_layout();
	plt.show();
	
	
if option == 2:

	E = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/EEFs/EEF_U0_NORTH.npy');
	E50 = E[:,0];
	E60 = E[:,1];
	E70 = E[:,2];

	NU = np.shape(E50)[0];
	U0 = np.linspace(-0.3,0.5,NU);

	plt.figure(1);
	plt.plot(U0,E50,'b-',label='50 days',linewidth=1.3);
	plt.plot(U0,E60,'r-',label='60 days',linewidth=1.3);
	plt.plot(U0,E70,'g-',label='70 days',linewidth=1.3);
	plt.axhline(0,color='k',ls='--');
	plt.axvline(0,color='k',ls='--');
	plt.xlim(U0[0],U0[NU-1]);
	plt.ylim(-0.1,1.7);
	plt.ylabel('EEF',fontsize=22);
	plt.xlabel('U0',fontsize=22);
	plt.legend(prop={'size':20});
	plt.tight_layout();
	plt.show();




