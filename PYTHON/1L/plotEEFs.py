# plotEEFs
#========================================

import numpy as np
import matplotlib.pyplot as plt

from inputFile_1L import *

#========================================

# Options are: 1 for y0 test; 2 for U0 test.  
option = 4;

if option == 1:
	# Must make sure input file has BG = GAUSSIAN/UNIFORM appropriately. 
	# Comment out U0 plotting terms if necessary.
	E = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/EEF_y0_GAUSSIAN.npy');
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

	E = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/EEF_U0_NORTH.npy');
	E50 = E[:,0];
	E60 = E[:,1];
	E70 = E[:,2];

	NU = np.shape(E50)[0];
	U0 = np.linspace(-0.3,0.5,NU);
	print(NU);
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

if option == 3:

	E50 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/EEF_50_512.npy');
	E60 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/EEF_60_512.npy');
	E70 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/EEF_70_512.npy');

	E50 = E50[:,0,0] - E50[:,0,1];
	E60 = E60[:,0,0] - E60[:,0,1];
	E70 = E70[:,0,0] - E70[:,0,1];

	NU = np.shape(E50)[0];

	plt.figure(1);
	plt.plot(E50,'b-',label='50 days',linewidth=1.3);
	plt.plot(E60,'r-',label='60 days',linewidth=1.3);
	plt.plot(E70,'g-',label='70 days',linewidth=1.3);
	plt.show();

if option == 4:

	E50 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/EEF_50_512.npy');
	E60 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/EEF_60_512.npy');
	E70 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/EEF_70_512.npy');

	E50 = E50[:,0,0] - E50[:,0,1];
	E60 = E60[:,0,0] - E60[:,0,1];
	E70 = E70[:,0,0] - E70[:,0,1];

	NU = np.shape(E50)[0];
	
	excl = 20;

	#plt.figure(1);
	#plt.plot(E50,'b-',label='50 days',linewidth=1.3);
	#plt.plot(E60,'r-',label='60 days',linewidth=1.3);
	#plt.plot(E70,'g-',label='70 days',linewidth=1.3);
	#plt.show();
	
	E50 = E50[excl:NU-excl];
	E60 = E60[excl:NU-excl];
	E70 = E70[excl:NU-excl];

	NU = np.shape(E50)[0];
	N2 = NU/2;
	x = np.linspace(1,N2,N2);
	
	E50_north = E50[NU-N2:NU];
	E50_south = E50[0:N2];
	E50_south = E50_south[::-1];

	E60_north = E60[NU-N2:NU];
	E60_south = E60[0:N2];
	E60_south = E60_south[::-1];

	E70_north = E70[NU-N2:NU];
	E70_south = E70[0:N2];
	E70_south = E70_south[::-1];

	plt.plot(E50_north,label='NORTH');
	plt.plot(E50_south,label='SOUTH');	
	plt.legend();
	plt.show();

	inv50 = E50_south/E50_north;
	inv60 = E60_south/E60_north;
	inv70 = E70_south/E70_north;
	
	d = 2;
	pc50 = np.polyfit(x,inv50,d);
	pc60 = np.polyfit(x,inv60,d);
	pc70 = np.polyfit(x,inv70,d);
	print(pc50,pc60,pc70);

	p50 = np.zeros(N2);
	p60 = np.zeros(N2);
	p70 = np.zeros(N2);
	for di in range(0,d+1):
		p50 = p50 + pc50[di] * x**(d-di);
		p60 = p60 + pc60[di] * x**(d-di);
		p70 = p70 + pc70[di] * x**(d-di);

	plt.subplot(131);
	plt.plot(inv50);
	plt.plot(p50);	
	plt.subplot(132);
	plt.plot(inv60);
	plt.plot(p60);
	plt.subplot(133);
	plt.plot(inv70);
	plt.plot(p70);
	#plt.plot(inv60);
	#plt.plot(inv70);
	plt.show();

	plt.plot(E50_south);	
	plt.plot(p50*E50_north);
	plt.show();




