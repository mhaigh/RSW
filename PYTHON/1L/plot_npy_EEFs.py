import output_read
import matplotlib.pyplot as plt
import numpy as np
from inputFile import *

#BG = 'U0=Gaussian'
BG = 'U0=16'
k = 'k=0';

opt = 'k';

# Read EEF arrays
if opt == 'r':
	EEF_0, uq_0, Uq_0, uQ_0, vq_0, vQ_0 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+BG+'/PV/EEF_PV_y0_r60.npy');
	EEF_1, uq_1, Uq_1, uQ_1, vq_1, vQ_1 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+BG+'/PV/EEF_PV_y0.npy');
	EEF_2, uq_2, Uq_2, uQ_2, vq_2, vQ_2 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+BG+'/PV/EEF_PV_y0_r120.npy');
	l0 = 'r0 = 60 km';
	l1 = 'r0 = 90 km';
	l2 = 'r0 = 120 km';
	# Note that varying r0 means that the lengths of the EEF array will differ
	for ri in range(0,3):
		exec('NN = len(EEF_' + str(ri) + ')');
		N_skip = (N - NN) / 2; 
		exec('y_forced_' + str(ri) + '= y_nd[N_skip:N-N_skip]');
		
elif opt == 'k':
	EEF_0, uq_0, Uq_0, uQ_0, vq_0, vQ_0 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+BG+'/PV/EEF_PV_y0_k0.npy');
	EEF_1, uq_1, Uq_1, uQ_1, vq_1, vQ_1 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+BG+'/PV/EEF_PV_y0.npy');
	EEF_2, uq_2, Uq_2, uQ_2, vq_2, vQ_2 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+BG+'/PV/EEF_PV_y0_k200.npy');
	l0 = 'k = 0';
	l1 = 'k = 100';
	l2 = 'k = 200';
elif opt == 'w':
	EEF_0, uq_0, Uq_0, uQ_0, vq_0, vQ_0 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+BG+'/PV/EEF_PV_y0_om50.npy');
	EEF_1, uq_1, Uq_1, uQ_1, vq_1, vQ_1 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+BG+'/PV/EEF_PV_y0.npy');
	EEF_2, uq_2, Uq_2, uQ_2, vq_2, vQ_2 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+BG+'/PV/EEF_PV_y0_om70.npy');
	l0 = 'T = 50 days';
	l1 = 'T = 60 days';
	l2 = 'T = 70 days';

doMom = False;
if doMom:
	EEF_u_50 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+k+'/'+BG+'/EEF_u_y0_om50.npy');
	EEF_u_50 = EEF_u_50[:,0] - EEF_u_50[:,1];
	EEF_v_50 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+k+'/'+BG+'/EEF_v_y0_om50.npy');
	EEF_v_50 = EEF_v_50[:,0] - EEF_v_50[:,1];
	
	EEF_u_60 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+k+'/'+BG+'/EEF_u_y0_om60.npy');
	EEF_u_60 = EEF_u_60[:,0] - EEF_u_60[:,1];
	EEF_v_60 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+k+'/'+BG+'/EEF_v_y0_om60.npy');
	EEF_v_60 = EEF_v_60[:,0] - EEF_v_60[:,1];

	EEF_u_70 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+k+'/'+BG+'/EEF_u_y0_om70.npy');
	EEF_u_70 = EEF_u_70[:,0] - EEF_u_70[:,1];
	EEF_v_70 = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+k+'/'+BG+'/EEF_v_y0_om70.npy');
	EEF_v_70 = EEF_v_70[:,0] - EEF_v_70[:,1];

NN = len(EEF_0);

if False:
	for j in range(1,100):
		if j%2 == 0:	
			EEF_50[j] = 0.5 * (EEF_50[j+1] + EEF_50[j-1]);
			EEF_60[j] = 0.5 * (EEF_60[j+1] + EEF_60[j-1]);
			EEF_70[j] = 0.5 * (EEF_70[j+1] + EEF_70[j-1]);
			EEF_50[NN-j] = 0.5 * (EEF_50[NN-j+1] + EEF_50[NN-j-1]);
			EEF_60[NN-j] = 0.5 * (EEF_60[NN-j+1] + EEF_60[NN-j-1]);
			EEF_70[NN-j] = 0.5 * (EEF_70[NN-j+1] + EEF_70[NN-j-1]);

N_skip = (N - NN) / 2; # Should always be even
y_forced = y_nd[N_skip:N-N_skip];

print(NN);
print(y_forced);

plt.figure(1);
if opt == 'r':
	plt.plot(y_forced_0,EEF_0,label=l0,linewidth=1.3);
	plt.plot(y_forced_1,EEF_1,label=l1,linewidth=1.3);
	plt.plot(y_forced_2,EEF_2,label=l2,linewidth=1.3);
else:
	#plt.plot(y_forced,EEF_0,label=l0,linewidth=1.3);
	plt.plot(y_forced,EEF_1,label=l1,linewidth=1.3);
	plt.plot(y_forced,EEF_2,label=l2,linewidth=1.3);
plt.xlim(-0.5,0.5);
plt.ylim(0.0,0.09);
plt.grid(b=True, which='both', color='0.65',linestyle='--');
#plt.title(BG+', '+k)
plt.legend();
plt.show();




