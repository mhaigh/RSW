import output_read
import matplotlib.pyplot as plt
import numpy as np

BG = 'U0=16'
k = 'k=0';

# This line reads the numpy array file to output the full EEF and each of its components.
EEF_50, uq_50, Uq_50, uQ_50, vq_50, vQ_50 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+k+'/'+BG+'/EEF_y0_om50.npy');
EEF_60, uq_60, Uq_60, uQ_60, vq_60, vQ_60 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+k+'/'+BG+'/EEF_y0_om60.npy');
EEF_70, uq_70, Uq_70, uQ_70, vq_70, vQ_70 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/'+k+'/'+BG+'/EEF_y0_om70.npy');

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

N = len(EEF_60);

if False:
	for j in range(1,100):
		if j%2 == 0:	
			EEF_50[j] = 0.5 * (EEF_50[j+1] + EEF_50[j-1]);
			EEF_60[j] = 0.5 * (EEF_60[j+1] + EEF_60[j-1]);
			EEF_70[j] = 0.5 * (EEF_70[j+1] + EEF_70[j-1]);
			EEF_50[N-j] = 0.5 * (EEF_50[N-j+1] + EEF_50[N-j-1]);
			EEF_60[N-j] = 0.5 * (EEF_60[N-j+1] + EEF_60[N-j-1]);
			EEF_70[N-j] = 0.5 * (EEF_70[N-j+1] + EEF_70[N-j-1]);
		

plt.plot(EEF_50);
plt.plot(EEF_60);
plt.plot(EEF_70);
plt.title(BG+', '+k)
plt.show();

plt.subplot(121);
plt.plot(EEF_u_50);
plt.plot(EEF_u_60);
plt.plot(EEF_u_70);
plt.title('EEF_u');
plt.subplot(122);
plt.plot(EEF_v_50);
plt.plot(EEF_v_60);
plt.plot(EEF_v_70);
plt.title('EEF_v');
plt.show();
