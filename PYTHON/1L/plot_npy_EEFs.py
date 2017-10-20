import output_read
import matplotlib.pyplot as plt

# This line reads the numpy array file to output the full EEF and each of its components.
EEF_50, uq_50, Uq_50, uQ_50, vq_50, vQ_50 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/k0/EEF_y0_om50.npy');
EEF_60, uq_60, Uq_60, uQ_60, vq_60, vQ_60 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/k0/EEF_y0_om60.npy');
EEF_70, uq_70, Uq_70, uQ_70, vq_70, vQ_70 = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/high_res/k0/EEF_y0_om70.npy');

plt.plot(EEF_50);
plt.plot(EEF_60);
plt.plot(EEF_70);
plt.show();


plt.subplot(321);
plt.plot(EEF_50);
plt.title('EEF');

plt.subplot(322);
plt.plot(uq_50);
plt.title('uq');

plt.subplot(323);
plt.plot(Uq_50);
plt.title('Uq');

plt.subplot(324);
plt.plot(uQ_50);
plt.title('uQ');

plt.subplot(325);
plt.plot(vq_50);
plt.title('vq');

plt.subplot(326);
plt.plot(vQ_50);
plt.title('vQ');
plt.show();
