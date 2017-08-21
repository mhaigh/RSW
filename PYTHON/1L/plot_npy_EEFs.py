import output_read
import matplotlib.pyplot as plt

EEF, uq, Uq,uQ, vq, vQ = output_read.npyReadEEF_y0_components('/home/mike/Documents/GulfStream/RSW/DATA/1L/EEFs/EEF_components_70.npy');

plt.subplot(321);
plt.plot(EEF);
plt.title('EEF');

plt.subplot(322);
plt.plot(uq);
plt.title('uq');

plt.subplot(323);
plt.plot(Uq);
plt.title('Uq');

plt.subplot(324);
plt.plot(uQ);
plt.title('uQ');

plt.subplot(325);
plt.plot(vq);
plt.title('vq');

plt.subplot(326);
plt.plot(vQ);
plt.title('vQ');
plt.show();
