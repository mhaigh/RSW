# diff_test.py
#====================================================

import diagnostics
import numpy as np
import matplotlib.pyplot as plt

#====================================================

#SCHEME = diagnostics.diff;
#SCHEME = diagnostics.diff_center_4th;
SCHEME = diagnostics.diff_fwd_2nd;

N=128;
L=3840.0;

x = np.linspace(0,L,N+1);
x = x[0:N];
dx = x[1] - x[0];
print(x[0],x[N-1]);

y0 = np.zeros((N,N));
y1 = np.zeros((N,N));
for j in range(0,N):
	y0[j,:] = np.cos(4 * np.pi * x[j] / L);	
	y1[:,j] = np.cos(4 * np.pi * x[j] / L);

y0 = SCHEME(y0,0,1,dx);
y1 = SCHEME(y1,1,1,dx);

plt.subplot(121);
plt.plot(y0[:,0]);
plt.subplot(122);
plt.plot(y1[0,:]);
plt.show();



