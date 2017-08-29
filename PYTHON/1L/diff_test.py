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

y = np.cos(4 * np.pi * x / L);

y1 = SCHEME(y,2,1,dx);

plt.plot(y1);
plt.show();



