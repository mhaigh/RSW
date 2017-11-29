import numpy as np
import matplotlib.pyplot as plt

N = 256;
Ln = 10.;

K = np.linspace(-Ln/2,Ln/2,(N+1));
L = np.linspace(-Ln/2,Ln/2,(N+1))

U = 0.01;

beta = 2.0e-11;

w = np.zeros((N+1,N+1));

for i in range(0,N+1):
	if K[i] == 0:
		w[i] = 0;
	else:
		w[i] = U * K[i] - beta / (K[i]);

plt.plot(K,w);
plt.show();



