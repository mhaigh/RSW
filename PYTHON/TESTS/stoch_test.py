# stoch_test.py
#=======================================================
#=======================================================

# Test file for producing a stochastic process and presenting it in frequency space

#=======================================================

import numpy as np
import matplotlib.pyplot as plt

#=======================================================

Nt = 10;

Po = np.random.poisson(lam=60,size=Nt-1);

# The times
T = np.zeros(Nt);
for ti in range(1,Nt):
	T[ti] = T[ti-1] + Po[ti-1];

S = np.zeros(T[Nt-1]);
for ti in range(0,len(S)):
	if ti in T:
		S[ti] = 1;

Stilde = np.fft.fft(S);

plt.subplot(121)
plt.plot(S);
plt.subplot(122);
plt.plot(Stilde);
plt.show();
