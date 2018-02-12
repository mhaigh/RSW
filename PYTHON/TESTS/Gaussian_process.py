# Gaussian_process.py
#=================================================================

# A function that simulates a stochastic process S = (S[0], S[1], S[2]...)
# which S[i] are Gaussian distributed random variables and the process has
# exponentially distributed correllation, i.e. <S[i],S[i+n]>=exp(-n/tau),
# where tau>0 is some time-scale.

#=================================================================

import numpy as np
import matplotlib.pyplot as plt

#=================================================================

N = 2*120;
period_days = 60.0
tau = 60.0

f = np.exp(-1/tau);
ff = np.sqrt(1-f**2);

T = np.linspace(0,period_days,N+1);
dt = T[1] - T[0];

S = np.zeros(N);

S[0] = 0#np.random.normal(0,1.0)

for ti in range(1,N):
	g = np.random.normal(0,1.0);
	S[ti] = f * S[ti-1] + ff * g;

np.save('time_series',S);

S_tilde = np.fft.fft(S);

plt.subplot(121);
plt.plot(S);
plt.subplot(122);
plt.plot(S_tilde);
plt.show();


