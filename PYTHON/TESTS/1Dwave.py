# 1Dwave
#=============================================================

# 1-D wave equation plunger test

#=============================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.sparse.linalg import eigsh

#=============================================================

#=============================================================
#=============================================================

# diff
def diff(f):
	
	N = len(f);
	df = np.zeros(N);

	for i in range(0,N-1):
		df[i] = f[i+1] - f[i];
	df[N-1] = f[0] - f[N-1];	

	return df

#=============================================================
#=============================================================

I = np.complex(0.0,1.0);

N = 512;

g = 9.81;
H = 4000.0;

# Space parameters
L = 3840.0 * 1000.0
x = np.linspace(-L/2,L/2,N+1);
x0 = 180.0 * 1000.0
K = np.fft.fftfreq(N,L/N);
dx = x[1] - x[0];

# Time parameters
period_days = 120.0;
period = 3600. * 24. * period_days;		# Periodicity of plunger (s)
omega = 1. / (period);          		# Frequency of plunger, once every 50 days (e-6) (s-1)
Nt = 200;								# Number of time samples
T = np.linspace(0,period,Nt+1);			# Array of time samples across one forcing period (s)
dt = T[1] - T[0];						# Size of the timestep (s)
ts = 100; 								# index at which the time-snapshot is taken
t = T[ts];
c = np.sqrt(g*H);	
time_dep = np.cos(2. * np.pi * omega * T);
# Plunger
A = 10.0e-7
F_dcts = np.zeros(N);
F_cts = np.zeros(N);
for i in range(0,N):
	if abs(x[i]) < x0:
		F_dcts[i] = A * np.cos((np.pi / 2) * x[i] / x0);
		F_cts[i] = 0.5 * A * (1 + np.cos(np.pi * x[i] / x0));

# Fourier transform of plunger
F_tilde_dcts = np.zeros(N);
F_tilde_cts = np.zeros(N);
a = np.pi / (2 * x0);
for i in range(1,N):
 	b = 2.0 * np.pi * K[i];
	F_tilde_dcts[i] = A * (2 * a * np.sin(a*x0) * np.cos(b*x0) - 2 * b * np.cos(a*x0) * np.sin(b*x0)) / (a**2 - b**2);
	F_tilde_cts[i] = A * np.sin(2*np.pi*K[i]*x0) / (2 * np.pi * K[i] * (1 - 4 * x0**2 * K[i]**2));

# Plot forcing terms
plt.subplot(121);
plt.plot(F_cts);
plt.subplot(122);
plt.plot(F_tilde_cts);
plt.show();

# Solution is given by
# u_tilde = - F_tilde / (4 * np.pi**2 * (omega**2 + c**2 * k**2));
u_tilde = np.zeros(N);
for i in range(0,N):
	u_tilde[i] = - F_tilde_cts[i] / (4 * np.pi**2 * (omega**2 + c**2 * K[i]**2));
u = np.zeros((N,Nt),dtype=complex);
for ti in range(0,Nt):
	u[:,ti] = np.fft.ifft(u_tilde) * time_dep[ti] / dx;

# Solution via finite differences instead
A = np.zeros((N,N));
a = - 4 * np.pi**2 * omega**2;	b = - c**2 / dx**2;
A[0,0] = a - 2.0 * b; A[0,1] = b; A[0,N-1] = b;	
A[N-1,N-1] = a - 2.0 * b; A[N-1,N-2] = b; A[N-1,0] = b;
for i in range(1,N-1):
	A[i,i] = a - 2.0 * b;
	A[i,i+1] = b;
	A[i,i-1] = b;
u = np.linalg.solve(A,F_cts);

# Eigenmodes and frequencies
vec = np.zeros((N,N),dtype=complex);	# N eigenmodes (first index), each N long (second index).
val = np.zeros(N,dtype=float);			# N eigenvalues/frequencies. omega = c * k		
for i in range(0,N):
	vec[:,i] = np.exp(2.0 * I*np.pi * K[i] * x[0:N]);
	val[i] = c * K[i];
val_period_days = 1. / (val * 24.0 * 3600.0);

# Decompose solution snapshot into eigenmodes
theta = np.linalg.solve(vec,u);
theta_abs = np.abs(theta);
dom_index = np.argsort(-theta_abs);	
vec = vec[:]		

# Create a projection of Nm most dominant modes
Nm = 5;
proj = theta[dom_index[0]] * vec[:,dom_index[0]];
for mi in range(1,Nm):
	di = dom_index[mi];
	proj = proj + theta[di] * vec[:,di]
print(val_period_days[dom_index[0:Nm]]);

# Now plot
u = np.real(u);
proj = np.real(proj);
plt.figure(1);
plt.plot(u);	
plt.plot(proj);
plt.show();

# Error of the solution
u_tt = - 4 * np.pi**2 * omega**2 * u;
u_xx = diff(diff(u)) / dx**2
e = u_tt - c**2 * u_xx - F_cts; 

# Plot error
plt.subplot(221);
plt.plot(u_tt);
plt.subplot(222);
plt.plot(c**2*u_xx);
plt.subplot(223);
plt.plot(F_cts*time_dep[ts]);
plt.subplot(224);
plt.plot(e);
plt.show();


