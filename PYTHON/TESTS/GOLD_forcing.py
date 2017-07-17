# Gold_forcing

import numpy as np
import matplotlib.pyplot as plt

OPT = 3;

N = 256;
y = np.linspace(1,256,N);
x = np.linspace(1,256,N);
Amp = 1.0;

j2 = int((y[N-1] - y[0]) / 2) + int(y[0]);
jSo = int((j2 - y[0]) / 2) + int(y[0]);
jNo = int((y[N-1] - j2)/ 2) + j2;


f = np.zeros(N);

if OPT == 1:
	for j in range(0,jSo+1):
		f[j] = Amp;
	for j in range(jSo+1,jNo):
		f[j] = Amp * (2 * (y[jNo] - y[j]) / (y[jNo] - y[jSo+1]) - 1);
	for j in range(jNo-1,N):
		f[j] = - Amp;

	plt.plot(f,y);
	plt.show();

if OPT == 2:
	for j in range(0,jSo+1):
		f[j] = Amp * (y[j] - y[0]) / (y[jSo] - y[0]);
	for j in range(jSo+1,jNo):
		f[j] = Amp * (2 * (y[jNo] - y[j]) - (y[jNo] - y[jSo])) / (y[jNo] - y[jSo]);
	for j in range(jNo-1,N):
		f[j] = - Amp * (y[N-1] - y[j]) / (y[N-1] - y[jNo-1]);

	plt.plot(f,y);
	plt.show();

if OPT == 3:
	a = 255;
	b = 255;
	beta = 0.8;
	Amp = 1.0;
	wind_asym1 = 0.5;

	taux1 = np.zeros(256);
	taux2 = np.zeros(256);
	for j in range(0,N):
		x0 = beta * y[j]
		Pr= (beta + Amp * (- beta + x0 / b)) * wind_asym1;			 # The asymmetry of the wind
		taux1[j] = 0.5 * (1.0 - np.cos(2.0 * np.pi * y[j] / 256));
		taux2[j] = 0.5 * (1.0 - np.cos(2.0 * np.pi * y[j] / 256) + Pr);

	i1 = np.argsort(-taux1);
	i2 = np.argsort(-taux2);
	print(i1[0],i2[0]);

	plt.subplot(121);
	plt.plot(taux1);
	plt.subplot(122);
	plt.plot(taux2);
	
	plt.show();

