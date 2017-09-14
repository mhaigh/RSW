# Gold_forcing

import numpy as np
import matplotlib.pyplot as plt

OPT = 5;

N = 256;
y = np.linspace(1,N,N);
x = np.linspace(1,N,N);
Amp = 1.0;

j2 = int((y[N-1] - y[0]) / 2) + int(y[0]);
jSo = int((j2 - y[0]) / 2) + int(y[0]);
jNo = int((y[N-1] - j2)/ 2) + j2;


f = np.zeros(N);

#==========================================================================


if OPT == -3:

	h1_restore = np.zeros(N);
	
	je = N; js = 0;
	j_50 = int((je - js) / 2) + js ;
	j_25 = int((j_50 - js) / 2) + js ; j_75 = int((je - j_50) / 2) + j_50;
	j_12 = int((j_25 - js) / 2) + js ; j_37 = int((j_50 - j_25) / 2) + j_25;
	j_62 = int((j_75 - j_50) / 2) + j_50 ; j_87 = int((je - j_75) / 2) + j_75;
	j_44 = int((j_50 - j_37) / 2) + j_37 ; j_56 = int((j_62 - j_50) / 2) + j_50;

	for j in range(js,je):
		h1_restore[j] = 500;
	for j in range(j_44,j_56):
		y = 2.0 * np.pi * (j - j_37) / (j_62 - j_37) 
		restore_shape = 40.0 * (y - y**3 / 6.0 + y**5 / 120.0 - y**7 / 5040.0 + y**9 / 362880.0 - y**11 / 39916800.0 + y**13 / 6227020800.0 - y**15 / 1307674368000 + y**17 / 3.55687428e14);
		h1_restore[j] = 500.0 + restore_shape
	for j in range(js,j_44):
		y = np.pi * (j - js) / (2 * (j_44 - js)); 
		restore_shape = 40.0 * (y - y**3 / 6.0 + y**5 / 120.0 - y**7 / 5040.0 + y**9 / 362880.0 - y**11 / 39916800.0 + y**13 / 6227020800.0 - y**15 / 1307674368000 + y**17 / 3.55687428e14);
		h1_restore[j] = 500.0 + restore_shape
	for j in range(j_56,je):
		y = np.pi * (j - je) / (2 * (je - j_56)); 
		restore_shape = 40.0 * (y - y**3 / 6.0 + y**5 / 120.0 - y**7 / 5040.0 + y**9 / 362880.0 - y**11 / 39916800.0 + y**13 / 6227020800.0 - y**15 / 1307674368000 + y**17 / 3.55687428e14);
		h1_restore[j] = 500.0 + restore_shape
	
	plt.plot(h1_restore);
	plt.show();

if OPT == -2:
	h1_restore = np.zeros(N);
	
	je = N; js = 0;
	j_50 = int((je - js) / 2) + js ;
	j_25 = int((j_50 - js) / 2) + js ; j_75 = int((je - j_50) / 2) + j_50;
	j_12 = int((j_25 - js) / 2) + js ; j_37 = int((j_50 - j_25) / 2) + j_25;
	j_62 = int((j_75 - j_50) / 2) + j_50 ; j_87 = int((je - j_75) / 2) + j_75;
	for j in range(js,je):
		h1_restore[j] = 500;
	for j in range(j_37,j_62):
		y =2.0 * np.pi * (j - j_37) / (j_62 - j_37) 
		restore_shape = 40.0 * (y - y**3 / 6.0 + y**5 / 120.0 - y**7 / 5040.0 + y**9 / 362880.0 - y**11 / 39916800.0 + y**13 / 6227020800.0 - y**15 / 1307674368000 + y**17 / 3.55687428e14);
		h1_restore[j] = 500.0 + restore_shape
	
	plt.plot(h1_restore);
	plt.show();



if OPT == -1:
	h1_restore = np.zeros(N);
	
	je = N; js = 0;
	j_50 = int((je - js) / 2) + js ;
	j_25 = int((j_50 - js) / 2) + js ; j_75 = int((je - j_50) / 2) + j_50;
	j_12 = int((j_25 - js) / 2) + js ; j_37 = int((j_50 - j_25) / 2) + j_25;
	j_62 = int((j_75 - j_50) / 2) + j_50 ; j_87 = int((je - j_75) / 2) + j_75;
	for j in range(js,je):
		h1_restore[j] = 500;
	for j in range(j_37,j_62):
		restore_shape = 40.0 * np.sin(2 * np.pi * (j - j_37) / (j_62 - j_37));
		h1_restore[j] = 500.0 + restore_shape
	
	plt.plot(h1_restore);
	plt.show();

if OPT == 0:
	h_restore1 = np.zeros(N);
	h_restore2 = np.zeros(N);
	for j in range(0,N):
		h_restore1[j] = 500 - 1 * np.tanh((j-N/2));
		h_restore2[j] = 500 - 1 * np.tanh((j-N/2)/10.0);

	plt.subplot(121);
	plt.plot(y,h_restore1);
	plt.subplot(122);
	plt.plot(y,h_restore2);
	plt.show();

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

# Igor's wind forcing
if OPT == 3:

	a = N-1;
	b = N-1;
	beta = 0.8;
	Amp = 1.0;
	wind_asym1 = 1.5;

	taux1 = np.zeros(N);
	taux2 = np.zeros(N);
	for j in range(0,N):
		x0 = beta * y[j]
		Pr= (beta + Amp * (- beta + x0 / b)) * wind_asym1;			 	# The asymmetry of the wind
		taux1[j] = 0.5 * (1.0 - np.cos(2.0 * np.pi * y[j] / N));
		taux2[j] = 0.5 * (1.0 - np.cos(2.0 * np.pi * y[j] / N) + Pr);

	i1 = np.argsort(-taux1);
	i2 = np.argsort(-taux2);
	
	print('Buoyancy shift = ' + str(int(N / 32.0)))
	print('Wind shift = ' + str(i2[0] - i1[0]));

	plt.subplot(121);
	plt.plot(taux1);
	plt.subplot(122);
	plt.plot(taux2);	
	plt.show();

# Pavel's wind forcing
if OPT == 4:
	A = 0.9;
	B = 0.2;
	
	tau0 = 0.8;

	W = np.zeros((N,N));

	for j in range(0,N):
		for i in range(0,N):
			if y[j] <= B * x[i]:
				W[j,i] = - (np.pi * tau0 * A / N) * np.sin(np.pi * (N + y[j]) / (N + B * x[i]));
			else:
				W[j,i] = (np.pi * tau0 / (A * N)) * np.sin(np.pi * (y[j] - B * x[i]) / (N - B * x[i]));

	plt.contourf(W);
	plt.colorbar();
	plt.show();


if OPT == 5:
	
	N = 400;
	y = np.linspace(1,N,N);
	x = np.linspace(1,N,N);

	a = N-1;
	b = N-1;

	m = 0.2
	yf = y[N/2+1];	
	
	beta = 0.8;
	Amp = 1.0;

	wind_asym1 = 1.5;
	
	theta = np.pi / 8;	

	taux = np.zeros((N,N));
	tauy = np.zeros((N,N));
	for j in range(0,N):
		for i in range(0,N):
			x0 = beta * y[j]
			#Pr= (beta + Amp * (- beta + x0 / b)) * wind_asym1;			 	# The asymmetry of the wind
			taux[j,i] = 0.5 * (1.0 + np.cos(2.0 * np.pi *(m * x[i] - y[j] - yf) / N));
			tauy[j,i] = m * taux[j,i];


	Nv = 15;
	yv = np.linspace(1,Nv,Nv);
	xv = np.linspace(1,Nv,Nv);

	taux_vec = np.zeros((Nv,Nv));
	tauy_vec = np.zeros((Nv,Nv));
	for j in range(0,Nv):
		for i in range(0,Nv):
			x0 = beta * yv[j]
			#Pr= (beta + Amp * (- beta + x0 / b)) * wind_asym1;			 	# The asymmetry of the wind
			taux_vec[j,i] = 0.5 * (1.0 + np.cos(1.0 * 2.0 * np.pi *(m * xv[i] - yv[j] - yf) / Nv));
			tauy_vec[j,i] = m * taux_vec[j,i];

	

	plt.subplot(121);	
	plt.quiver(xv,yv,taux_vec,tauy_vec);
	plt.subplot(122);
	plt.contourf(taux);
	plt.colorbar();
	plt.show();
	

if OPT == 6:
	
	N = 20;
	y = np.linspace(1,N,N);
	x = np.linspace(1,N,N);

	a = N-1;
	b = N-1;
	beta = 0.8;
	Amp = 1.0;
	wind_asym1 = 1.5;
	
	theta = np.pi / 8;	

	taux = np.zeros((N,N));
	for j in range(0,N):
		x0 = beta * y[j]
		Pr= (beta + Amp * (- beta + x0 / b)) * wind_asym1;			 	# The asymmetry of the wind
		taux[j,:] = 0.5 * (1.0 - np.cos(2.0 * np.pi * y[j] / N) + Pr);
		
	Rtaux = np.zeros((N,N));
	Rtauy = np.zeros((N,N));
	for j in range(0,N):
		Rtaux[j,:] = np.cos(theta) * taux[j,:];
		Rtauy[j,:] = np.sin(theta) * taux[j,:];

	plt.subplot(121);
	plt.contourf(Rtaux);
	plt.subplot(122);	
	plt.contourf(Rtauy);
	plt.show();
	
	plt.subplot(121);
	plt.quiver(x,y,taux,np.zeros((N,N)));
	plt.subplot(122);
	plt.quiver(x,y,Rtaux,Rtauy);
	plt.show();
		
	





	
