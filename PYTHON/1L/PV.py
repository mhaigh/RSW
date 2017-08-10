# PV
# A module containing PV-related functions
#=====================================================================

import numpy as np
import matplotlib.pyplot as plt

from diagnostics import diff, extend

#=====================================================================

# vortcity
# Calculate potential vorticity
def vort(u_nd,v_nd,eta_nd,u_full,eta_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd):
	
	RV_full = np.zeros((N,N,Nt));
	RV_prime = np.zeros((N,N,Nt));
	for ti in range(0,Nt):
		# Define the relative vorticities (RV_full=RV_BG+RV_prime, can always check this numerically)
		RV_full[:,:,ti] = diff(v_nd[:,:,ti],1,1,dx_nd) - diff(u_full[:,:,ti],0,0,dy_nd);
		RV_prime[:,:,ti] = diff(v_nd[:,:,ti],1,1,dx_nd) - diff(u_nd[:,:,ti],0,0,dy_nd);
	RV_BG = - diff(U0_nd,2,0,dy_nd);	# This is defined outside the loop as it has no time-dependence.
	
	# Now define two of the PVs
	PV_full = np.zeros((N,N,Nt));
	PV_BG = np.zeros(N);
	for j in range(0,N):
		PV_BG[j] = (RV_BG[j] + f_nd[j]) / H0_nd[j];
		for i in range(0,N):
			for ti in range(0,Nt):
				PV_full[j,i,ti] = (RV_full[j,i,ti] + f_nd[j]) / eta_full[j,i,ti];
		
	# Two options to define the PV induced in the forced system: (1) PV_full-PV_BG or (2) use the algebraic def. given in the report.
	PV_prime = np.zeros((N,N,Nt));
	for j in range(1,N-1):
		for i in range(0,N):
			for ti in range(0,Nt):
				PV_prime[j,i,ti] = PV_full[j,i,ti] - PV_BG[j];		# Option 1
	#PV_prime = np.zeros((N,N,Nt));	# Option 2 - keep this commented out, just use it as a check.
	#for j in range(0,N):
	#	for i in range(0,N):
	#		for ti in range(0,Nt):
	#			PV_prime[j,i,ti] = (RV_full[j,i,ti] - f[j]) / (eta_full[j,i,ti]) - (f[j] - RV_BG[j]) / (H0_nd[j]);

	return PV_prime, PV_full, PV_BG

#====================================================

# footprint_1L
# A function that calculates the PV footprint of the 1L SW solution as produced by RSW_1L.py
def footprint_1L(u_full,v_nd,eta_full,PV_full,U0_nd,U,Umag,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,AmpF_nd,FORCE,r0,nu,BG,Fpos,ts,period_days,N,Nt,GAUSS):
# This code calculates the PV footprint, i.e. the PV flux convergence defined by
# P = -(div(u*q,v*q)-div((u*q)_av,(v*q)_av)), where _av denotees the time average.
# We will take the time average over one whole forcing period, and subtract this PV flux
# convergence from the PV flux convergence at the end of one period.
# To calculate footprints, we use the full PV and velocity profiles.

	qu = PV_full * u_full;		# Zonal PV flux
	qv = PV_full * v_nd;		# Meridional PV flux

	# Next step: taking appropriate derivatives of the fluxes. To save space, we won't define new variables, but instead overwrite the old ones.
	# From these derivatives we can calculate the 
	P = - diff(qu[:,:,0],1,1,dx_nd) - diff(qv[:,:,0],0,0,dy_nd);		# Initialise the footprint with the first time-step
	for ti in range(1,Nt):
		P[:,:] = P[:,:] - diff(qu[:,:,ti],1,1,dx_nd) - diff(qv[:,:,ti],0,0,dy_nd);
	P = P / Nt;
	#P_av = np.trapz(P,T_nd[:Nt],dt_nd,axis=2) / T_nd[Nt-1];
	
	# Normalisation
	P = P / AmpF_nd**2;
	
	# We are interested in the zonal average of the footprint
	P_xav = np.trapz(P,x_nd[:N],dx_nd,axis=1);
	
	P = extend(P);
	
	return P, P_xav

#====================================================

# footprintComponents
def footprintComponents(u_nd,v_nd,eta_nd,PV_prime,PV_BG,U0_nd,x_nd,dx_nd,dy_nd,N,Nt):
# A function that calculates the PV footprint of the 1L SW solution in terms of its components, allowing for analysis.
# The function calculates the following terms: (1) uq, (2) Uq, (3) uQ, (4) UQ, (5) vq and (6) vQ. (UQ has zero zonal derivative.)
# The zonal/meridional derivative of the zonal/meridional PV flux is taken, averaged over one forcing period.
# Lastly, the zonal averages are calculated and everything useful returned.

	uq = u_nd * PV_prime;
	vq = v_nd * PV_prime;
	UQ = U0_nd * PV_BG;
	# Other components need to be defined in a loop
	uQ = np.zeros((N,N,Nt));
	Uq = np.zeros((N,N,Nt));
	vQ = np.zeros((N,N,Nt));
	for ti in range(0,Nt):
		for i in range(0,N):
			uQ[:,i,ti] = u_nd[:,i,ti] * PV_BG[:];
			Uq[:,i,ti] = U0_nd[:] * PV_prime[:,i,ti];
			vQ[:,i,ti] = v_nd[:,i,ti] * PV_BG[:];

	# Derivatives (no need to edit operate on UQ) and time-averaging
	uq_av = diff(uq[:,:,0],1,1,dx_nd);
	uQ_av = diff(uQ[:,:,0],1,1,dx_nd);
	Uq_av = diff(Uq[:,:,0],1,1,dx_nd);
	vQ_av = diff(vQ[:,:,0],0,0,dy_nd);
	vq_av = diff(vq[:,:,0],0,0,dy_nd)

	for ti in range(1,Nt):
		uq_av = uq_av + diff(uq[:,:,ti],1,1,dx_nd);
		uQ_av = uQ_av + diff(uQ[:,:,ti],1,1,dx_nd);
		Uq_av = Uq_av + diff(Uq[:,:,ti],1,1,dx_nd);
		vQ_av = vQ_av + diff(vQ[:,:,ti],0,0,dy_nd);
		vq_av = vq_av + diff(vq[:,:,ti],0,0,dy_nd);

	uq_av = uq_av / Nt;
	uQ_av = uQ_av / Nt;
	Uq_av = Uq_av / Nt;
	vq_av = vq_av / Nt;
	vQ_av = vQ_av / Nt;

	# Zonal averaging
	uq_xav = np.trapz(uq_av,x_nd[:N],dx_nd,axis=1);
	uQ_xav = np.trapz(uQ_av,x_nd[:N],dx_nd,axis=1);
	Uq_xav = np.trapz(Uq_av,x_nd[:N],dx_nd,axis=1);
	vq_xav = np.trapz(vq_av,x_nd[:N],dx_nd,axis=1);
	vQ_xav = np.trapz(vQ_av,x_nd[:N],dx_nd,axis=1);

	return uq_av, uQ_av, Uq_av, UQ, vq_av, vQ_av, uq_xav, uQ_xav, Uq_xav, vq_xav, vQ_xav;

#====================================================

# footprintComponentsPlot
# A function that plots the footprint components.
def footprintComponentsPlot(P,P_xav,x_nd,y_nd,ts,T_nd,dx_nd,dy_nd,N,Nt):

	plt.subplot(121);
	plt.contourf(v_nd[:,:,ts]);
	plt.subplot(122);
	plt.contourf(PV_prime[:,:,ts]);
	plt.show();		

	plt.figure(1);

	plt.subplot(321);
	plt.contourf(uq[:,:,ts]);
	plt.title('uq');
	plt.colorbar();

	plt.subplot(322);
	plt.contourf(uQ[:,:,ts]);
	plt.title('uQ');
	plt.colorbar();

	plt.subplot(323);
	plt.contourf(Uq[:,:,ts]);
	plt.title('Uq');
	plt.colorbar();

	plt.subplot(324);
	plt.plot(UQ,y_nd);
	plt.title('UQ');

	plt.subplot(325);
	plt.contourf(vq[:,:,ts]);
	plt.title('vq');
	plt.colorbar();

	plt.subplot(326);
	plt.contourf(vQ[:,:,ts]);
	plt.title('vQ');
	plt.colorbar();

	plt.show();

	# Now overwrite the values with their derivatives
	uq = diff(uq[:,:,ts],1,1,dx_nd);
	uQ = diff(uQ[:,:,ts],1,1,dx_nd);
	Uq = diff(Uq[:,:,ts],1,1,dx_nd);
	vQ = diff(vQ[:,:,ts],0,0,dy_nd);
	vq = np.zeros((N,N,Nt));
	for ti in range(0,Nt):
			vq[:,:,ti] = diff(vq[:,:,ti],0,0,dy_nd);	

	plt.figure(2);

	plt.subplot(321);
	plt.contourf(uq);
	plt.title('uq');
	plt.colorbar();

	plt.subplot(322);
	plt.contourf(uQ);
	plt.title('uQ');
	plt.colorbar();

	plt.subplot(323);
	plt.contourf(Uq);
	plt.title('Uq');
	plt.colorbar();

	plt.subplot(324);
	plt.contourf(vq[:,:,ts]);
	plt.title('vq');
	plt.colorbar();

	plt.subplot(325);
	plt.contourf(vQ);
	plt.title('vQ');
	plt.colorbar();

	plt.show();

	# It can be seen that vQ and uQ are relatively small. Let's look at zonal averages instead.
	uq_av = np.trapz(uq,x_nd[:N],dx_nd,axis=1);
	vQ_av = np.trapz(vQ,x_nd[:N],dx_nd,axis=1);
	uQ_av = np.trapz(uQ,x_nd[:N],dx_nd,axis=1);
	Uq_av = np.trapz(Uq,x_nd[:N],dx_nd,axis=1);
	vq_av = np.zeros((N,Nt));
	for ti in range(0,Nt):
		vq_av[:,ti] = np.trapz(vq[:,:,ti],x_nd[:N],dx_nd,axis=1);
		
	
	plt.figure(3);

	plt.subplot(321);
	plt.contourf(uq);
	plt.title('uq');
	plt.colorbar();
	plt.subplot(322);	
	plt.plot(uq_av,y_nd);

	plt.subplot(323);
	plt.contourf(vq[:,:,ts]);
	plt.title('vq');
	plt.colorbar();
	plt.subplot(324);	
	plt.plot(vq_av,y_nd);

	plt.subplot(325);
	plt.contourf(Uq);
	plt.title('Uq');
	plt.colorbar();
	plt.subplot(326);	
	plt.plot(Uq_av,y_nd);
	
	plt.show()

	plt.figure(4);
	plt.subplot(221);
	plt.contourf(vq[:,:,20]);
	plt.subplot(222);
	plt.contourf(vq[:,:,40]);
	plt.subplot(223);
	plt.contourf(vq[:,:,60]);
	plt.subplot(224);
	plt.contourf(vq[:,:,100]);
	plt.show()
	
	
#====================================================

# equivEddyFlux
# A function that calculates the equivalent eddy flux, given a zonally averaged footprint
def EEF(P_xav,y_nd,y0_nd,dy_nd,omega_nd,N):
# The code works by calculating six integrals (three each side of the forcing) that make up each component of the equivalent eddy flux:
# int1_north/south = int_{y > / < y0} P_xav dy;
# int2_north/south = int_{y > / < y0} |y| |P_xav| dy;
# int3_north/south = int_{y > / < y0} |P_xav| dy.

	# Define two y arrays, with all gridpoints north and south of the forcing location.
	# Define two corresponding P_xav arrays
	y_north = [];
	y_south = [];
	P_north = [];
	P_south = [];
	for j in range(0,N):
		if y_nd[j] > y0_nd:
			y_north.append(y_nd[j]);
			P_north.append(P_xav[j]);
		else:
			y_south.append(y_nd[j]);
			P_south.append(P_xav[j]);

	# Convert all to numpy arrays.
	P_north = np.array(P_north);	
	P_south = np.array(P_south);
	y_north = np.array(y_north);
	y_south = np.array(y_south);

	Pabs = abs(P_xav);
	yabs = abs(y_nd);
	
	Pabs_north = abs(P_north);
	Pabs_south = abs(P_south);
	yabs_north = abs(y_north);
	yabs_south = abs(y_south);

	# Now calculate the 6 integrals
	int1_north = np.trapz(P_north,y_north,dy_nd);
	int1_south = np.trapz(P_south,y_south,dy_nd);
	int2_north = np.trapz(Pabs_north*yabs_north,y_north,dy_nd);
	int2_south = np.trapz(Pabs_south*yabs_south,y_south,dy_nd);
	int3_north = np.trapz(Pabs_north,y_north,dy_nd);
	int3_south = np.trapz(Pabs_south,y_south,dy_nd);	

	EEF_north = int1_north * int2_north / int3_north;
	EEF_south = int1_south * int2_south / int3_south;
	
	EEF = (EEF_north - EEF_south) * omega_nd;	
	
	return EEF, EEF_north, EEF_south;

	
 
	
	



