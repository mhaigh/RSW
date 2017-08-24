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
	#P = P / AmpF_nd**2;
	
	# We are interested in the zonal average of the footprint
	P_xav = np.trapz(P,x_nd[:N],dx_nd,axis=1);
	
	P = extend(P);
	
	return P, P_xav

#====================================================

# footprintComponents
def footprintComponents(u_nd,v_nd,eta_nd,PV_prime,PV_BG,U0_nd,AmpF_nd,x_nd,dx_nd,dy_nd,N,Nt):
# A function that calculates the PV footprint of the 1L SW solution in terms of its components, allowing for analysis.
# The function calculates the following terms: (1) uq, (2) Uq, (3) uQ, (4) UQ, (5) vq and (6) vQ. (UQ has zero zonal derivative.)
# The zonal/meridional derivative of the zonal/meridional PV flux is taken, averaged over one forcing period.
# Lastly, the zonal averages are calculated and everything useful returned.

	uq = u_nd * PV_prime;
	vq = v_nd * PV_prime;
	#plt.contourf(vq[:,:,10]);
	#plt.show();
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

	# Derivatives (no need to operate on UQ) and time-averaging
	uq_tav = diff(uq[:,:,0],1,1,dx_nd);
	uQ_tav = diff(uQ[:,:,0],1,1,dx_nd);
	Uq_tav = diff(Uq[:,:,0],1,1,dx_nd);
	vQ_tav = diff(vQ[:,:,0],0,0,dy_nd);
	vq_tav = diff(vq[:,:,0],0,0,dy_nd);

	for ti in range(1,Nt):
		uq_tav = uq_tav + diff(uq[:,:,ti],1,1,dx_nd);
		uQ_tav = uQ_tav + diff(uQ[:,:,ti],1,1,dx_nd);
		Uq_tav = Uq_tav + diff(Uq[:,:,ti],1,1,dx_nd);
		vQ_tav = vQ_tav + diff(vQ[:,:,ti],0,0,dy_nd);
		vq_tav = vq_tav + diff(vq[:,:,ti],0,0,dy_nd);

	# Divide by the number of time samples and simultaneously normalise by forcing amplitude and apply the negation
	uq_tav = - uq_tav / Nt;
	uQ_tav = - uQ_tav / Nt;
	Uq_tav = - Uq_tav / Nt;
	vq_tav = - vq_tav / Nt;
	vQ_tav = - vQ_tav / Nt;

	# Normalisation by AmpF_nd not needed if normalised quanities are passed into the function.

	P_tav = uq_tav + uQ_tav + Uq_tav + vq_tav + vQ_tav;

	# Zonal averaging 
	uq_xav = np.trapz(uq_tav,x_nd[:N],dx_nd,axis=1);
	uQ_xav = np.trapz(uQ_tav,x_nd[:N],dx_nd,axis=1);
	Uq_xav = np.trapz(Uq_tav,x_nd[:N],dx_nd,axis=1);
	vq_xav = np.trapz(vq_tav,x_nd[:N],dx_nd,axis=1);
	vQ_xav = np.trapz(vQ_tav,x_nd[:N],dx_nd,axis=1);
	
	#P_xav = np.trapz(P_tav,x_nd[:N],dx_nd,axis=1);
	P_xav = uq_xav + uQ_xav + Uq_xav + vq_xav + vQ_xav;
	# Tests confirm that the multiple approaches for calculating P_xav and P_tav yield the same results.

	return P_tav, uq_tav, uQ_tav, Uq_tav, UQ, vq_tav, vQ_tav, P_xav, uq_xav, uQ_xav, Uq_xav, vq_xav, vQ_xav;

#====================================================

# EEF
def EEF(P_xav,y_nd,y0_nd,dy_nd,omega_nd,N):
# A function that calculates the equivalent eddy flux, given a zonally averaged footprint
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
	y_north = np.array(y_north);
	y_south = np.array(y_south);
	P_north = np.array(P_north);	
	P_south = np.array(P_south);

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

	EEF_north = (int1_north * int2_north / int3_north) * omega_nd;
	EEF_south = (int1_south * int2_south / int3_south) * omega_nd;
	
	EEF_array = np.array([EEF_north, EEF_south]);	
	
	return EEF_array;

#====================================================

# EEF_components
def EEF_components(P_xav,uq_xav,uQ_xav,Uq_xav,vq_xav,vQ_xav,y_nd,y0_nd,dy_nd,omega_nd,N):
# This function works in the same way as the EEF function, but instead takes as input each individual component of the zonally averaged footprint
# (of which there are five) and returns the EEF contribution in the north and south from each one.  
# The five footprint components are (1) uq_xav, (2) Uq_xav, (3) uQ_xav, (4) vq_xav and (5) vQ_xav (UQ doesn't contribute - zero zonal derivative).

	# Define two y arrays, with all gridpoints north and south of the forcing location
	# and define 12 corresponding EEF component arrays. Two arrays are needed for the full footprint,
	# so that we can calculate the normalisation integrals - the contribution from each footprint component
	# is normalised by the same normalisation constant: norm1 / norm2.

	# First initialise these arrays as empty lists.
	# y
	y_north = [];
	y_south = [];
	# uq
	uq_north = [];
	uq_south = [];
	# Uq
	Uq_north = [];
	Uq_south = [];
	# uQ
	uQ_north = [];
	uQ_south = [];
	# vq
	vq_north = [];
	vq_south = [];
	# vQ
	vQ_north = [];
	vQ_south = [];
	# P
	P_north = [];
	P_south = [];

	# Now assign the values to the appropriate lists
	for j in range(0,N):
		if y_nd[j] > y0_nd:
			y_north.append(y_nd[j]);
			uq_north.append(uq_xav[j]);
			Uq_north.append(Uq_xav[j]);
			uQ_north.append(uQ_xav[j]);
			vq_north.append(vq_xav[j]);
			vQ_north.append(vQ_xav[j]);
			P_north.append(P_xav[j]);
		else:
			y_south.append(y_nd[j]);
			uq_south.append(uq_xav[j]);
			Uq_south.append(Uq_xav[j]);
			uQ_south.append(uQ_xav[j]);
			vq_south.append(vq_xav[j]);
			vQ_south.append(vQ_xav[j]);
			P_south.append(P_xav[j])

	# Convert all to numpy arrays.
	# y
	y_north = np.array(y_north);
	y_south = np.array(y_south);
	# uq
	uq_north = np.array(uq_north);	
	uq_south = np.array(uq_south);
	# Uq
	Uq_north = np.array(Uq_north);	
	Uq_south = np.array(Uq_south);
	# uQ
	uQ_north = np.array(uQ_north);	
	uQ_south = np.array(uQ_south);
	# vq
	vq_north = np.array(vq_north);	
	vq_south = np.array(vq_south);
	# vQ
	vQ_north = np.array(vQ_north);	
	vQ_south = np.array(vQ_south);
	# P
	P_north = np.array(P_north);	
	P_south = np.array(P_south);
	
	# Absolute values are needed for normalisation - only the full footprint P is need here though.
	yabs = abs(y_nd);
	Pabs = abs(P_xav);

	Pabs_north = abs(P_north);
	Pabs_south = abs(P_south);
	yabs_north = abs(y_north);
	yabs_south = abs(y_south);

	# The normalisation constants: norm_north & norm_south (to multiply the integrals of uq_north/south,...)
	norm1_north = np.trapz(Pabs_north*yabs_north,y_north,dy_nd);
	norm2_north = np.trapz(Pabs_north,y_north,dy_nd);
	norm1_south = np.trapz(Pabs_south*yabs_south,y_south,dy_nd);
	norm2_south = np.trapz(Pabs_south,y_south,dy_nd);
	norm_north = norm1_north / norm2_north;
	norm_south = norm1_south / norm2_south;	

	# Now integrate uq_north/south etc. (overwrite the original variables, not needed), multiply by normalisation constant and forcing frequency.
	# uq
	uq_north = np.trapz(uq_north,y_north,dy_nd) * norm_north * omega_nd;
	uq_south = np.trapz(uq_south,y_south,dy_nd) * norm_south * omega_nd;
	# Uq
	Uq_north = np.trapz(Uq_north,y_north,dy_nd) * norm_north * omega_nd;
	Uq_south = np.trapz(Uq_south,y_south,dy_nd) * norm_south * omega_nd;
	# uQ
	uQ_north = np.trapz(uQ_north,y_north,dy_nd) * norm_north * omega_nd;
	uQ_south = np.trapz(uQ_south,y_south,dy_nd) * norm_south * omega_nd;
	# vq
	vq_north = np.trapz(vq_north,y_north,dy_nd) * norm_north * omega_nd;
	vq_south = np.trapz(vq_south,y_south,dy_nd) * norm_south * omega_nd;
	# vQ
	vQ_north = np.trapz(vQ_north,y_north,dy_nd) * norm_north * omega_nd;
	vQ_south = np.trapz(vQ_south,y_south,dy_nd) * norm_south * omega_nd;

	EEF_north = uq_north + Uq_north + uQ_north + vq_north + vQ_north;
	EEF_south = uq_south + Uq_south + uQ_south + vq_south + vQ_south;	

	# Define a single array to be returned by the function, containing all necessary information.	
	EEF_array = np.zeros((6,2));
	EEF_array[0,:] = [EEF_north, EEF_south]; EEF_array[1,:] = [uq_north,uq_south];
	EEF_array[2,:] = [Uq_north, Uq_south]; EEF_array[3,:] = [uQ_north,uQ_south];
	EEF_array[4,:] = [vq_north, vq_south]; EEF_array[5,:] = [vQ_north,vQ_south];

#= np.array((,[uq_north,uq_south],[Uq_north,Uq_south],[uQ_north,uQ_south],[vq_north,vq_south],[vQ_north,vQ_south]));

	return EEF_array;

	
 
	
	



