# PV
# A module containing PV-related functions
#=====================================================================

import numpy as np
import matplotlib.pyplot as plt

from diagnostics import diff, extend, timeAverage

#=====================================================================

# fluxes
def fluxes(u_nd,v_nd):
# Calculates momentum flux terms.
# For the purpose of EEFs of momentum, we only need three flux terms:
# uu, uv, and vv.
	
	uu = u_nd * u_nd;
	uv = u_nd * v_nd;
	vv = v_nd * v_nd;

	return uu, uv, vv

#====================================================

# footprint
def footprint(uu,uv,vv,x_nd,T_nd,dx_nd,dy_nd,N,Nt):
# A function that calculates the momentum footprint of the 1L SW solution as produced by RSW.py	
	
	# Time-averaging
	uu = timeAverage(uu,T_nd,Nt);
	uv = timeAverage(uv,T_nd,Nt);
	vv = timeAverage(vv,T_nd,Nt);

	# Two footprint terms to calculate
	Mu = - diff(uu,1,1,dx_nd) - diff(uv,0,0,dy_nd);
	Mv = - diff(uv,1,1,dx_nd) - diff(vv,0,0,dy_nd);

	Mu = extend(Mu);
	Mv = extend(Mv);
		
	# We are interested in the zonal average of the footprint
	Mu_xav = np.trapz(Mu,x_nd,dx_nd,axis=1);
	Mv_xav = np.trapz(Mv,x_nd,dx_nd,axis=1);

	return Mu, Mv, Mu_xav, Mv_xav

#====================================================

# EEF_mom
def EEF_mom(Mu_xav,Mv_xav,y_nd,y0_nd,dy_nd,omega_nd,N):
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
	yabs = abs(y_nd - y0_nd);
	
	Pabs_north = abs(P_north);
	Pabs_south = abs(P_south);
	yabs_north = abs(y_north - y0_nd);
	yabs_south = abs(y_south - y0_nd);

	# Now calculate the 6 integrals
	int1_north = np.trapz(P_north,y_north,dy_nd);
	int1_south = np.trapz(P_south,y_south,dy_nd);

	norm1_north = np.trapz(Pabs_north*yabs_north,y_north,dy_nd);
	norm1_south = np.trapz(Pabs_south*yabs_south,y_south,dy_nd);
	norm2_north = np.trapz(Pabs_north,y_north,dy_nd);
	norm2_south = np.trapz(Pabs_south,y_south,dy_nd);

	EEF_north = (int1_north * norm1_north / norm2_north) * omega_nd;
	EEF_south = (int1_south * norm1_south / norm2_south) * omega_nd;
	
	EEF_array = np.array([EEF_north, EEF_south]);

	return EEF_array;


 
	
	



