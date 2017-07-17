# buoy.py
# A module containing buoyancy-related functions
#=====================================================================

import numpy as np
import matplotlib.pyplot as plt

from diagnostics import diff, extend

#=====================================================================

# footprint
# A function that calculates the PV footprint of the 1L SW solution as produced by RSW_1L.py
def footprint(u_full,v_nd,eta_full,U0_nd,U,Umag,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,AmpF_nd,FORCE,r0,nu,BG,Fpos,ts,period_days,N,Nt,GAUSS):
# This code calculates the PV footprint, i.e. the PV flux convergence defined by
# P = -(div(u*q,v*q)-div((u*q)_av,(v*q)_av)), where _av denotees the time average.
# We will take the time average over one whole forcing period, and subtract this PV flux
# convergence from the PV flux convergence at the end of one period.
# To calculate footprints, we use the full PV and velocity profiles.

	bu = eta_full * u_full;		# Zonal buoy flux
	bv = eta_full * v_nd;		# Meridional buoy flux

	# Next step: taking appropriate derivatives of the fluxes. To save space, we won't define new variables, but instead overwrite the old ones.
	# From these derivatives we can calculate the 
	P = - diff(bu[:,:,0],1,1,dx_nd) - diff(bv[:,:,0],0,0,dy_nd);		# Initialise the footprint with the first time-step
	for ti in range(1,Nt):
		P[:,:] = P[:,:] - diff(bu[:,:,ti],1,1,dx_nd) - diff(bv[:,:,ti],0,0,dy_nd);
	P = P / Nt;
	#P_av = np.trapz(P,T_nd[:Nt],dt_nd,axis=2) / T_nd[Nt-1];
	
	# Normalisation
	#P = P / AmpF_nd**2;
	
	# We are interested in the zonal average of the footprint
	P_xav = np.trapz(P,x_nd[:N],dx_nd,axis=1);
	
	P = extend(P);

	PLOT = True;
	if PLOT:
		plt.figure(1,figsize=[6,13]);
		plt.subplot(321)
		plt.contourf(bu[:,:,0]);
		plt.title('bu');
		plt.colorbar();
		plt.subplot(322);
		plt.contourf(bv[:,:,0]);
		plt.title('bv')
		plt.colorbar();
		plt.subplot(323)
		plt.contourf(-diff(bu[:,:,0],1,1,dx_nd));
		plt.title('-(bu)_x');
		plt.colorbar();
		plt.subplot(324);
		plt.contourf(-diff(bv[:,:,0],0,0,dy_nd));
		plt.title('-(bv)_y');
		plt.colorbar();
		plt.subplot(325)
		plt.contourf(-diff(bu[:,:,0],1,1,dx_nd)-diff(bv[:,:,0],0,0,dy_nd))
		plt.title('P instant')
		plt.colorbar();
		plt.subplot(326)
		plt.contourf(P)
		plt.title('P');
		plt.colorbar();
		plt.show();

	return P, P_xav
