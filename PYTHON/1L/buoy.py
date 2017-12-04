# buoy.py
# A module containing buoyancy-related functions
#=====================================================================

import numpy as np
import matplotlib.pyplot as plt

from diagnostics import diff, extend, timeAverage

#=====================================================================

# footprint
# A function that calculates the thickness footprint of the 1L SW solution as produced by RSW_1L.py
def footprint(u_full,v_nd,eta_full,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,N,Nt):

	bu = eta_full * u_full;		# Zonal buoy flux
	bv = eta_full * v_nd;		# Meridional buoy flux

	bu = timeAverage(bu,T_nd,Nt);
	bv = timeAverage(bv,T_nd,Nt);

	# Calculate the footprint.
	P = - diff(bu,1,1,dx_nd) - diff(bv,0,0,dy_nd);

	P = extend(P);

	# We are interested in the zonal average of the footprint
	P_xav = np.trapz(P,x_nd,dx_nd,axis=1);

	PLOT = True;
	if PLOT:
		plt.figure(1,figsize=[6,13]);
		plt.subplot(321)
		plt.contourf(bu);
		plt.title('bu');
		plt.colorbar();
		plt.subplot(322);
		plt.contourf(bv);
		plt.title('bv')
		plt.colorbar();
		plt.subplot(323)
		plt.contourf(-diff(bu,1,1,dx_nd));
		plt.title('-(bu)_x');
		plt.colorbar();
		plt.subplot(324);
		plt.contourf(-diff(bv,0,0,dy_nd));
		plt.title('-(bv)_y');
		plt.colorbar();
		plt.subplot(325)
		plt.contourf(P)
		plt.title('P');
		plt.colorbar();
		plt.subplot(326);
		plt.plot(P_xav,y_nd);
		plt.tight_layout();
		plt.show();

	return P, P_xav


