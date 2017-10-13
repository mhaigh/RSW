# diff_test.py
#====================================================

import PV
import diagnostics
import numpy as np
import matplotlib.pyplot as plt
import plotting

from inputFile import *

#====================================================

# This test code computes the various components of the footprint at different stages of the operation.
# TEST1: Do zonal derivative and averaging cause uq to have small contribution?
# TEST2: What causes vq behaviour?

#====================================================
# Load the solution, this has already been normalised.
u_nd = np.load('u_nd.npy');
v_nd = np.load('v_nd.npy');
eta_nd = np.load('eta_nd.npy');
	
eta_full = np.zeros((N,N,Nt));
u_full = np.zeros((N,N,Nt));
for j in range(0,N):
	eta_full[j,:,:] = eta_nd[j,:,:] + H0_nd[j];
	u_full[j,:,:] = u_nd[j,:,:] + U0_nd[j];

PV_prime, PV_full, PV_BG = PV.potentialVorticity(u_nd,v_nd,eta_nd,u_full,eta_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd);
uq, Uq, uQ, UQ, vq, vQ = PV.fluxes(u_nd,v_nd,U0_nd,PV_prime,PV_BG,N,Nt);

P, P_uq, P_uQ, P_Uq, P_vq, P_vQ, P_xav, P_uq_xav, P_uQ_xav, P_Uq_xav, P_vq_xav, P_vQ_xav = PV.footprintComponents(uq,Uq,uQ,vq,vQ,x_nd,T_nd,dx_nd,dy_nd,N,Nt);
plotting.footprintComponentsPlot(uq,Uq,uQ,vq,vQ,P,P_uq,P_Uq,P_uQ,P_vq,P_vQ,P_xav,P_uq_xav,P_uQ_xav,P_Uq_xav,P_vq_xav,P_vQ_xav,x_nd,y_nd,N,Nt);
plotting.plotPrimaryComponents(P_uq,P_vq,P_uq_xav,P_vq_xav,x_nd,y_nd,FORCE,BG,Fpos,N);

# Why does the contribution from u^prime * q^prime disappear?
# Is it because <d/dx(u^prime * q^prime)> = c * (u_prime * q^prime)?
# This is periodic in x, and only need be evaluated at the endpoints in x.


TEST1 = False;
if TEST1:
	uq_tav = diagnostics.timeAverage(uq,T_nd,Nt); 		# Take the time-average
	uq_tav_xav = np.zeros(N);
	for j in range(0,N):
		uq_tav_xav[j] = - (uq_tav[j,N-1] - uq_tav[j,0]);

	plt.plot(uq_tav_xav);
	plt.plot(P_uq_xav);
	plt.show();

# We conclude:
# The zonal derivative and the zonal average operators cancel each other out.
# So the contribution from the zonal PV flux is just the difference between uq at either end of the domain,
# but the domain is periodic so this is small!
# The average of the derivative of a period function is zero.

TEST2 = True;
if TEST2:	

	v_nd = v_nd * AmpF_nd;
	vq_full = v_nd * PV_full;
	#for j in range(0,N):
#		vQ[j,:,:] = v_nd[j,:,:] * PV_BG[j];

	plt.subplot(131);
	plt.contourf(vq_full[:,:,ts]);
	plt.colorbar();
	plt.subplot(132);
	plt.contourf(vQ[:,:,ts]);
	plt.colorbar();
	plt.subplot(133);
	plt.contourf(vq[:,:,ts]);
	plt.colorbar();
	plt.show();






















	
			
