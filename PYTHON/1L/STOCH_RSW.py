# STOCH_RSW.py

#=======================================================

# This code solves the single-layer shallow water equations (centered-in-space finite difference), with external forcing terms on each of the three equations.
# The equations are solved in a beta-plane zonally periodic channel, with no-normal flow BCs at the northern and southern boundaries.
# The model includes simple linear Rayleigh drag on the ocean's bottom and viscosity.
# Also included is a latitude-dependent zonal BG flow and corresponding geostrophic BG sea-surface height, around which the equations are linearised.
# The original governing equations are simplified by implementing the zonal Fourier transform and assuming time-periodicity.

# This means solving a system of the form:
# a1 * u + a2 * u_yy + a3 * v + a4 * eta = Ro * F1,
# b1 * u + b2 * v + b3 * v_yy + b4 * eta_y = Ro * F2,
# c1 * u + c2 * v + c3 * v_y + c4 * eta = F3,
# where (u,v) is the horizontal velocity vector, eta is the interface height, Fi are the forcings with correpsonding amplitude alphai
# and ai, bi, ci are k- and y-dependent coefficients.

#====================================================

import sys

import numpy as np

import diagnostics
import PV
import buoy
import forcing_1L
import solver
import output
import energy
import plotting

from inputFile import *

# RSW SOLVER STOCHASTIC 
#====================================================
#====================================================

# Forcing
if FORCE_TYPE == 'CTS':
	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_cts(x_nd,y_nd,K_nd,y0_nd,r0_nd,N,FORCE,AmpF_nd,f_nd,f0_nd,dx_nd,dy_nd);
elif FORCE_TYPE == 'DCTS':
	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_dcts(x,y,K,y0,r0,N,FORCE,AmpF,g,f,f0,dx,dy);
else:
	sys.exit('ERROR: Invalid forcing option selected.');
#plotting.forcingPlot_save(x_grid,y_grid,F3_nd[:,0:N],FORCE,BG,Fpos,N);

#F1_nd, F2_nd, F3_nd = forcing_1L.forcingInv(Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,x_nd,y_nd,dx_nd,N);
#F1_nd, F2_nd = forcing_1L.F12_from_F3(F3_nd,f_nd,dx_nd,dy_nd,N);
#F3_nd = forcing_1L.F3_from_F1(F1_nd,f_nd,y_nd,dy_nd,N);
#plotting.forcingPlots(x_nd[0:N],y_nd,Ro*F1_nd,Ro*F2_nd,F3_nd,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N);

#sys.exit();

u_nd = np.zeros((N,N,Nt),dtype=complex);
v_nd = np.zeros((N,N,Nt),dtype=complex);
eta_nd = np.zeros((N,N,Nt),dtype=complex);

S = np.load('time_series.npy');
plt.plot(S);
plt.show();
Om = np.fft.fftfreq(Nt,dt_nd);
S_tilde = np.fft.fft(S) / dt_nd;
plt.plot(S_tilde);
plt.show();

for wi in range(1,Nt):
	print(wi);
	omega_nd = Om[wi];
	# Coefficients
	a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
	# Solver
	if BC == 'NO-SLIP':
		solution = solver.NO_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,S_tilde[wi]*Ro*Ftilde1_nd,S_tilde[wi]*Ro*Ftilde2_nd,S_tilde[wi]*Ftilde3_nd,N,N2);
	if BC == 'FREE-SLIP':
		solution = solver.FREE_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,S_tilde[wi]*Ro*Ftilde1_nd,S_tilde[wi]*Ro*Ftilde2_nd,S_tilde[wi]*Ftilde3_nd,N,N2);

	u_nd[:,:,wi], v_nd[:,:,wi], eta_nd[:,:,wi] = solver.extractSols(solution,N,N2,BC);

u_nd, v_nd, eta_nd = solver.SPEC_TO_PHYS_STOCH(u_nd,v_nd,eta_nd,T_nd,dx_nd,Om,N);

u_nd = np.real(u_nd);
v_nd = np.real(v_nd);
eta_nd = np.real(eta_nd);

# Normalise all solutions by the (non-dimensional) forcing amplitude. 
#u_nd = u_nd / AmpF_nd;
#v_nd = v_nd / AmpF_nd;
#eta_nd = eta_nd / AmpF_nd;

# In order to calculate the vorticities/energies of the system, we require full (i.e. BG + forced response) u and eta
eta_full = np.zeros((N,N,Nt));
u_full = np.zeros((N,N,Nt));
for j in range(0,N):
	eta_full[j,:,:] = eta_nd[j,:,:] + H0_nd[j];
	u_full[j,:,:] = u_nd[j,:,:] + U0_nd[j];

np.save('u_nd.npy',u_nd);
np.save('v_nd.npy',v_nd);
np.save('eta_nd.npy',eta_nd);

sys.exit();


# Soltuion Plots
if plotSol:
	plotting.solutionPlots(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,False);
	plotting.solutionPlots_save(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,True);
	#plotting.solutionPlotsDim(x,y,u,v,eta,ts,L,FORCE,BG,Fpos,N);


