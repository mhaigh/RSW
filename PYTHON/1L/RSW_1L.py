# RSW_visc_1L
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

from inputFile_1L import *

# 1L SW Solver
#====================================================
#====================================================

# Forcing
if FORCE_TYPE == 'CTS':
	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_cts(x,y,K,y0,r0,N,FORCE,AmpF,g,f,f0,U,L,dx,dy);
elif FORCE_TYPE == 'DCTS':
	F1_nd, F2_nd, F3_nd, Ftilde1_nd, Ftilde2_nd, Ftilde3_nd = forcing_1L.forcing_dcts(x,y,K,y0,r0,N,FORCE,AmpF,g,f,f0,U,L,dx,dy);
else:
	sys.exit('ERROR: Invalid forcing option selected.');

#F1_nd, F2_nd = forcing_1L.F12_from_F3(F3_nd,f_nd,dx_nd,dy_nd,N);
#F3_nd = forcing_1L.F3_from_F1(F1_nd,f_nd,y_nd,dy_nd,N);
#diagnostics.forcingPlots(x_nd,y_nd,Ro*F1_nd,F2_nd,F3_nd,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N);


# Coefficients
a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
# Solver
if BC == 'NO-SLIP':
	solution = solver.NO_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2);
if BC == 'FREE-SLIP':
	solution = solver.FREE_SLIP_SOLVER2(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2);

utilde_nd, vtilde_nd, etatilde_nd = solver.extractSols(solution,N,N2,BC);
u_nd, v_nd, eta_nd = solver.SPEC_TO_PHYS(utilde_nd,vtilde_nd,etatilde_nd,T_nd,dx_nd,omega_nd,N);

u_nd = np.real(u_nd);
v_nd = np.real(v_nd);
eta_nd = np.real(eta_nd);

# Normalise all solutions by the (non-dimensional) forcing amplitude. 
u_nd = u_nd / AmpF_nd;
v_nd = v_nd / AmpF_nd;
eta_nd = eta_nd / AmpF_nd;

# In order to calculate the vorticities/energies of the system, we require full (i.e. BG + forced response) u and eta
eta_full = np.zeros((N,N,Nt));
u_full = np.zeros((N,N,Nt));
for j in range(0,N):
	eta_full[j,:,:] = eta_nd[j,:,:] + H0_nd[j];
	u_full[j,:,:] = u_nd[j,:,:] + U0_nd[j];


#np.save('u_nd.npy',u_nd);
#np.save('v_nd.npy',v_nd);
#np.save('eta_nd.npy',eta_nd);

#====================================================

# Energy

if doEnergy:
	KE_BG, KE_BG_tot, PE_BG, PE_BG_tot = energy.energy_BG(U0_nd,H0_nd,Ro,y_nd,dy_nd,N);
	E_BG_tot = KE_BG_tot + PE_BG_tot;
	#print (KE_BG_tot, PE_BG_tot);
	tii = 10
	KE, KE_tot = energy.KE(u_full[:,:,tii],v_nd[:,:,tii],eta_full[:,:,tii],x_nd,y_nd,dx_nd,dy_nd,N);
	PE, PE_tot = energy.PE(eta_full[:,:,tii],Ro,x_nd,y_nd,dx_nd,dy_nd,N);
	E = KE + PE;
	E_tot = KE_tot + PE_tot;

	plt.contourf(KE);
	plt.show();
	#print(E_tot-E_BG_tot);

#====================================================

# Error - if calculated, should be done before real part of solution is taken
if errorPhys:
	e1, e2, e3 = diagnostics.error(u_nd,v_nd,eta_nd,dx_nd,dy_nd,dt_nd,U0_nd,H0_nd,Ro,gamma_nd,Re,f_nd,F1_nd,F2_nd,F3_nd,T_nd,ts,omega_nd,N);
	e = np.sqrt((e1**2 + e2**2 + e3**2) / 3.0);
	print 'Error = ' + str(e) + '. Error split = ' + str(e1) + ', ' + str(e2) + ', ' + str(e3);
if errorSpec:
	error_spec = np.zeros((3,N));	# An array to save the spectral error at each wavenumber for each equation.
	for i in range(0,N):
		error_spec[:,i] = diagnostics.specError(utilde_nd[:,i],vtilde_nd[:,i],etatilde_nd[:,i],Ftilde1_nd[:,i],Ftilde2_nd[:,i],Ftilde3_nd[:,i],a1[:,i],a2,a3,a4[i],\
b4,c1[:,i],c2,c3,c4[:,i],f_nd,Ro,K_nd[i],H0_nd,y_nd,dy_nd,N);
	for eq in range(0,3):
		error = sum(error_spec[eq,:]) / N;
		print('Error' + str(int(eq+1)) + '=' + str(error));

#====================================================

# Dimensional solutions
u = u_nd * U;
v = v_nd * U;
eta = eta_nd * chi;

#====================================================

# PV and PV footprints
#====================================================

# Calculate PV fields, footprints and equivalent eddy fluxes (EEFs)
if doPV:
	PV_prime, PV_full, PV_BG = PV.potentialVorticity(u_nd,v_nd,eta_nd,u_full,eta_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd);
	uq, Uq, uQ, UQ, vq, vQ = PV.fluxes(u_nd,v_nd,U0_nd,PV_prime,PV_BG,N,Nt);
	# Keep these next two lines commented out unless testing effects of normalisation.
	# uq, Uq, uQ, UQ, vq, vQ = uq/AmpF_nd**2, Uq/AmpF_nd**2, uQ/AmpF_nd**2, UQ/AmpF_nd**2, vq/AmpF_nd**2, vQ/AmpF_nd**2;
	# PV_prime, PV_full = PV_prime/AmpF_nd, PV_full/AmpF_nd;
	if doFootprints:
		if footprintComponents: 
			P, P_uq, P_uQ, P_Uq, P_vq, P_vQ, P_xav, P_uq_xav, P_uQ_xav, P_Uq_xav, P_vq_xav, P_vQ_xav = PV.footprintComponents(uq,Uq,uQ,UQ,vq,vQ,x_nd,dx_nd,dy_nd,N,Nt);
			#plotting.footprintComponentsPlot(uq,Uq,uQ,vq,vQ,P,P_uq,P_Uq,P_uQ,P_vq,P_vQ,P_xav,P_uq_xav,P_uQ_xav,P_Uq_xav,P_vq_xav,P_vQ_xav,x_nd,y_nd,N,Nt);
			plotting.plotPrimaryComponents(P_uq,P_vq,P_uq_xav,P_vq_xav,x_nd,y_nd,FORCE,BG,Fpos,N);
		else: 
			P, P_xav = PV.footprint(uq,Uq,uQ,UQ,vq,vQ,x_nd,T_nd,dx_nd,dy_nd,N,Nt);			
		if doEEFs:
			if footprintComponents:
				EEF_array = PV.EEF_components(P_xav,P_uq_xav,P_uQ_xav,P_Uq_xav,P_vq_xav,P_vQ_xav,y_nd,y0_nd,dy_nd,omega_nd,N);
				# This returns EEF_array, an array with the following structure:
				# EEF_array = ([EEF_north,EEF_south],[uq_north,uq_south],[Uq_north,Uq_south],[uQ_north,uQ_south],[vq_north,vq_south],[vQ_north,vQ_south]).
				EEF_north = EEF_array[0,0]; EEF_south = EEF_array[0,1];
			else:
				EEF_array = PV.EEF(P_xav,y_nd,y0_nd,dy_nd,omega_nd,N);
				EEF_north = EEF_array[0]; EEF_south = EEF_array[1];
			print(EEF_north, EEF_south);
			
# Buoyancy footprints
#====================================================

# Should these be zero, according to conservation of mass?
#Pb, Pb_xav = buoy.footprint(u_full,v_nd,eta_full,U0_nd,U,Umag,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,AmpF_nd,FORCE,r0,nu,BG,Fpos,ts,period_days,N,Nt,GAUSS);

#====================================================

u_nd = diagnostics.extend(u_nd);
v_nd = diagnostics.extend(v_nd);
eta_nd = diagnostics.extend(eta_nd);

#====================================================

#output.ncSave(utilde_nd,vtilde_nd,etatilde_nd,u_nd,v_nd,eta_nd,x_nd,y_nd,K_nd,T_nd,PV_full,PV_prime,PV_BG,Pq,EEFq,N,Nt);

# Plots
#====================================================
#====================================================

# Call the function that plots the forcing in physical and physical-spectral space.
if plotForcing:
	plotting.forcingPlots(x_nd,y_nd,F1_nd,F2_nd,F3_nd,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N);
	#forcing_1L.forcingInv(Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,x_nd,y_nd,dx_nd,N); # For diagnostic purposes

# Background state plots (inc. BG SSH, BG flow, BG PV)
if plotBG:
	plotting.bgPlots(y_nd,H0_nd,U0_nd,PV_BG);

# Soltuion Plots
if plotSol:
	plotting.solutionPlots(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,True);
	plotting.solutionPlots_save(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,True);
	#plotting.solutionPlotsDim(x,y,u,v,eta,ts,L,FORCE,BG,Fpos,N);

# Plots of PV and zonally averaged PV
if plotPV:
	plotting.pvPlots(PV_full,PV_prime,x_nd,y_nd);
	plotting.pvPlots_save(PV_full,PV_prime,P,P_xav,x_nd,y_nd,ts,FORCE,BG,Fpos,N);
if plotPV_av:
	plotting.PV_avPlots(x_nd,y_nd,PV_prime,PV_BG,PV_full,ts,FORCE,BG,Fpos,N);

# Plots of footprints - may need to edit the source code at times.
if plotFootprint:
	plotting.footprintPlots(x_nd,y_nd,P,P_xav,Fpos,BG,GAUSS,FORCE,nu,r0,period_days,U0_nd,U,N);

# Phase and amplitude
if plotPhaseAmp:
	plotting.solutionPlotsAmp(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);
	plotting.solutionPlotsPhase(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);


