# RESTART.py
#=======================================================

# This code works in the same way as the master file RSW_1L.py, but does not solve the SW model.
# Instead it takes saved solution (in .npy or .nc format), and then can perform all necessary diagnostics.
# Should ensure that the input file is defined in the same way as it was for originally producing the saved npy files.

#=======================================================

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

#=======================================================

#u_nd = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/PAPER1/UNIFORM/u_U0=32.npy');
#v_nd = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/PAPER1/UNIFORM/v_U0=32.npy');
#eta_nd = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/PAPER1/UNIFORM/eta_U0=32.npy');

u_nd = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/PAPER1/GAUSSIAN/y0=0/u_y0=0.npy');
v_nd = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/PAPER1/GAUSSIAN/y0=0/v_y0=0.npy');
eta_nd = np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/PAPER1/GAUSSIAN/y0=0/eta_y0=0.npy');

#u_nd = np.load('/home/mike/Documents/GulfStream/RSW/PYTHON/1L/u_nd.npy');
#v_nd = np.load('/home/mike/Documents/GulfStream/RSW/PYTHON/1L/v_nd.npy');
#eta_nd = np.load('/home/mike/Documents/GulfStream/RSW/PYTHON/1L/eta_nd.npy');

#====================================================


# In order to calculate the vorticities/energies of the system, we require full (i.e. BG + forced response) u and eta
eta_full = np.zeros((N,N,Nt));
u_full = np.zeros((N,N,Nt));
for j in range(0,N):
	eta_full[j,:,:] = eta_nd[j,:,:] + H0_nd[j];
	u_full[j,:,:] = u_nd[j,:,:] + U0_nd[j];

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
			P, P_uq, P_uQ, P_Uq, P_vq, P_vQ, P_xav, P_uq_xav, P_uQ_xav, P_Uq_xav, P_vq_xav, P_vQ_xav = PV.footprintComponents(uq,Uq,uQ,vq,vQ,x_nd,T_nd,dx_nd,dy_nd,N,Nt);
			#plotting.footprintComponentsPlot(uq,Uq,uQ,vq,vQ,P,P_uq,P_Uq,P_uQ,P_vq,P_vQ,P_xav,P_uq_xav,P_uQ_xav,P_Uq_xav,P_vq_xav,P_vQ_xav,x_nd,y_nd,N,Nt);
			#plotting.plotPrimaryComponents(P_uq,P_vq,P_uq_xav,P_vq_xav,x_nd,y_nd,FORCE,BG,Fpos,N);
		else: 
			P, P_xav = PV.footprint(uq,Uq,uQ,UQ,vq,vQ,x_nd,T_nd,dx_nd,dy_nd,N,Nt);	
		np.save('P.npy',P);		
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
	#plotting.solutionPlots(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,True);
	plotting.solutionPlots_save(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,True);
	#plotting.solutionPlotsDim(x,y,u,v,eta,ts,L,FORCE,BG,Fpos,N);

# Plots of PV and zonally averaged PV
#plotting.pvPlots(PV_full,PV_prime,x_nd,y_nd);
plotting.pvPlots_save(PV_full,PV_prime,P,P_xav,x_nd,y_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,True);
if plotPV_av:
	plotting.PV_avPlots(x_nd,y_nd,PV_prime,PV_BG,PV_full,ts,FORCE,BG,Fpos,N);

# Plots of footprints - may need to edit the source code at times.
if plotFootprint:
	plotting.footprintPlots(x_nd,y_nd,P,P_xav,Fpos,BG,GAUSS,FORCE,nu,r0,period_days,U0_nd,U,N);

# Phase and amplitude
if plotPhaseAmp:
	plotting.solutionPlotsAmp(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);
	plotting.solutionPlotsPhase(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);


