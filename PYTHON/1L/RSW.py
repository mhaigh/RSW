# RSW.py
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
import matplotlib.pyplot as plt

from core import solver, PV, momentum, thickness, energy, diagnostics
from output import plotting

from inputFile import *

# 1L SW Solver
#====================================================
#====================================================

def RSW_main():
	# Forcing

	#plotting.forcingPlot_save(x_grid,y_grid,F3_nd[:,0:N],FORCE,BG,Fpos,N);

	#F1_nd, F2_nd, F3_nd = forcing.forcingInv(Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,x_nd,y_nd,dx_nd,N);
	#F1_nd, F2_nd = forcing.F12_from_F3(F3_nd,f_nd,dx_nd,dy_nd,N);
	#F3_nd = forcing.F3_from_F1(F1_nd,f_nd,y_nd,dy_nd,N);
	#plotting.forcingPlots(x_nd[0:N],y_nd,Ro*F1_nd,Ro*F2_nd,F3_nd,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N);

	#sys.exit();
	# Coefficients
	a1,a2,a3,a4,b4,c1,c2,c3,c4 = solver.SOLVER_COEFFICIENTS(Ro,Re,K_nd,f_nd,U0_nd,H0_nd,omega_nd,gamma_nd,dy_nd,N);
	# Solver
	if BC == 'NO-SLIP':
		solution = solver.NO_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2);
	if BC == 'FREE-SLIP':
		solution = solver.FREE_SLIP_SOLVER(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2);
		#solution = solver.FREE_SLIP_SOLVER4(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ro*Ftilde3_nd,N,N2)

	utilde_nd, vtilde_nd, etatilde_nd = solver.extractSols(solution,N,N2,BC);
	u_nd, v_nd, eta_nd = solver.SPEC_TO_PHYS(utilde_nd,vtilde_nd,etatilde_nd,T_nd,dx_nd,omega_nd,N);

	#np.save('u_nd_complex.npy',u_nd[:,:,ts]);
	#np.save('v_nd_complex.npy',v_nd[:,:,ts]);
	#np.save('eta_nd_complex.npy',eta_nd[:,:,ts]);

	#plotting.solutionPlotsPhase(x_grid,y_grid,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);
	#plotting.solutionPlotsAmp(x_grid,y_grid,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);

	#sys.exit();

	u_nd = np.real(u_nd);
	v_nd = np.real(v_nd);
	eta_nd = np.real(eta_nd);
	
	# Normalise all solutions by the (non-dimensional) forcing amplitude. 
	u_nd = u_nd / AmpF_nd;
	v_nd = v_nd / AmpF_nd;
	eta_nd = eta_nd / AmpF_nd;

	# Correlation.
	# Central half?
	cs = N / 4; 
	ce = N - N / 4;
	corr = diagnostics.arrayCorrTime(u_nd[cs:ce,cs:ce,:],v_nd[cs:ce,cs:ce,:]);
	print corr

	# In order to calculate the vorticities/energies of the system, we require full (i.e. BG + forced response) u and eta.
	eta_full = np.zeros((N,N,Nt));
	u_full = np.zeros((N,N,Nt));
	for j in range(0,N):
		eta_full[j,:,:] = eta_nd[j,:,:] + H0_nd[j];
		u_full[j,:,:] = u_nd[j,:,:] + U0_nd[j];

	#np.save('u_nd.npy',u_nd);
	#np.save('v_nd.npy',v_nd);
	#np.save('eta_nd.npy',eta_nd);

	#plt.subplot(121);
	#plt.pcolor(x_grid,y_grid,u_full[:,:,ts],cmap='bwr');
	#plt.colorbar();
	#plt.subplot(122);
	#plt.pcolor(x_grid,y_grid,eta_full[:,:,ts],cmap='bwr');
	#plt.colorbar();
	#plt.show();

	#sys.exit();

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

	# Momentum footprints
	#====================================================
	
	if doMomentum:
		uu, uv, vv = momentum.fluxes(u_nd,v_nd);
		Mu, Mv, Mu_xav, Mv_xav = momentum.footprint(uu,uv,vv,x_nd,T_nd,dx_nd,dy_nd,N,Nt);
		#plotting.MomFootprints(Mu,Mv,Mu_xav,Mv_xav);
		
		Mu = Mu / np.max(abs(Mu));
		Mv = Mv / np.max(abs(Mv));	
	
		plt.subplot(221);
		plt.pcolor(x_grid,y_grid,Mu,cmap='bwr', vmin=-1., vmax=1.);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');

		plt.subplot(222);
		plt.plot(Mu_xav,y_nd);

		plt.subplot(223);
		plt.pcolor(x_grid,y_grid,Mv,cmap='bwr', vmin=-1., vmax=1.);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');

		plt.subplot(224);
		plt.plot(Mv_xav,y_nd);

		plt.tight_layout();
		plt.show();

		EEF_u, EEF_v = momentum.EEF_mom(Mu_xav,Mv_xav,y_nd,y0_nd,y0_index,dy_nd,omega_nd,N);
		
		print(EEF_u, EEF_v);


	# PV and PV footprints
	#====================================================

	# Calculate PV fields, footprints and equivalent eddy fluxes (EEFs)
	if doPV:
		PV_prime, PV_full, PV_BG = PV.potentialVorticity(u_nd,v_nd,eta_nd,u_full,eta_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd,Ro);
		PV_prime1, PV_prime2, PV_prime3 = PV.potentialVorticity_linear(u_nd,v_nd,eta_nd,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd,Ro);
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
			if doEEFs:
				if footprintComponents:
					EEF_array = PV.EEF_components(P_xav,P_uq_xav,P_uQ_xav,P_Uq_xav,P_vq_xav,P_vQ_xav,y_nd,y0_nd,y0_index,dy_nd,omega_nd,N);
					# This returns EEF_array, an array with the following structure:
					# EEF_array = ([EEF_north,EEF_south],[uq_north,uq_south],[Uq_north,Uq_south],[uQ_north,uQ_south],[vq_north,vq_south],[vQ_north,vQ_south]).
					EEF_north = EEF_array[0,0]; EEF_south = EEF_array[0,1];
				else:
					EEF_array = PV.EEF(P_xav,y_nd,y0_nd,y0_index,dy_nd,N);
					EEF_north = EEF_array[0]; EEF_south = EEF_array[1];
					EEF = EEF_north - EEF_south;
				print(EEF_north - EEF_south, EEF_north, EEF_south);
			
	# Buoyancy footprints
	#====================================================
	
	if doThickness:
		# Should these be zero, according to conservation of mass?
		Pb, Pb_xav = thickness.footprint(u_full,v_nd,eta_full,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,N,Nt);

	#====================================================

	#output.ncSave(utilde_nd,vtilde_nd,etatilde_nd,u_nd,v_nd,eta_nd,x_nd,y_nd,K_nd,T_nd,PV_full,PV_prime,PV_BG,Pq,EEFq,N,Nt);

	nrows = 4	
	#plt.figure(1,)
	fig, axes = plt.subplots(nrows=nrows, ncols=3,figsize=(22,7*nrows))
	for row in range(0,nrows):
		plotting.fp_PV_plot(PV_prime,P,P_xav,x_nd,y_nd,ts,N,x_grid,y_grid,row,nrows)
	plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=1.);
	plt.savefig('fig1.png');
	plt.show();
	sys.exit();

	#====================================================

	# Take relevant derivatives
	v_y = np.zeros((N,N,Nt));
	u_y = np.zeros((N,N,Nt));
	u_yy = np.zeros((N,N,Nt));
	for ti in range(0,Nt):
		v_y[:,:,ti] = diagnostics.diff(v_nd[:,:,ti],0,0,dy_nd);
		u_y[:,:,ti] = diagnostics.diff(u_nd[:,:,ti],0,0,dy_nd);
		u_yy[:,:,ti] = diagnostics.diff(u_y[:,:,ti],0,0,dy_nd);
	
	uv1 = v_y * u_y;
	uv2 = v_nd * u_yy;
	uv3 = v_nd * u_y


	plt.subplot(131);
	plt.contourf(x_nd[0:N],y_nd,v_nd[:,:,ts]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.subplot(132);
	plt.contourf(x_nd[0:N],y_nd,u_yy[:,:,ts]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.subplot(133);
	plt.contourf(x_nd[0:N],y_nd,uv2[:,:,ts]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.show()

	plt.subplot(221);
	plt.contourf(x_nd[0:N],y_nd,uv1[:,:,20]);
	plt.colorbar();
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.subplot(222);
	plt.contourf(x_nd[0:N],y_nd,uv1[:,:,100]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.subplot(223);
	plt.contourf(x_nd[0:N],y_nd,uv2[:,:,20]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.subplot(224);
	plt.contourf(x_nd[0:N],y_nd,uv2[:,:,100]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.show();

	# Define initial footprint contributions (include SSH terms later)
	P1 = diagnostics.timeAverage(uv1,T_nd,Nt);
	P2 = diagnostics.timeAverage(uv2,T_nd,Nt);
	P3 = diagnostics.timeAverage(uv3,T_nd,Nt);
	

	P1 = diagnostics.extend(P1);
	P2 = diagnostics.extend(P2);
	P3 = diagnostics.extend(P3);

	plt.subplot(121);
	plt.contourf(x_nd,y_nd,P1);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.subplot(122);	
	plt.contourf(x_nd,y_nd,P2);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.show();

	# Account for H0_nd terms
	#H0_y = diagnostics.diff(H0_nd,2,0,dy_nd);
	#for i in range(0,N):
	#	P1[:,i] = P1[:,i] / H0_nd[:];
	#	P2[:,i] = P2[:,i] / H0_nd[:];
	#	P3[:,i] = P3[:,i] * H0_y[:] / H0_nd[:]**2;
	


	P1 = np.trapz(P1,x_nd,dx_nd,axis=1);
	P2 = np.trapz(P2,x_nd,dx_nd,axis=1);
	P3 = np.trapz(P3,x_nd,dx_nd,axis=1);

	plt.subplot(121);
	#plt.plot(P1,label='P1');
	plt.plot(P2,label='P2');
	plt.legend();
	plt.subplot(122);
	plt.plot(H0_nd*P_xav);
	plt.plot(P1+P2+P3);
	plt.show();

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
		plotting.solutionPlots(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,False);
		#plotting.solutionPlots_save(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,True);
		#plotting.solutionPlotsDim(x,y,u,v,eta,ts,L,FORCE,BG,Fpos,N);
	
	# Plots of PV and zonally averaged PV
	if doPV:
		if plotPV:
			#plotting.pvPlots(PV_full,PV_prime,x_nd,y_nd);
			plotting.pvPlots_save(PV_full,PV_prime,P,P_xav,x_nd,y_nd,ts,FORCE,BG,Fpos,N,U0_str,x_grid,y_grid,True);
		if plotPV_av:
			plotting.PV_avPlots(x_nd,y_nd,PV_prime,PV_BG,PV_full,ts,FORCE,BG,Fpos,N);
		if doFootprints:
			if plotFootprint:
				plotting.footprintPlots(x_nd,y_nd,P,P_xav,Fpos,BG,FORCE,nu,r0,period_days,U0_nd,U,N);
	
	# Phase and amplitude
	if plotPhaseAmp:
		plotting.solutionPlotsAmp(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);
		plotting.solutionPlotsPhase(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);

if __name__ == '__main__':
	RSW_main();


