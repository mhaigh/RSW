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
	uBC, etaBC = solver.BC_COEFFICIENTS(Ro,Re,f_nd,H0_nd,dy_nd,N);
	solution = solver.FREE_SLIP_SOLVER2(a1,a2,a3,a4,f_nd,b4,c1,c2,c3,c4,uBC,etaBC,Ro*Ftilde1_nd,Ro*Ftilde2_nd,Ftilde3_nd,N,N2);

utilde_nd, vtilde_nd, etatilde_nd = solver.extractSols(solution,N,N2,BC);
u_nd, v_nd, eta_nd = solver.SPEC_TO_PHYS(utilde_nd,vtilde_nd,etatilde_nd,T_nd,dx_nd,omega_nd,N);

#====================================================

# Error - if calculated, should be done before real part of solution is taken
if errorPhys:
	e1, e2, e3 = diagnostics.error(u_nd,v_nd,eta_nd,dx_nd,dy_nd,dt_nd,U0_nd,H0_nd,Ro,gamma_nd,Re,f_nd,F1_nd,F2_nd,F3_nd,T_nd,ts,omega_nd,N);
	print 'ERROR: ' + str(e1) + ', ' + str(e2) + ', ' + str(e3);
if errorSpec:
	e1_spec = diagnostics.specError(utilde_nd,vtilde_nd,etatilde_nd,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,a1,a2,a3,a4,b4,c1,c2,c3,c4,Ro,K_nd,H0_nd,y_nd,dy_nd,N,1);
	print e1_spec;

#====================================================

u_nd = np.real(u_nd);
v_nd = np.real(v_nd);
eta_nd = np.real(eta_nd);

# In order to calculate the vorticities of the system, we require full (i.e. BG + forced response) u and eta
eta_full = np.zeros((N,N,Nt));
u_full = np.zeros((N,N,Nt));
for j in range(0,N):
	eta_full[j,:,:] = eta_nd[j,:,:] + H0_nd[j];
	u_full[j,:,:] = u_nd[j,:,:] + U0_nd[j];

# Dimensional solutions
u = u_nd * U;
v = v_nd * U;
eta = eta_nd * chi;

#====================================================

# PV and PV footprints
#====================================================

# Calculate PV fields, footprints and equivalent eddy fluxes (EEFs)
if doPV:
	PV_prime, PV_full, PV_BG = PV.vort(u_nd,v_nd,eta_nd,u_full,eta_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd);
	if doFootprints: 
		P, P_xav = PV.footprint_1L(u_full,v_nd,eta_full,PV_full,U0_nd,U,Umag,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,AmpF_nd,FORCE,r0,nu,BG,Fpos,ts,period_days,N,Nt,GAUSS); 
		#PV.footprintComponents_1L(u_nd,v_nd,eta_nd,PV_prime,PV_BG,U0_nd,U,Umag,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,AmpF_nd,FORCE,r0,nu,BG,Fpos,ts,period_days,N,Nt,GAUSS); # For diagnostic purposes
		if doEEFs:
			EEFq_total, EEFq_north, EEFq_south = PV.equivEddyFlux(P_xav,y_nd,y0_nd,dy_nd,omega_nd,N);
			EEFq = np.array([EEFq_total,EEFq_north,EEFq_south]);
			print(EEFq_total, EEFq_north, EEFq_south);

			
# Buoyancy footprints
#====================================================

# Should these be zero, according to conservation of mass?
#Pb, Pb_xav = buoy.footprint(u_full,v_nd,eta_full,U0_nd,U,Umag,x_nd,y_nd,T_nd,dx_nd,dy_nd,dt_nd,AmpF_nd,FORCE,r0,nu,BG,Fpos,ts,period_days,N,Nt,GAUSS);

#====================================================

# Last step: normalise all solutions by the (non-dimensional) forcing amplitude. 
#u_nd = u_nd / AmpF_nd;
#v_nd = v_nd / AmpF_nd;
#eta_nd = eta_nd / AmpF_nd;
#if doPV:
#	PV_prime = PV_prime / AmpF_nd;

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
	diagnostics.forcingPlots(x_nd,y_nd,F1_nd,F2_nd,F3_nd,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N);
	#forcing_1L.forcingInv(Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,x_nd,y_nd,dx_nd,N); # For diagnostic purposes

# Background state plots (inc. BG SSH, BG flow, BG PV)
if plotBG:
	diagnostics.bgPlots(y_nd,H0_nd,U0_nd,PV_BG);

# Soltuion Plots
if plotSol:
	diagnostics.solutionPlots(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);
	#diagnostics.solutionPlotsDim(x,y,u,v,eta,ts,L,FORCE,BG,Fpos,N);

# Plots of PV and zonally averaged PV
if plotPV:
	pvPlots(PV_full,PV_prime,P,x_nd,y_nd);
if plotPV_av:
	diagnostics.PV_avPlots(x_nd,y_nd,PV_prime,PV_BG,PV_full,ts,FORCE,BG,Fpos,N);

# Plots of footprints - may need to edit the source code at times.
if plotFootprint:
	diagnostics.footprintPlots(x_nd,y_nd,P,P_xav,Fpos,BG,GAUSS,FORCE,nu,r0,period_days,U0_nd,U,N);

# Phase and amplitude
if plotPhaseAmp:
	diagnostics.solutionPlotsAmp(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);
	diagnostics.solutionPlotsPhase(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N);


