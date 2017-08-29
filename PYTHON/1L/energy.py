# energy
#=========================================================

# A set of energy-related functions

#=========================================================

import numpy as np
import matplotlib.pyplot as plt

#=========================================================

# energy_BG
def energy_BG(U0_nd,H0_nd,Ro,y_nd,dy_nd,N):
# A function that calculates the total energy E = KE + PE of the background state - 
# a useful quantity for when calculating the forcing induced energy or energy of eigenmodes.

	# First kinetic energy, KE = 0.5 * U**2 * H (dimensionless, no v term)
	KE_BG = 0.5 * U0_nd**2 * H0_nd;					# KE in terms of y (uniform in x)
	KE_BG_tot = np.trapz(KE_BG,y_nd,dy_nd);		# Total KE in the 2D domain (x-length of domain is 1)
	
	# Next the potential energy, PE = 0.5 * H**2 / Ro
	PE_BG = 0.5 * H0_nd**2 / Ro;				# PE in terms of y
	PE_BG_tot = np.trapz(PE_BG,y_nd,dy_nd);		# Total PE in the 2D domain

	return KE_BG, KE_BG_tot, PE_BG, PE_BG_tot;

#=========================================================

# KE_from_spec
def KE_from_spec(u_tilde,v_tilde,eta_tilde,k_nd,x_nd,y_nd,Nt,N,output):
# A function that takes the spectral representation of the solution at only one wavenumber k, indexed by i,
# and calculates the kinetic energy (KE) at that wavenumber.
# Because of the linearity of the system, the KE outputted from this function can be summed over all wavenumbers
# to produce the total KE.
# All values are dimensionless - see documentation for scaling arguments.

# Options for output:
# 1. 'av': outputs only the KE temporal average (as a function of x and y);
# 2. 'av_tot': outputs the temporal and spatial average;
# 2. 'full': outputs the full KE, as a function of x, y and t.
 
	I = np.complex(0,1);
	
	# We want the KE at a discrete set of times over one forcing/solution/eigenmode period.
	# Note that this calculation can be made independent of this period/frequency;
	# instead we only need to sample 'times' from a interval of unit length, with Nt entries,
	# the frequency omega cancels out with the period T.

	u = np.zeros((N,N,Nt));
	v = np.zeros((N,N,Nt));
	eta = np.zeros((N,N,Nt));

	omega_t = np.linspace(0,1,Nt);

	for ti in range(0,Nt):
		for j in range(0,N):
			u[j,:,ti] = np.real(u_tilde[j] * np.exp(2 * np.pi * (k_nd * x_nd[0:N] - omega_t[ti])));
			v[j,:,ti] = np.real(v_tilde[j] * np.exp(2 * np.pi * (k_nd * x_nd[0:N] - omega_t[ti])));
			eta[j,:,ti] = np.real(eta_tilde[j] * np.exp(2 * np.pi * (k_nd * x_nd[0:N] - omega_t[ti])));
	
	# 1. Temporally averaged KE
	if output == 'av':
		KE_av = 0.5 * (u[:,:,0]**2 + v[:,:,0]**2) * eta[:,:,0];
		for ti in range(1,Nt):
			KE_av = KE_av + 0.5 * (u[:,:,ti]**2 + v[:,:,ti]**2) * eta[:,:,ti];
		KE_av = KE_av / Nt;
		return KE_av;

	# 2. Temporal and spatial average of KE
	elif output == 'av_tot':
		dx_nd = x_nd[1] - x_nd[0];
		dy_nd = y_nd[1] - y_nd[0];
		KE = 0.5 * (u[:,:,0]**2 + v[:,:,0]**2) * eta[:,:,0];
		KE_av_tot = np.trapz(np.trapz(KE,x_nd[0:N],dx_nd,1),y_nd,dy_nd,0);
		for ti in range(1,Nt):
			KE = 0.5 * (u[:,:,ti]**2 + v[:,:,ti]**2) * eta[:,:,ti];
			KE_av_tot = KE_av_tot + np.trapz(np.trapz(KE,x_nd[0:N],dx_nd,1),y_nd,dy_nd,0);
		KE_av_tot = KE_av_tot / Nt;
		return KE_av_tot;
	
	# 3. Time-dependent KE
	elif output == 'full':
		KE = np.zeros((N,N,Nt));
		for ti in range(0,Nt):
			KE[:,:,ti] = 0.5 * (u[:,:,ti]**2 + v[:,:,ti]**2) * eta[:,:,ti];
		return KE;

	else:
		import sys
		sys.exit('Invalid output selection; must be "av", "av_tot" or "full".');

#=========================================================

# PE_from_spec
def PE_from_spec(eta_tilde,Ro,k_nd,x_nd,y_nd,Nt,N,output):
# A function that takes the spectral representation of the solution at only one wavenumber k, indexed by i,
# and calculates the potential energy (PE) at that wavenumber.
# Because of the linearity of the system, the PE outputted from this function can be summed over all wavenumbers
# to produce the total PE.
# All values are dimensionless - see documentation for scaling arguments.

# Options for output:
# 1. 'av': outputs only the PE temporal average (as a function of x and y);
# 2. 'av_tot': outputs the temporal and spatial average;
# 2. 'full': outputs the full PE, as a function of x, y and t.
 
	I = np.complex(0,1);
	
	# We want the KE at a discrete set of times over one forcing/solution/eigenmode period.
	# Note that this calculation can be made independent of this period/frequency;
	# instead we only need to sample 'times' from a interval of unit length, with Nt entries,
	# the frequency omega cancels out with the period T.

	eta = np.zeros((N,N,Nt));

	omega_t = np.linspace(0,1,Nt);

	for ti in range(0,Nt):
		for j in range(0,N):
			eta[j,:,ti] = np.real(eta_tilde[j] * np.exp(2 * np.pi * (k_nd * x_nd[0:N] - omega_t[ti])));
	
	# 1. Temporally averaged PE
	if output == 'av':
		PE_av = 0.5 * eta[:,:,0] / Ro;
		for ti in range(1,Nt):
			PE_av = PE_av + 0.5 * (u[:,:,ti]**2 + v[:,:,ti]**2) * eta[:,:,ti];
		PE_av = PE_av / Nt;
		return PE_av;

		# 1. Temporal and spatial average of PE
	elif output == 'av_tot':
		dx_nd = x_nd[1] - x_nd[0];
		dy_nd = y_nd[1] - y_nd[0];
		PE = 0.5 * eta[:,:,0]**2 / Ro;
		PE_av_tot = np.trapz(np.trapz(PE,x_nd[0:N],dx_nd,1),y_nd,dy_nd,0);
		for ti in range(1,Nt):
			PE = 0.5 * eta[:,:,ti]**2 / Ro;
			PE_av_tot = PE_av_tot + np.trapz(np.trapz(PE,x_nd[0:N],dx_nd,1),y_nd,dy_nd,0);
		PE_av_tot = PE_av_tot / Nt;
		return PE_av_tot;
	
	# 3. Time-dependent PE
	elif output == 'full':
		PE = np.zeros((N,N,Nt));
		for ti in range(0,Nt):
			PE[:,:,ti] = 0.5 * eta[:,:,ti]**2 / Ro;
		return PE;

	else:
		import sys
		sys.exit('Invalid output selection; must be "av", "av_tot" or "full".');

#=========================================================

# E_from_spec
def E_from_spec(u_tilde,v_tilde,eta_tilde,Ro,k_nd,x_nd,y_nd,Nt,N,output):
# A function that takes the spectral representation of the solution at only one wavenumber k, indexed by i,
# and calculates the energy (both KE and PE) at that wavenumber.
# Because of the linearity of the system, the energy outputted from this function can be summed over all wavenumbers
# to produce the energy.
# All values are dimensionless - see documentation for scaling arguments.

# Options for output:
# 1. 'av': outputs only the temporal average (as a function of x and y);
# 2. 'av_tot': outputs the temporal and spatial average;
# 2. 'full': outputs the full energy, as a function of x, y and t.
 
	I = np.complex(0,1);
	
	# We want the KE at a discrete set of times over one forcing/solution/eigenmode period.
	# Note that this calculation can be made independent of this period/frequency;
	# instead we only need to sample 'times' from a interval of unit length, with Nt entries,
	# the frequency omega cancels out with the period T.

	u = np.zeros((N,N,Nt));
	v = np.zeros((N,N,Nt));
	eta = np.zeros((N,N,Nt));
	
	omega_t = np.linspace(0,1,Nt);

	for ti in range(0,Nt):
		for j in range(0,N):
			u[j,:,ti] = np.real(u_tilde[j] * np.exp(2 * np.pi * I * (k_nd * x_nd[0:N] - omega_t[ti])));
			v[j,:,ti] = np.real(v_tilde[j] * np.exp(2 * np.pi * I *(k_nd * x_nd[0:N] - omega_t[ti])));
			eta[j,:,ti] = np.real(eta_tilde[j] * np.exp(2 * np.pi * I * (k_nd * x_nd[0:N] - omega_t[ti])));
	
	#plt.contourf(u[:,:,10]);
	#plt.colorbar();
	#plt.show();

	# 1. Temporally averaged KE and PE
	if output == 'av':
		KE_av = 0.5 * (u[:,:,0]**2 + v[:,:,0]**2) * eta[:,:,0];
		PE_av = 0.5 * eta[:,:,0] / Ro;
		for ti in range(1,Nt):
			KE_av = KE_av + 0.5 * (u[:,:,ti]**2 + v[:,:,ti]**2) * eta[:,:,ti];
			PE_av = PE_av + 0.5 * (u[:,:,ti]**2 + v[:,:,ti]**2) * eta[:,:,ti];
		KE_av = KE_av / Nt;
		PE_av = PE_av / Nt;
		return KE_av, PE_av;

		# 1. Temporal and spatial average of KE and PE
	elif output == 'av_tot':
		dx_nd = x_nd[1] - x_nd[0];
		dy_nd = y_nd[1] - y_nd[0];
		KE = 0.5 * (u[:,:,0]**2 + v[:,:,0]**2) * eta[:,:,0];
		KE_av_tot = np.trapz(np.trapz(KE,x_nd[0:N],dx_nd,1),y_nd,dy_nd,0);
		PE = 0.5 * eta[:,:,0]**2 / Ro;
		PE_av_tot = np.trapz(np.trapz(PE,x_nd[0:N],dx_nd,1),y_nd,dy_nd,0);
		for ti in range(1,Nt):
			KE = 0.5 * (u[:,:,ti]**2 + v[:,:,ti]**2) * eta[:,:,ti];
			KE_av_tot = KE_av_tot + np.trapz(np.trapz(KE,x_nd[0:N],dx_nd,1),y_nd,dy_nd,0);
			PE = 0.5 * eta[:,:,ti]**2 / Ro;
			PE_av_tot = PE_av_tot + np.trapz(np.trapz(PE,x_nd[0:N],dx_nd,1),y_nd,dy_nd,0);
		KE_av_tot = KE_av_tot / Nt;
		PE_av_tot = PE_av_tot / Nt;
		return KE_av_tot, PE_av_tot;
	
	# 3. Time-dependent KE and PE
	elif output == 'full':
		KE = np.zeros((N,N,Nt));
		PE = np.zeros((N,N,Nt));
		for ti in range(0,Nt):
			KE[:,:,ti] = 0.5 * (u[:,:,ti]**2 + v[:,:,ti]**2) * eta[:,:,ti];
			PE[:,:,ti] = 0.5 * eta[:,:,ti]**2 / Ro;
		return KE, PE;

	else:
		import sys
		sys.exit('Invalid output selection; must be "av", "av_tot" or "full".');

#=========================================================

# KE
def KE(u_full,v_nd,eta_full,x_nd,y_nd,dx_nd,dy_nd,N):
# Outputs the KE at a moment in time, by taking as input as time-snapshot of the solution u,v,eta.

	KE = 0.5 * (u_full**2 + v_nd**2) * eta_full;

	KE_tot = np.trapz(np.trapz(KE,x_nd[0:N],dx_nd,axis=1),y_nd,dy_nd);

	return KE, KE_tot;

#=========================================================

# PE
def PE(eta_full,Ro,x_nd,y_nd,dx_nd,dy_nd,N):
# Outputs the KE at a moment in time, by taking as input as time-snapshot of the solution u,v,eta.

	PE = 0.5 * eta_full**2 / Ro;

	PE_tot = np.trapz(np.trapz(PE,x_nd[0:N],dx_nd,axis=1),y_nd,dy_nd);

	return PE, PE_tot;




	
	



def ENER(b):
	u_nd = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/' + str(FORCE) + '/' + str(BG) + '/u_nd_' + str(Fpos) + str(N) + '.npy');
	v_nd = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/' + str(FORCE) + '/' + str(BG) + '/v_nd_' + str(Fpos) + str(N) + '.npy');
	eta_nd = np.load('/home/mike/Documents/GulfStream/Code/DATA/1L/' + str(FORCE) + '/' + str(BG) + '/eta_nd_' + str(Fpos) + str(N) + '.npy');
	
	I = np.complex(0,1);		# Define I = sqrt(-1)
	
	u_nd = extend(u_nd);
	v_nd = extend(v_nd);
	eta_nd = extend(eta_nd);
	
	# In order to calculate the energy forced into the system, we require the energy of the full system
	etaFull = np.zeros((N,N+1,Nt));
	uFull = np.zeros((N,N+1,Nt));
	for j in range(0,N):
		etaFull[j,:,:] = eta_nd[j,:,:] + H0_nd[j];
		uFull[j,:,:] = u_nd[j,:,:] + U0_nd[j];

	# Initialise the energy arrays (Full = perturbation + BG state)
	KEphysFull = np.zeros((N,N+1,Nt));
	KEtotFull = np.zeros(Nt);
	PEphysFull = np.zeros((N,N+1,Nt));
	PEtotFull = np.zeros(Nt);
	# Now calculate the KE and PE of the full system (i.e. including BG state) at each time.
	for ti in range(0,Nt):
		KEphysFull[:,:,ti] = 0.5 * (etaFull[:,:,ti]) * (uFull[:,:,ti]**2 + v_nd[:,:,ti]**2);	# v_nd = vFull
		KEtotFull[ti] = np.trapz(np.trapz(KEphysFull[:,:,ti],x,dx,1),y,dy,0);
		PEphysFull[:,:,ti] = 0.5 * g * etaFull[:,:,ti]**2;										# To keep it nd, there is no multiplicatio  by g
		PEtotFull[ti] = np.trapz(np.trapz(PEphysFull[:,:,ti],x,dx,1),y,dy,0);

	# Calculate the energies of the steady BG state
	KE_BG = np.zeros(N);
	PE_BG = np.zeros(N);
	KE_BG[0:N] = 0.5 * H0_nd[0:N] * U0_nd[0:N]**2;
	PE_BG[0:N] = 0.5 * g * H0_nd[0:N]**2;
	# Division by Umag to normalise the BG energies

	# Initialise the perturbation energy arrays
	KEphys = np.zeros((N,N+1,Nt));
	KEtot = np.zeros(Nt);
	PEphys = np.zeros((N,N+1,Nt));
	PEtot = np.zeros(Nt);
	for j in range(0,N):
		KEphys[j,:,:] = KEphysFull[j,:,:] - KE_BG[j];
		PEphys[j,:,:] = PEphysFull[j,:,:] - PE_BG[j];

	# A more direct alternative for calculating the energies - a useful check.
	# These use the faster algebraic definition given in the report and do not calculate energies using TOT-BG.
	#for j in range(0,N):
	#	PEphys[j,:,:] = 0.5 * g * eta_nd[j,:,:] * (eta_nd[j,:,:] + 2 * H0_nd[j]);
	#	KEphys[j,:,:] = 0.5 * etaFull[j,:,:] * (u_nd[j,:,:]**2 + 2 * u_nd[j,:,:] * U0_nd[j] + v_nd[j,:,:]**2) + 0.5 * eta_nd[j,:,:] * U0_nd[j]**2;

	# Normalise the energy by the forcing amplitude
	KEphys = KEphys / AmpF_nd**3;		
	PEphys = PEphys / AmpF_nd**2;

	KE_FullTot = np.zeros(Nt);
	PE_FullTot = np.zeros(Nt);

	# Calculates the extra PE and KE in the system at each time throughout the period of the forcing
	for ti in range(0,Nt):
		KEtot[ti] = np.trapz(np.trapz(KEphys[:,:,ti],x_nd,dx_nd,1),y_nd,dy_nd,0);
		PEtot[ti] = np.trapz(np.trapz(PEphys[:,:,ti],x_nd,dx_nd,1),y_nd,dy_nd,0);
		KE_FullTot[ti] = np.trapz(np.trapz(KEphysFull[:,:,ti],x_nd,dx_nd,1),y_nd,dy_nd,0);
		PE_FullTot[ti] = np.trapz(np.trapz(PEphysFull[:,:,ti],x_nd,dx_nd,1),y_nd,dy_nd,0);

	
	#=========================================================
	
	# PLOTS
	
	#=========================================================

	# Snapshot of PE and KE perturbation
	plt.figure(1)
	plt.subplot(121)
	plt.contourf(x_nd,y_nd,KEphys[:,:,ts]);
	plt.title('KE');
	plt.xticks((-1./2,0,1./2));
	plt.yticks((-1./2,0,1./2));
	plt.colorbar()
	plt.subplot(122)
	plt.contourf(x_nd,y_nd,PEphys[:,:,ts]);
	plt.title('PE');
	plt.xticks((-1./2,0,1./2));
	plt.yticks((-1./2,0,1./2));
	plt.colorbar()
	plt.show()
	
	# Time behaviour of KE and PE perturbations
	plt.figure(2)
	fig, ax1 = plt.subplots()
	ax1.plot(T_nd[:Nt], KEtot, 'b-',linewidth=2)
	ax1.set_xlabel('TIME',fontsize=18)
	# Make the y-axis label, ticks and tick labels match the line color.
	ax1.set_ylabel('KE', color='b',fontsize=18)
	ax2 = ax1.twinx()
	ax2.plot(T_nd[:Nt], PEtot, 'r-',linewidth=2)
	ax2.set_ylabel('PE', color='r',fontsize=18)
	#ax3 = ax1.twinx()
	#ax3.plot(T_nd,0.5*max(KEtot)*np.sin(T_nd),'k--',label='F amp',linewidth=2)
	fig.tight_layout()
	plt.savefig('/home/mike/Documents/GulfStream/Code/IMAGES/1L/' + str(FORCE) + '/' + str(BG) + '/energy_' + str(Fpos) + str(N) + '.png')
	plt.show()

	# Time behaviour of full KE and PE
	plt.figure(2)
	fig, ax1 = plt.subplots()
	ax1.plot(T_nd[:Nt], KE_FullTot, 'b-',linewidth=2)
	ax1.set_xlabel('TIME',fontsize=18)
	# Make the y-axis label, ticks and tick labels match the line color.
	ax1.set_ylabel('KE', color='b',fontsize=18)
	ax1.tick_params('y', colors='b')
	ax2 = ax1.twinx()
	ax2.plot(T_nd[:Nt], PE_FullTot, 'r-',linewidth=2)
	ax2.set_ylabel('PE', color='r',fontsize=18)
	ax2.tick_params('y', colors='r')
	#ax3 = ax1.twinx()
	#ax3.plot(T_nd,0.5*max(KEtot)*np.sin(T_nd),'k--',label='F amp',linewidth=2)
	fig.tight_layout()
	plt.show()



