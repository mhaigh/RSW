# RESTART.py
#=======================================================

# This code works in the same way as the master file RSW_1L.py, but does not solve the SW model.
# Instead it takes saved solution (in .npy or .nc format), and then can perform all necessary diagnostics.
# Should ensure that the input file is defined in the same way as it was for originally producing the saved npy files.

#=======================================================

import sys

import numpy as np

from core import diagnostics, PV, forcing, solver, energy
from output import plotting_bulk

from inputFile import *

#=======================================================

test = 'sigma' # U0 or sigma

#samples = ['16']
#samples = ['-08','00','08','16']
samples = ['0'] 				# sigma samples

samples2 = [0]

ns = len(samples)
ns2 = len(samples2)

# Initialise u,v,eta,PV,P
u = np.zeros((N,N,ns),dtype=complex);
v = np.zeros((N,N,ns),dtype=complex);
h = np.zeros((N,N,ns),dtype=complex);
q = np.zeros((N,N,ns));
P = np.zeros((N,N+1,ns));
P_xav = np.zeros((N,ns))

#=======================================================

# First define all necessary terms to be plotted.
for si in range(0,ns):
	
	sample = samples[si]

	# For each U0 or sigma value, we need to first load the saved solutions
	if test == 'U0':
		path = '/media/mike/Seagate Expansion Drive/Documents/GulfStream/RSW/DATA/1L/PAPER1/UNIFORM/'
		U0_load = str(sample)
		u_tmp = np.load(path + 'u_U0='+U0_load+'.npy')[:,:,ts]
		v_tmp = np.load(path + 'v_U0='+U0_load+'.npy')[:,:,ts]
		h_tmp = np.load(path + 'eta_U0='+U0_load+'.npy')[:,:,ts]
		P[:,:,si] = np.load(path + 'P_U0='+U0_load+'.npy')
		P_xav[:,si] = np.trapz(P[:,:,si],x_nd,dx_nd,axis=1);

		# Last step: redefine U0 and H0 for each sample
		U0, H0 = BG_state.BG_uniform(float(sample)/100.,Hflat,f0,beta,g,y,N);
		U0_nd = U0 / U
		H0_nd = H0 / chi	
	else:
		path = '/media/mike/Seagate Expansion Drive/Documents/GulfStream/RSW/DATA/1L/PAPER1/GAUSSIAN/'

		u_tmp = np.load(path + 'u_y0='+sample+'.npy')[:,:,ts]
		v_tmp = np.load(path + 'v_y0='+sample+'.npy')[:,:,ts]
		h_tmp = np.load(path + 'eta_y0='+sample+'.npy')[:,:,ts]
		P[:,:,si] = np.load(path + 'P_y0='+sample+'.npy')
		P_xav[:,si] = np.trapz(P[:,:,si],x_nd,dx_nd,axis=1);


	# Calculate full flows.
	h_full = np.zeros((N,N))
	u_full = np.zeros((N,N))
	for j in range(0,N):
		h_full[j,:] = h_tmp[j,:] + H0_nd[j]
		u_full[j,:] = u_tmp[j,:] + U0_nd[j]

	# Snapshot of PV
	#q[:,:,si] = PV.PV_instant(u_tmp,v_tmp,h_tmp,u_full,h_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd,Ro)

	# Now we have snapshot of solution, snapshot of PV, and the footprint.
	u[:,:,si] = u_tmp;	v[:,:,si] = v_tmp;	h[:,:,si] = h_tmp;
	
print('Plotting...')

#=======================================================

# Second, create all plots.

# Solutions
if True:
	fig, axes = plt.subplots(nrows=ns,ncols=3,figsize=(22,7*ns))

	for si in range(0,ns):
		sample = samples[si]
		if test == 'U0':
			string = r'$U_{0} = ' + str(float(sample)/100) + '$'	
		else:
			if sample == '0':
				string = r'$y_{0}=0$'
			elif sample == '-1':
				string = r'$y_{0}=-\sigma$'
			else:
				string = r'$y_{0}=' + sample + '\sigma$'
		U0_str = 'U0 = ' + str(U0)
		plotting_bulk.plotSolutions(u[:,:,si],v[:,:,si],h[:,:,si],N,x_grid,y_grid,si,ns,string)
		plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=1.0);
		plt.savefig('fig0.png');

# Phase & Amp
if False:
	fig, axes = plt.subplots(nrows=ns,ncols=3,figsize=(22,2*7*ns2))

	for si in range(0,ns2):
		sample = samples[samples2[si]]
		if test == 'U0':
			string = r'$U_{0} = ' + str(float(sample)/100) + '$'
		U0_str = 'U0 = ' + str(U0)
		plotting_bulk.plotSolutionsAmpPhase(u[:,:,si],v[:,:,si],h[:,:,si],N,x_grid,y_grid,si,ns2,string,fig)
		plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.8);
		plt.savefig('fig1.png');


# Footprints
if False:
	fig, axes = plt.subplots(nrows=ns,ncols=3,figsize=(22,7*ns))

	for si in range(0,ns):
		sample = samples[si]
		if test == 'U0':
			string = r'$U_{0} = ' + str(float(sample)/100) + '$'
		U0_str = 'U0 = ' + str(U0)
		plotting_bulk.fp_PV_plot(q[:,:,si],P[:,:,si],P_xav[:,si],N,x_grid,y_grid,y_nd,si,ns,string)
		plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=1.0);
		plt.savefig('fig2.png');





			



