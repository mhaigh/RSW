#PLOT.py
#=======================================================

# This module plots sols and footprint from saved data.

#=======================================================

import sys

import numpy as np
import diagnostics
import output
import plotting
import PV

from inputFile import *

#=======================================================

# UNIFORM


U0_name = 'U0=08';
U0_str = r'$U_{0}=0.08$';
#U0_name = 'y0=15sigma';
#U0_str = r'$y_{0}=1.5\sigma$';

#path = '/home/mike/Documents/GulfStream/RSW/DATA/1L/Paper1/';
path = '/media/mike/Seagate Expansion Drive/Documents/GulfStream/RSW/DATA/1L/PAPER1/UNIFORM/';
#path = '/media/mike/Seagate Expansion Drive/Documents/GulfStream/RSW/DATA/1L/PAPER1/GAUSSIAN/' + U0_name + '/';


u_nd = np.load(path+'u_'+U0_name+'.npy');
v_nd = np.load(path+'v_'+U0_name+'.npy');
eta_nd = np.load(path+'eta_'+U0_name+'.npy');
#P = np.load(path+'P_'+U0_name+'.npy');

#u_nd = np.load(path+'u_nd_complex_' + U0_name + '.npy');
#v_nd = np.load(path+'v_nd_complex_' + U0_name + '.npy');
#eta_nd = np.load(path+'eta_nd_complex_' + U0_name + '.npy');

#plotting.solutionPlotsPhase(x_grid,y_grid,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,U0_name,U0_str,N);
#plotting.solutionPlotsAmp(x_grid,y_grid,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,U0_name,U0_str,N);

#sys.exit();

#=======================================================

# GAUSSIAN

#U0_name = 'y0=0';
#U0_str = r'$y_{0}=-\sigma$';
#U0_str = r'$y_{0}=0$';

#path = '/home/mike/Documents/GulfStream/RSW/DATA/1L/PAPER1/GAUSSIAN/'
#u_nd = np.load(path + '/' + U0_name + '/u_' + U0_name + '.npy');
#v_nd = np.load(path + '/' + U0_name + '/v_' + U0_name + '.npy');
#eta_nd = np.load(path + '/' + U0_name + '/eta_' + U0_name + '.npy');
#P = np.load(path + '/' + U0_name + '/P_' + U0_name + '.npy');
#u_nd = np.load('/home/mike/Documents/GulfStream/RSW/PYTHON/1L/u_nd.npy');
#v_nd = np.load('/home/mike/Documents/GulfStream/RSW/PYTHON/1L/v_nd.npy');
#eta_nd = np.load('/home/mike/Documents/GulfStream/RSW/PYTHON/1L/eta_nd.npy');
#P = np.load('/home/mike/Documents/GulfStream/RSW/PYTHON/1L/P.npy');

#=======================================================


#ss = 0;
#if ss == 1:
#	u_nd = u_nd/0.01;
#	v_nd = v_nd/0.01;
#	eta_nd = eta_nd/0.01;
#	#P = P * 1.0e-8;
#	np.save(path+'u_U0='+U0_name+'.npy',u_nd);
#	np.save(path+'v_U0='+U0_name+'.npy',v_nd);
#	np.save(path+'eta_U0='+U0_name+'.npy',eta_nd);
#	#np.save(path+'P_U0='+U0_name+'.npy',P);

if False:
	eta_full = np.zeros((N,N,Nt));
	u_full = np.zeros((N,N,Nt));
	for j in range(0,N):
		eta_full[j,:,:] = eta_nd[j,:,:] + H0_nd[j];
		u_full[j,:,:] = u_nd[j,:,:] + U0_nd[j];
	PV_prime, PV_full, PV_BG = PV.potentialVorticity(u_nd,v_nd,eta_nd,u_full,eta_full,H0_nd,U0_nd,N,Nt,dx_nd,dy_nd,f_nd,Ro);
#P_xav = np.trapz(P,x_nd,dx_nd,axis=1);

#=======================================================

#plotting.pvPlots_save(PV_full,PV_prime,x_nd,y_nd,ts,FORCE,BG,Fpos,N,U0_str,x_grid,y_grid,U0_name);
plotting.solutionPlots_save(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,U0_str,x_grid,y_grid,True);
#plotting.footprintPlots_save(P,P_xav,x_nd,y_nd,ts,FORCE,BG,Fpos,N,U0_str,x_grid,y_grid,U0_name);


#EEF_array = PV.EEF(P_xav,y_nd,y0_nd,y0_index,dy_nd,N);
#EEF_north = EEF_array[0]; EEF_south = EEF_array[1];
#print(EEF_north, EEF_south);
#print(EEF_north-EEF_south);




#=======================================================


