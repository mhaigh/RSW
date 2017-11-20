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

path = '/home/mike/Documents/GulfStream/RSW/DATA/1L/PAPER1/UNIFORM/'
U0_name = 'U0=16';
U0_str = r'$U_{0}=0.16$';

u_nd = np.load(path+'u_'+U0_name+'.npy');
v_nd = np.load(path+'v_'+U0_name+'.npy');
eta_nd = np.load(path+'eta_'+U0_name+'.npy');
P = np.load(path+'P_'+U0_name+'.npy');

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


ss = 0;
if ss == 1:
	u_nd = u_nd*0.01;
	v_nd = v_nd*0.01;
	eta_nd = eta_nd*0.01;
	P = P * 1.0e-8;
	np.save(path+'u_U0='+U0_name+'.npy',u_nd);
	np.save(path+'v_U0='+U0_name+'.npy',v_nd);
	np.save(path+'eta_U0='+U0_name+'.npy',eta_nd);
	np.save(path+'P_U0='+U0_name+'.npy',P);

#sys.exit();
P_xav = np.trapz(P,x_nd,dx_nd,axis=1);

#=======================================================

#plotting.solutionPlots_save(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,U0_str,x_grid,y_grid,True);
#plotting.footprintPlots_save(P,P_xav,x_nd,y_nd,ts,FORCE,BG,Fpos,N,U0_str,x_grid,y_grid);

EEF_array = PV.EEF(P_xav,y_nd,y0_nd,y0_index,dy_nd,N);
EEF_north = EEF_array[0]; EEF_south = EEF_array[1];
print(EEF_north, EEF_south);
print(EEF_north-EEF_south);




#=======================================================


