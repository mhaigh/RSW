# output.py
#=======================================================

# This file contains modules to be called by RSW_1L.py.
# It saves data in netcdf format using the netCDF4 module. 

#=======================================================

import numpy as np
import netCDF4 as nc

#=======================================================

# ncSave
def ncSave(utilde_nd,vtilde_nd,etatilde_nd,u_nd,v_nd,eta_nd,x_nd,y_nd,K_nd,T_nd,PV_FULL,PV_PRIME,PV_bg,Pq,EEFq,N,Nt):
# A function that saves the output of RSW_1L.py in netcdf format.
# Saves the solutions (physical and spectral by default) and depending on whether or not
# they were calculated, also saves the PV, footprints, and EEF.
# The last thing to be saved are the BG flow U0, forcing radius r0, kinematic viscosity nu, forcing period.

	# Initialise the nc file
	RSW1L = nc.Dataset('RSW1L.nc','w',format='NETCDF4');
		
	# Create dimensions
	x_dim = RSW1L.createDimension('x_dim',N+1);	
	y_dim = RSW1L.createDimension('y_dim',N);
	k_dim = RSW1L.createDimension('k_dim',N);
	t_dim = RSW1L.createDimension('t_dim',Nt);

	# Initialise dimension variables...
	x = RSW1L.createVariable('x','f8',('x_dim',));
	y = RSW1L.createVariable('y','f8',('y_dim',));
	k = RSW1L.createVariable('k','f8',('k_dim',));
	t = RSW1L.createVariable('t','f8',('t_dim',));
	# ...and assign the data.
	x[:] = x_nd;
	y[:] = y_nd;
	t[:] = T_nd[0:Nt];

	# Initialise solution variables...
	u = RSW1L.createVariable('u','f8',('y_dim','x_dim','t_dim',));
	v = RSW1L.createVariable('v','f8',('y_dim','x_dim','t_dim',));
	eta = RSW1L.createVariable('eta','f8',('y_dim','x_dim','t_dim',));
	utilde_real = RSW1L.createVariable('utilde_real','f4',('k_dim','y_dim',));
	vtilde_real = RSW1L.createVariable('vtilde_real','f4',('k_dim','y_dim',));
	etatilde_real = RSW1L.createVariable('etatilde_real','f4',('k_dim','y_dim',));	
	utilde_imag = RSW1L.createVariable('utilde_imag','f4',('k_dim','y_dim',));
	vtilde_imag = RSW1L.createVariable('vtilde_imag','f4',('k_dim','y_dim',));
	etatilde_imag = RSW1L.createVariable('etatilde_imag','f4',('k_dim','y_dim',));
	
	# ...and assign the data to the variables
	u[:,:,:] = u_nd;
	v[:,:,:] = v_nd;
	eta[:,:,:] = eta_nd;
	utilde_real[:,:] = np.real(utilde_nd);  
	vtilde_real[:,:] = np.real(vtilde_nd);
	etatilde_real[:,:] = np.real(etatilde_nd);
	utilde_imag[:,:] = np.imag(utilde_nd);  
	vtilde_imag[:,:] = np.imag(vtilde_nd);
	etatilde_imag[:,:] = np.imag(etatilde_nd);

	# Some variables (PV, footprint, EEF) are conditional on their existence in RSW_1L.py
	# Here we initialise them, and assign data to them.
	if PV_FULL != None:
		PV_full = RSW1L.createVariable('PV','f4',('t_dim','x_dim','y_dim',));
		PV_prime = RSW1L.createVariable('PV_prime','f4',('t_dim','x_dim','y_dim',));
		PV_BG = RSW1L.createVariable('PV_BG','f4',('y_dim',));
		#==
		PV_full[:,:,:] = PV_FULL;
		PV_prime[:,:,:] = PV_prime;
		PV_BG[:] = PV_bg;
	if Pq != None:
		P = RSW1L.createVariable('P','f4',('x_dim','y_dim',));
		#==
		P[:,:] = Pq;		
	if EEFq != None:
		EEF = RSW1L.createVariable('EEF','f4');
		#==
		EEF[:] = EEFq;

	RSW1L.close();

#=======================================================

# ncSaveEigmodes
def ncSaveEigmodes(u_modes,v_modes,eta_modes,val,y_nd,k,N,dim):
# A function that saves the output of RSW_1L.py in netcdf format.
# Saves the solutions (physical and spectral by default) and depending on whether or not
# they were calculated, also saves the PV, footprints, and EEF.
# The last thing to be saved are the BG flow U0, forcing radius r0, kinematic viscosity nu, forcing period.

	file_name = 'RSW1L_Eigmodes_k' + str(k) + '.nc';

	# Initialise the nc file
	RSW1L_Eigmodes = nc.Dataset(file_name,'w',format='NETCDF4');
		
	# Create dimensions
	y_dim = RSW1L_Eigmodes.createDimension('y_dim',N);
	omega_dim = RSW1L_Eigmodes.createDimension('omega_dim',dim)

	# Initialise dimension variables...
	y = RSW1L_Eigmodes.createVariable('y','f8',('y_dim',));
	omega = RSW1L_Eigmodes.createVariable('omega','f8',('omega_dim',));
	# ...and assign the data.
	y[:] = y_nd;
	omega[:] = val;
	
	# Initialise solution variables...
	u_vec = RSW1L.createVariable('u_vec','f8',('y_dim','omega_dim',));
	v_vec = RSW1L.createVariable('v_vec','f8',('y_dim','omega_dim',));
	eta_vec = RSW1L.createVariable('eta_vec','f8',('y_dim','omega_dim',));
	
	# ...and assign the data to the variables
	u_vec[:,:] = u_modes;
	v_vec[:,:] = v_modes;
	eta_vec[:,:] = eta_modes;

	RSW1L_Eigmodes.close();

