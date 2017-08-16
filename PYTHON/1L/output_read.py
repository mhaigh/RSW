# output_read.py
#=======================================================

# This files read netcdf data, as produced by output.py. 

#=======================================================

import numpy as np
import netCDF4 as nc

#=======================================================

def ncReadEigenmodes(ncFile):

	I = np.complex(0,1);
	
	Eigenmodes = nc.Dataset(ncFile);	# Read the NETCDF file

	# Extract the eigenmodes and reconstruct from real and imag. parts.
	vec_tmp = Eigenmodes.variables['vec'][:,:,:];
	vec = vec_tmp[:,:,0] + I * vec_tmp[:,:,1];

	# Extract the eigenvalues
	val_tmp = Eigenmodes.variables['omega'][:,:];
	val = val_tmp[:,0] + I * val_tmp[:,1]; 

	# And extract the pseudo-wavenumber = count
	count = Eigenmodes.variables['count'][:];

	return val, vec, count


#=======================================================

def ncReadEEF_y0_components(ncFile):

	dataset = nc.Dataset(ncFile);	# Read the NETCDF file

	EEF = dataset.variables['EEF'][:,:];
	uq = dataset.variables['uq'][:,:];
	Uq = dataset.variables['Uq'][:,:];
	uQ = dataset.variables['uQ'][:,:];
	vq = dataset.variables['vq'][:,:];
	vQ = dataset.variables['vQ'][:,:];

	return EEF, uq, Uq,uQ, vq, vQ;
	

	
