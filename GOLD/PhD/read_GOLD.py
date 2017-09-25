# read_GOLD.py
#=======================================================

# This file reads GOLD ouput and plots it 

#=======================================================

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

#=======================================================

ncFile =  '/home/mike/cluster/gold/prog__0004_003.nc'
#ncFile =  '/home/mike/cluster/gold5/RESTART/GOLD.res_Y0074_D238_S00000.nc'

ncFile = nc.Dataset(ncFile);	# Read the NETCDF file

# Extract the eigenmodes and reconstruct from real and imag. parts.
PV = ncFile.variables['PV'][:,:,:,:];
print(np.shape(PV));
nt = np.shape(PV)[0];
PV = PV[nt-1,0,:,:];

print(np.shape(PV));
plt.contourf(PV);
plt.show();



