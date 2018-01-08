# inputFile_2L
#=======================================================
#=======================================================
# File of input parameters for the 1L RSW plunger code

import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
#=======================================================
#=======================================================


# Domain
#=======================================================

BC = 'FREE-SLIP';			# Two boundary condition choices at north and south boundaries: NO-SLIP or FREE-SLIP.

N = 128; 			# Number of gridpoints in each direction.
	
Lx = 3840000.		# Zonal lengthscale (m)
Ly = 3840000.		# Meridional lengthscale (m)

H1flat = 1000.		# Depth of upper layer for ocean at rest (i.e. without BG flow SSH adjustment) (m)		
H2flat = 3000.		# Depth of lower layer for ocean at rest

f0 = 0.83e-4; # 0.0001214  		# Base value of Coriolis parameter (s-1)
beta = 2.0e-11;     		# Planetary vorticity gradient (m-1 s-1)

# Other physical parameters
#=======================================================

g = 9.81;		# Acceleration due to gravity (m s-2)
gamma = 4e-8;	# Frictional coefficient (s-1)
nu = 100.;		# Kinematic viscosity (m2 s-1)

rho1 = 1020.;    # Density of upper layer (kg m-3)
rho2 = 1030.;	# Density of lower layer (kg m-3) (should be greater than rho1)

# Background flow
#=======================================================

# Define a y-dependent zonal BG flow in each layer. From U1 and U2, the code prescribes the BG geostrophic interface heights
# so that the BG steady state satisfies the geostrophic equations.

# Upper layer
BG1 = 'NONE';			# Options: UNIFORM, GAUSSIAN, NONE.

Umag = 0.16; # 0.8 for Gaussian BG flow.

#sigma = 0.02 * 3840000.0;	# Jet width
#JET_POS = 'CENTER';

# Lower layer;
BG2 = 'NONE';			# ONLY NONE FOR NOW

# Forcing
#=======================================================

FORCE1 = 'BALANCED';       	# 'BALANCED' for geostrophically balanced forcing, 
							# 'VORTICITY' for forcing on the momentum eqautions only,
							# 'BUOYANCY' for forcing on continuity equation only 'USER'.
FORCE2 = 'NONE'				# NONE


Fpos = 'CENTER';			# 4 choices for positioning of plunger, 'NORTH', 'CENTER' and 'SOUTH'

# Instead of defining the forcing amplitude in the forcing module, we define it here as other codes require this value for normalisation
r0 = 90. * 1000.;	# (m)
AmpF = 1.0e-7;

# Time parameters
#=======================================================

period_days = 60.;						# Periodicity of plunger (days)
Nt = 200;								# Number of time samples
ts = Nt-1; 								# index at which the time-snapshot is taken

#======================================================================================================================================================================================================
#======================================================================================================================================================================================================

# Only edit from here to alter USER defined routines.
# Remainder initialises RSW.py

#=======================================================

# Initialisation

#=======================================================
#=======================================================

# Domain
#=======================================================

N2 = N-2;

L = Lx;

y = np.linspace(-Ly/2,Ly/2,N);				# Array of all grid points in physical space.
x = np.linspace(-Lx/2,Lx/2+Lx/(N-1),N+1);	# Zonal grid points	set has an extra point, first and last points are duplicates of each other (periodicity).

# Remove dead points in the y-direction
yd = np.zeros(N2);
for j in range(0,N2):
    yd[j] = y[j+1];

dy = y[1] - y[0];
dx = x[1] - x[0];

y_grid, x_grid = np.mgrid[slice(-Ly/2,Ly/2+dy,dy),slice(-Lx/2,Lx/2+2.*dx,dx)];

K = np.fft.fftfreq(N,Lx/N); 		 # Array of x-gridpoints in wavenumber space

f = f0 + beta * y;      # Coriolis frequency (s-1)
   
# Other physical parameters
#=======================================================

rho1_nd = rho1 / rho2;
rho2_nd = (rho2 - rho1) / rho2;

# Background flow
#=======================================================

if BG2 = 'NONE':
	if BG1 = 'UNIFORM':	
		U1, U2, H1, H2 = BG_state.BG_uniform_none(Umag1,Hflat1,Hflat2,f0,beta,g,y,N);
	elif BG1 == 'GAUSSIAN':
		U1, U2, H1, H2 = BG_state.BG_Gaussian_none(Umag1,sigma,JET_POS,Hflat1,Hflat2,f0,beta,g,y,Ly,N);
	elif BG1 = 'NONE:
		U1, U2, H1, H2 = BG_state.BG_none_none(Hflat1,Hflat2);	
	else:
		raise ValueError('Invalid BG flow option selected');

# Code not updated to include nonzero lower layer BG flows yet.

# Forcing
#=======================================================

if Fpos == 'NORTH':
	y0 = Ly / 4;
if Fpos == 'CENTER':
	y0 = 0;
if Fpos == 'SOUTH':
	y0 = - Ly / 4;
if Fpos == 'USER':
	y0 = y[N/2+3];

# Be careful here to make sure that the plunger is forcing boundary terms.

# Time parameters
#=======================================================

period = 3600. * 24. * period_days;		# Periodicity of plunger (s)
omega = 1. / (period);          		# Frequency of plunger, once every 50 days (e-6) (s-1)
T = np.linspace(0,period,Nt+1);			# Array of time samples across one forcing period (s)
dt = T[1] - T[0];						# Size of the timestep (s)
t = T[ts];								# Time of the snapshot

#======================================================================================================================================================================================================
#======================================================================================================================================================================================================

# Non-dimensionalisation
#=======================================================
#=======================================================
# Here we denote all non-dimensional parameters by '_nd'.

# To find the characteristic velocity U and depthscale H, we find the spatial-average of the 
# geostrophic BG state U0 and H0.
#=======================================================

U = 1.0;				# Typical U value
H = H1flat + H2flat;	# Motionless ocean depth

chi = f0 * U * Ly / g;		# The scaling for eta0 and eta1

T_adv = Lx / U;				# Advective timescale

# Defining dimensionless parameters
#=======================================================

Lx_nd = Lx / Ly;		# In case Lx and Ly are chosen to be different, we still scale by the same length, Ly
Ly_nd = Ly / Ly;				

y_nd = y / Ly;    
x_nd = x / Ly;
yd_nd = yd / Ly;

dy_nd = y_nd[1] - y_nd[0];
dx_nd = x_nd[1] - x_nd[0];

K_nd = K * Lx;		# The same as: K_nd = np.fft.fftfreq(N2,dx_nd)

H1_nd = H1 / chi;	# BG SSH1 scales the same way as eta0.
U1_nd = U1 / U;
H2_nd = H2 / chi;	
U2_nd = U2 / U;

f0_nd = 1.0;				# =f0/f0      		 
beta_nd = beta * Ly / f0;
f_nd = f / f0;			# The same as: f_nd = f0_nd + beta_nd * y_nd      

gamma_nd = gamma / f0			# Simply scaled by the base Coriolis frequency

omega_nd = omega * T_adv;      # All time variable scale advectively, i.e. T_adv~L/U
t_nd = t / T_adv;
T_nd = T / T_adv;
dt_nd = dt / T_adv;

AmpF_nd = AmpF * g / (f0 * U**2);

# Note that gravity g and kinematic viscosity aren't scaled, but rather used to define some extra dimensionless parameters

# Important dimensionless numbers
#=======================================================

Ro = U / (f0 * Ly); 			# Rossby number: measures inertial forces relative to rotational ones.
Re = Ly * U / nu;				# Reynolds number: measures inertial forces relative to viscous ones.
Ld1 = np.sqrt(g * H1flat) / f0;	# Rossby def radius.
Ld2 = np.sqrt(g* H2flat) / f0;	# Rossby def radius.

# ======================================================

# Print essential parameters to the terminal
print('SCALES: U = ' + str(U) + ', H = ' + str(H) + ', T = ' + str(T_adv));
print('FORCING PERIOD = ' + str(period_days) + ' DAYS');
print(str(FORCE1) + ' FORCING WITH ' + str(BG1) + ' BG FLOW')
print('Ro = ' + str(Ro));
print('Re = ' + str(Re));
print('Ld1 = ' + str(Ld1));
print('N = ' + str(N));

