# initialise.py

#=======================================================
#=======================================================
# Module consisting of intiialisation functions.
# Parameters are collated into appropriate classes.

import numpy as np

from inputFile import *

#=======================================================
#=======================================================

class Grid(object):
	# Grid
	def __init__(self,Lx,Ly,N,Hflat):
		self.Lx = Lx; 
		self.Ly = Ly;
		self.N = N;
		self.N2 = N-2;
		self.x = np.linspace(-Lx/2,Lx/2+Lx/(N-1),N+1);	# Array of all grid points in physical space.
		self.y = np.linspace(-Ly/2,Ly/2,N);				# Zonal grid points	set has an extra point, first and last points are duplicates of each other (periodicity).
		self.dx = self.x[1] - self.x[0];
		self.dy = self.y[1] - self.y[0];
		self.y_grid, self.x_grid = np.mgrid[slice(-Ly/2,Ly/2+dy,dy),slice(-Lx/2,Lx/2+dx,dx)];

class Grid_nd(object):
	# Grid
	def __init__(self,Lx,Ly,N,Hflat):
		self.Lx = Lx; 
		self.Ly = Ly;
		self.N = N;
		self.N2 = N-2;
		self.x = np.linspace(-Lx/2,Lx/2+Lx/(N-1),N+1);	# Array of all grid points in physical space.
		self.y = np.linspace(-Ly/2,Ly/2,N);				# Zonal grid points	set has an extra point, first and last points are duplicates of each other (periodicity).
		self.dx = self.x[1] - self.x[0];
		self.dy = self.y[1] - self.y[0];
		self.y_grid, self.x_grid = np.mgrid[slice(-Ly/2,Ly/2+dy,dy),slice(-Lx/2,Lx/2+dx,dx)];
