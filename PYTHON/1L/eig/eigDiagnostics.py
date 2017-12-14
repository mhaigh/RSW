# eigDiagnostics.py
# File containing functions to be called by the master script EIG.py.
#====================================================

import sys

import numpy as np
import matplotlib.pyplot as plt

from diagnostics import diff
from diagnostics import extend

import plotly.graph_objs as go
import plotly.plotly as py

#====================================================

# eigPlot
def eigPlots(u_proj,v_proj,eta_proj,u_nd,v_nd,eta_nd,x_nd,y_nd,x_grid,y_grid,sol):
	
	ulim = np.max(abs(u_nd));
	vlim = np.max(abs(v_nd));
	etalim = np.max(abs(eta_nd));

	if sol:

		plt.figure(1,figsize=[21,10]);

		plt.subplot(231);
		plt.pcolor(x_grid,y_grid,u_proj, cmap='bwr', vmin=-ulim, vmax=ulim);
		plt.text(0.3,0.45,'u proj',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.subplot(232);
		plt.pcolor(x_grid,y_grid,v_proj, cmap='bwr', vmin=-vlim, vmax=vlim);
		plt.text(0.3,0.45,'v proj',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.subplot(233);
		plt.pcolor(x_grid,y_grid,eta_proj, cmap='bwr', vmin=-etalim, vmax=etalim);
		plt.text(0.3,0.45,'eta proj',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();
	
		plt.subplot(234);
		plt.pcolor(x_grid,y_grid,u_nd, cmap='bwr', vmin=-ulim, vmax=ulim);
		plt.text(0.3,0.45,'u',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.subplot(235);
		plt.pcolor(x_grid,y_grid,v_nd, cmap='bwr', vmin=-vlim, vmax=vlim);
		plt.text(0.3,0.45,'v',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.subplot(236);
		plt.pcolor(x_grid,y_grid,eta_nd, cmap='bwr', vmin=-etalim, vmax=etalim);
		plt.text(0.3,0.45,'eta',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.tight_layout();
		plt.show();

	# Code never executed, unless we want contourf plots.
	elif 1==0:

		u_nd = extend(u_nd);
		v_nd = extend(v_nd);
		eta_nd = extend(eta_nd);

		plt.figure(1,figsize=[21,6]);
		plt.subplot(231)
		plt.contourf(x_nd,y_nd,u_proj,vmin=-ulim,vmax=ulim);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.clim(-ulim,ulim);
		plt.colorbar();

		plt.subplot(232)
		plt.contourf(x_nd,y_nd,v_proj,vmin=-vlim,vmax=vlim);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.xlabel('x');
		plt.ylabel('y');
		plt.clim(-vlim,vlim);
		plt.colorbar();

		plt.subplot(233)
		plt.contourf(x_nd,y_nd,eta_proj,vmin=-etalim,vmax=etalim);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.clim(-etalim,etalim);
		plt.colorbar();
	
		plt.subplot(234);
		plt.contourf(x_nd,y_nd,u_nd,vmin=-ulim,vmax=ulim);
		plt.text(0.3,0.45,'u',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-ulim,ulim);
		plt.colorbar();

		plt.subplot(235);
		plt.contourf(x_nd,y_nd,v_nd,vmin=-vlim,vmax=vlim);
		plt.text(0.3,0.45,'v',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-vlim,vlim);
		plt.colorbar();

		plt.subplot(236);
		plt.contourf(x_nd,y_nd,eta_nd,vmin=-etalim,vmax=etalim);
		plt.text(0.3,0.45,'eta',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-etalim,etalim);
		plt.colorbar();

		plt.tight_layout();
		plt.show();


	else:	

		plt.figure(1,figsize=[21,6]);
		plt.subplot(131)
		plt.contourf(x_nd,y_nd,u_proj);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.subplot(132)
		plt.contourf(x_nd,y_nd,v_proj);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.subplot(133)
		plt.contourf(x_nd,y_nd,eta_proj);
		plt.xticks((-1./2,0,1./2));
		plt.yticks((-1./2,0,1./2));	
		plt.xlabel('x');
		plt.ylabel('y');
		plt.colorbar();
		plt.tight_layout();
		plt.show();

#====================================================

# scatterWeight
def scatterWeight(k,l,theta,theta_abs_tot,dom_index,Nm,Nk_neg,Nk_pos,Fpos):
	
	Nk = Nk_neg + Nk_pos + 1;

	dim = Nk * Nm;

	# Produce a set of weights normalised by theta_abs_tot, which is either
	# the total weight of all wavenumbers or just the first Nm wavenumbers.
	theta_normalised = np.zeros(dim);
	for i in range(0,Nk):
		for j in range(0,Nm):
			theta_normalised[i*Nm+j] = np.abs(theta[j,i] / theta_abs_tot[i]);

	theta_max = max(theta_normalised);
	
	# Calculate the performance of the decompisition at each wavenumber; 1 is perfect.
	perf = np.zeros(Nk);
	for i in range(0,Nk):
		perf[i] = sum(theta_normalised[i*Nm:(i+1)*Nm]);

	# Some points may be repeated, the next few lines allow repeated points to be plotted as different shapes.
	# First arrange all points as list of tuples
	L = [];
	L_theta = [];
	for ii in range(0,dim):
		L.append(tuple([k[ii],l[ii]]));
		L_theta.append(tuple([k[ii],l[ii],theta_normalised[ii]]));

	# Turn this list into a set and count how many time each element in the set S appears in the list L.	
	S = set(L);
	S_theta = set(L_theta);
	F = {};
	for ii in list(S):	# Loop through each element in the set S
		F[ii] = L.count(ii);	# And count how many times it appears in L

	# We also need a dictionary containing the weights of each tuple (or maximum weight if a tuple is repeated).
 	F_theta = {};
	for wi in range(0,dim):
		ii = L[wi];
		if ii in F_theta:							# Check if this key is already in the dictionary.
			if theta_normalised[wi] > F_theta[ii]:	# If so, check if new theta is larger than the previous...
				F[ii] = theta_normalised[wi];		# and if so, overwrite the old theta value.
		else:
			F_theta[ii] = theta_normalised[wi];		# Create a new key, with corresponding value theta

	# Now loop through all the tuples (keys) in the dictionary, storing each tuple in an array depending on the number of times it is repeated.
	k1 = []; l1 = []; theta1 = [];
	k2 = []; l2 = []; theta2 = [];
	k3 = []; l3 = []; theta3 = [];
	k4 = []; l4 = []; theta4 = [];
	for ii in list(S):
		if F[ii] == 1:
			k1.append(ii[0]); l1.append(ii[1]); theta1.append(F_theta[ii]);
		elif F[ii] == 2:
			k2.append(ii[0]); l2.append(ii[1]); theta2.append(F_theta[ii]);
		elif F[ii] == 3:
			k3.append(ii[0]); l3.append(ii[1]); theta3.append(F_theta[ii]);
		else:
			k4.append(ii[0]); l4.append(ii[1]); theta4.append(F_theta[ii]);
		
	k1 = np.array(k1); l1 = np.array(l1); theta1 = np.array(theta1);
	k2 = np.array(k2); l2 = np.array(l2); theta2 = np.array(theta2);
	k3 = np.array(k3); l3 = np.array(l3); theta3 = np.array(theta3);
	k4 = np.array(k4); l4 = np.array(l4); theta4 = np.array(theta4);

	colors1 = theta1; colors2 = theta2;
	colors3 = theta3; colors4 = theta4;
	
	l_max = max(l);
	if l_max % 2 == 0:
		y_max = l_max + 2;
	else:
		y_max = l_max + 1;

	#print(colors2);
	
	y_ticks = np.linspace(0,y_max,y_max/2+1);

	cm = 'YlOrRd';

	plt.scatter(k1,l1,c=colors1,cmap=cm,vmin=0,vmax=theta_max,marker='s',s=50);
	plt.scatter(k2,l2,c=colors2,cmap=cm,vmin=0,vmax=theta_max,marker='o',s=50);
	plt.scatter(k3,l3,c=colors3,cmap=cm,vmin=0,vmax=theta_max,marker='^',s=50);
	plt.scatter(k4,l4,c=colors4,cmap=cm,vmin=0,vmax=theta_max,marker='*',s=50);
	plt.ylim([0,y_max+4]);	
	plt.yticks(y_ticks);
	plt.colorbar();#ticks=np.linspace(0,theta_max,10)
	plt.grid();
	for i in range(-Nk_neg,Nk_pos+1):
		plt.text(i-0.25,y_max+0.8,str(round(perf[i],2)),fontsize=14,rotation=90);
	plt.title('Eigenmode decomposition - ' + str(Fpos),fontsize=22);
	plt.text(Nk_pos-4,y_max-6,'Weight',fontsize=22);
	plt.xlabel('k',fontsize=22);
	plt.ylabel('l',fontsize=22);
	plt.show();	
	

#====================================================

# scatterPeriod
def scatterPeriod(k,l,p,dom_index,Nm,Nk_neg,Nk_pos,Fpos):
	
	Nk = Nk_neg + Nk_pos + 1;

	dim = Nk * Nm;

	p = abs(p);	# Interested in period regardless of direction.
	pmin = 0;
	pmax = 2*60.0;#np.max(p);

	# Some points may be repeated, the next few lines allow repeated points to be plotted as different shapes.
	# First arrange all points as list of tuples
	L = [];
	L_p = [];
	for ii in range(0,dim):
		L.append(tuple([k[ii],l[ii]]));
		L_p.append(tuple([k[ii],l[ii],p[ii]]));

	# Turn this list into a set and count how many time each element in the set S appears in the list L.	
	S = set(L);
	S_p = set(L_p);
	F = {};
	for ii in list(S):			# Loop through each element in the set S...
		F[ii] = L.count(ii);	# and count how many times it appears in L

	# We also need a dictionary containing the periods of each tuple.
 	F_p = {};
	for wi in range(0,dim):
		ii = L[wi];
		if ii in F_p:				# Check if this key is already in the dictionary.
			F[ii] = p[wi];		
		else:
			F_p[ii] = p[wi];		# Create a new key, with corresponding value p[wi].

	# Now loop through all the tuples (keys) in the dictionary, storing each tuple in an array depending on the number of times it is repeated.
	k1 = []; l1 = []; p1 = [];
	k2 = []; l2 = []; p2 = [];
	k3 = []; l3 = []; p3 = [];
	k4 = []; l4 = []; p4 = [];
	for ii in list(S):
		if F[ii] == 1:
			k1.append(ii[0]); l1.append(ii[1]); p1.append(F_p[ii]);
		elif F[ii] == 2:
			k2.append(ii[0]); l2.append(ii[1]); p2.append(F_p[ii]);
		elif F[ii] == 3:
			k3.append(ii[0]); l3.append(ii[1]); p3.append(F_p[ii]);
		else:
			k4.append(ii[0]); l4.append(ii[1]); p4.append(F_p[ii]);
		
	k1 = np.array(k1); l1 = np.array(l1); p1 = np.array(p1);
	k2 = np.array(k2); l2 = np.array(l2); p2 = np.array(p2);
	k3 = np.array(k3); l3 = np.array(l3); p3 = np.array(p3);
	k4 = np.array(k4); l4 = np.array(l4); p4 = np.array(p4);

	colors1 = p1; colors2 = p2;
	colors3 = p3; colors4 = p4;
	
	l_max = max(l);
	if l_max % 2 == 0:
		y_max = l_max + 2;
	else:
		y_max = l_max + 1;

	#print(colors2);
	
	y_ticks = np.linspace(0,y_max,y_max/2+1);

	#cm = 'Set1'
	cm = plt.cm.get_cmap('Set2',8);

	plt.scatter(k1,l1,c=colors1,cmap=cm,vmin=pmin,vmax=pmax,marker='s',s=50);
	plt.scatter(k2,l2,c=colors2,cmap=cm,vmin=pmin,vmax=pmax,marker='o',s=50);
	plt.scatter(k3,l3,c=colors3,cmap=cm,vmin=pmin,vmax=pmax,marker='^',s=50);
	plt.scatter(k4,l4,c=colors4,cmap=cm,vmin=pmin,vmax=pmax,marker='*',s=50);
	plt.ylim([0,y_max+4]);	
	plt.yticks(y_ticks);
	plt.colorbar();#ticks=np.linspace(0,theta_max,10)
	plt.grid();
	plt.title('Eigenmode decomposition - ' + str(Fpos),fontsize=22);
	plt.text(Nk_pos-6,y_max,'Period (days)',fontsize=22);
	plt.xlabel('k',fontsize=22);
	plt.ylabel('l',fontsize=22);
	plt.show();	

#====================================================

# orderEigenmodes
def orderEigenmodes(vec,val,N,VECS):
# A function that takes the set of eigenmodes, given by vec, and orders them according to the number of zero crossings.
# When two or more eigenmodes cross zeros the same amount of times, they are ordered by their frequency, smallest first.
	
	dim = np.size(val);

	if VECS:
		vec = np.array(vec);
		u_vec = vec[0,:,:];
		v_vec = vec[1,:,:];
		eta_vec = vec[2,:,:];

	else:
		u_vec = vec[0:N,:];		# Extract the eigenmodes.
		v_vec = vec[N:2*N,:];
		eta_vec = vec[2*N:3*N,:];		

	# Initialise a counter for the number of zero crossings. 
	count = np.zeros((dim),dtype=int);
	for wi in range(0,dim):
		for j in range(1,N):
			if (u_vec[j-1,wi] >= 0 and u_vec[j,wi] <= 0) or (u_vec[j-1,wi] <= 0 and u_vec[j,wi] >= 0):
				count[wi] = count[wi] + 1;

	i_count = np.argsort(count);

	return count, i_count;

#====================================================

# orderEigenmodes2
def orderEigenmodes2(vec,val,N,VECS):
# A function that takes the set of eigenmodes, given by vec, and orders them according to the number of zero crossings.
# When two or more eigenmodes cross zeros the same amount of times, they are ordered by their frequency, smallest first.
	
	dim = np.size(val);

	if VECS:
		vec = np.array(vec);
		u_vec = vec[0,:,:];
		v_vec = vec[1,:,:];
		eta_vec = vec[2,:,:];

	else:
		u_vec = vec[0:N,:];		# Extract the eigenmodes.
		v_vec = vec[N:2*N,:];
		eta_vec = vec[2*N:3*N,:];		

	# Initialise a counter for the number of zero crossings. 
	count = np.zeros((dim),dtype=int);
	for wi in range(0,dim):
		u_abs = np.abs(u_vec[:,wi]);
		u_av = np.mean(u_abs);
		for j in range(1,N):
			if (u_abs[j-1] >= u_av and u_abs[j] <= u_av) or (u_abs[j-1] <= u_av and u_abs[j] >= u_av):
				count[wi] = count[wi] + 1;
	
	for wi in range(0,dim):
		count[wi] = int(np.floor(count[wi] / 2));

	i_count = np.argsort(count);

	return count, i_count;

#====================================================

# updateCount
def updateCount(count,wii):
# A function that takes raw input to manually update the zero-crossings count of an eigenvector.
# Type end to quit updating.

	update_count = raw_input('-->');		# The first step updates count, but i_count will no longer match.									
	if update_count != '' and update_count != 'end':	
		print('updated');
		count_new = int(update_count);	# Stores the new count, to be used to rearrange the vectors.
 	elif update_count == 'end':
		count_new = count;
		wii = wii + dim; 	# End the loop, and don't update any more modes.	
	else:
		count_new = count;

	return count_new, wii;

#====================================================

# vec2vecs
def vec2vecs(vec,N,dim,BC):
# A function to take the full eigenvector, vec = (u,v,eta), and 
# extract u, v, eta, depending on the boundary condition.

	if BC == 'FREE-SLIP':
		u_vec = vec[0:N,:];
		v_vec = np.zeros((N,dim),dtype=complex);
		v_vec[1:N-1,:] = vec[N:2*N-2,:];
		eta_vec = vec[2*N-2:3*N-2,:];
	elif BC == 'NO-SLIP':
		u_vec = np.zeros((N,dim),dtype=complex);
		u_vec[1:N-1,:] = vec[0:N-2,:];
		v_vec = np.zeros((N,dim),dtype=complex);
		v_vec[1:N-1,:] = vec[N:2*N-2,:];
		eta_vec = vec[2*N-4:3*N-4,:];	
		return val, u_vec, v_vec, eta_vec;		
	else:
		sys.exit('ERROR: choose valid BC');

	return u_vec, v_vec, eta_vec;

#===================================================

	
	