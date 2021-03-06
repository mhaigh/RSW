# plotting.py
# File containing numerous plotting functions.
#====================================================

import numpy as np
import matplotlib.pyplot as plt

#====================================================

# plotForcing
def plotForcing(x_grid,y_grid,F1_nd,F2_nd,F3_nd,F6_nd):
# This function plots the four typically nonzero forcing terms in the 2L SW system.

	plt.subplot(221);
	plt.pcolor(x_grid,y_grid,F1_nd);
	plt.text(0.3,0.45,'F1',fontsize=22);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));	
	plt.ylabel('y',fontsize=16);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.colorbar();

	plt.subplot(222);
	plt.pcolor(x_grid,y_grid,F2_nd);
	plt.text(0.3,0.45,'F2',fontsize=22);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));	
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.colorbar();

	plt.subplot(223);
	plt.pcolor(x_grid,y_grid,F3_nd);
	plt.text(0.3,0.45,'F3',fontsize=22);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));	
	plt.xlabel('x',fontsize=16);
	plt.ylabel('y',fontsize=16);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.colorbar();

	plt.subplot(224);
	plt.pcolor(x_grid,y_grid,F6_nd);
	plt.text(0.3,0.45,'F1',fontsize=22);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));	
	plt.ylabel('y',fontsize=16);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.colorbar();

	plt.tight_layout();
	plt.show();

#====================================================

# solutionPlots
def solutionPlots(x,y,x_grid,y_grid,u,v,h,ts,N,contour):

	# Take absolute maximum.
	u1lim = np.max(abs(u[:,:,ts,0]));
	u2lim = np.max(abs(u[:,:,ts,1]));
	v1lim = np.max(abs(v[:,:,ts,0]));
	v2lim = np.max(abs(v[:,:,ts,1]));
	h1lim = np.max(abs(h[:,:,ts,0]));
	h2lim = np.max(abs(h[:,:,ts,1]));

	# Normalise each solution, extract snapshot.	
	u1 = u[:,:,ts,0] / u1lim;
	u2 = u[:,:,ts,1] / u2lim;
	v1 = v[:,:,ts,0] / v1lim;
	v2 = v[:,:,ts,1] / v2lim;
	h1 = h[:,:,ts,0] / h1lim;
	h2 = h[:,:,ts,1] / h2lim;	

	if contour:
		
		# Another good option is 'seismic' (slightly darker shades).
		plt.figure(1,figsize=(22.,13.));

		plt.subplot(231);
		plt.contourf(x[0:N],y,u1);
		plt.text(0.3,0.45,'u1',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');

		plt.subplot(232);
		plt.contourf(x[0:N],y,v1);
		plt.text(0.3,0.45,'v1',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');

		plt.subplot(233);
		plt.contourf(x[0:N],y,h1);
		plt.text(0.3,0.45,'h1',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.colorbar();

		plt.subplot(234);
		plt.contourf(x[0:N],y,u2);
		plt.text(0.3,0.45,'u2',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');

		plt.subplot(235);
		plt.contourf(x[0:N],y,v2);
		plt.text(0.3,0.45,'v1',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');

		plt.subplot(236);
		plt.contourf(x[0:N],y,h2);
		plt.text(0.3,0.45,'h1',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.colorbar();

		plt.tight_layout();
		#plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/' + str(Fpos) + '_'  + str(N) + '.png');
		plt.show();


	else:

		plt.figure(1,figsize=(22.,13.));

		# 1
		plt.subplot(231);
		plt.pcolor(x_grid, y_grid, u1, cmap='bwr', vmin=-1.,vmax=1.);
		plt.text(0.3,0.45,'u',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.subplot(232);
		plt.pcolor(x_grid, y_grid, v1, cmap='bwr', vmin=-1.,vmax=1.);
		plt.text(0.3,0.45,'v',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.subplot(233);
		plt.pcolor(x_grid, y_grid, h1, cmap='bwr', vmin=-1.,vmax=1.);
		plt.text(0.3,0.45,'eta',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		# 2
		plt.subplot(234);
		plt.pcolor(x_grid, y_grid, u2, cmap='bwr', vmin=-1.,vmax=1.);
		plt.text(0.3,0.45,'u',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.subplot(235);
		plt.pcolor(x_grid, y_grid, v2, cmap='bwr', vmin=-1.,vmax=1.);
		plt.text(0.3,0.45,'v',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.subplot(236);
		plt.pcolor(x_grid, y_grid, h2, cmap='bwr', vmin=-1.,vmax=1.);
		plt.text(0.3,0.45,'eta',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.tight_layout();
		#plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/' + str(Fpos) + '_'  + str(N) + '.png');
		plt.show();


		#cmap='coolwarm'

#====================================================

# solutionPlots_save
def solutionPlots_save(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,U0_str,x_grid,y_grid,div):
# Function that saves plots of the solutions (including PV) separately.

	ulim = np.max(abs(u_nd[:,:,ts]));
	vlim = np.max(abs(v_nd[:,:,ts]));
	etalim = np.max(abs(eta_nd[:,:,ts]));

	u_nd = u_nd / ulim;
	v_nd = v_nd / vlim;
	eta_nd = eta_nd / etalim

	u_nd = extend(u_nd);
	v_nd = extend(v_nd);
	eta_nd = extend(eta_nd);

	u_str = 'max=' + str(round(ulim,4));
	v_str = 'max=' + str(round(vlim,4));
	eta_str = 'max=' + str(round(etalim,4));	

	if div:
			
		plt.pcolor(x_grid, y_grid, u_nd[:,:,ts], cmap='bwr', vmin=-1., vmax=1.);
		plt.text(0.4,0.4,r'$u^{\prime}$',fontsize=26);
		plt.text(-0.45,0.4,U0_str,color='k',fontsize=22);
		plt.text(-0.45,-0.4,u_str,color='k',fontsize=18);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.tight_layout();
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/u_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();

		plt.pcolor(x_grid, y_grid, v_nd[:,:,ts], cmap='bwr', vmin=-1., vmax=1.);
		plt.text(0.4,0.4,r'$v^{\prime}$',fontsize=26);
		#plt.text(-0.45,0.4,U0_str,color='k',fontsize=22);
		plt.text(-0.45,-0.4,v_str,color='k',fontsize=18);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2),fontsize=0);	
		plt.xlabel('x',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.tight_layout();
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/v_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();

		plt.pcolor(x_grid, y_grid, eta_nd[:,:,ts], cmap='bwr', vmin=-1., vmax=1.);
		plt.text(0.4,0.4,r'$\eta^{\prime}$',fontsize=26);
		#plt.text(-0.45,0.4,U0_str,color='k',fontsize=22);
		plt.text(-0.45,-0.4,eta_str,color='k',fontsize=18);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2),fontsize=0);	
		plt.xlabel('x',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.colorbar();
		plt.tight_layout();
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/eta_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();



	else:

		plt.figure(1)
		plt.contourf(x_nd,y_nd,u_nd[:,:,ts]);
		plt.text(0.4,0.4,r'$u^{\prime}$',fontsize=26);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=18);
		plt.ylabel('y',fontsize=18);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-ulim,ulim);
		plt.colorbar();
		plt.tight_layout();
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/u_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();

		plt.figure(2)
		plt.contourf(x_nd,y_nd,v_nd[:,:,ts]);
		plt.text(0.4,0.4,r'$v^{\prime}$',fontsize=26);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.xlabel('x',fontsize=18);
		plt.ylabel('y',fontsize=18);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-vlim,vlim);
		plt.colorbar();
		plt.tight_layout();
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/v_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();

		plt.figure(3)
		plt.contourf(x_nd,y_nd,eta_nd[:,:,ts]);
		plt.text(0.4,0.4,r'$\eta^{\prime}$',fontsize=26);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.xlabel('x',fontsize=18);
		plt.ylabel('y',fontsize=18);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-etalim,etalim);
		plt.colorbar();
		plt.tight_layout();
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/eta_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();
		
#====================================================

# solutionPlotsAmp
# Plots of amplitude 
def solutionPlotsAmp(x_grid,y_grid,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,U0_name,U0_str,N):

	ulim = np.max(abs(u_nd));
	vlim = np.max(abs(v_nd));
	etalim = np.max(abs(eta_nd));

	u_nd = u_nd / ulim;
	v_nd = v_nd / vlim;
	eta_nd = eta_nd / etalim

	u_nd = extend(u_nd);
	v_nd = extend(v_nd);
	eta_nd = extend(eta_nd);

	plt.pcolor(x_grid, y_grid, np.absolute(u_nd), vmin=0., vmax=0.5);
	plt.text(0.4,0.4,r'$u^{\prime}$',fontsize=26,color='w');
	plt.text(-0.45,0.4,U0_str,fontsize=22,color='w');
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));	
	plt.xlabel('x',fontsize=16);
	plt.ylabel('y',fontsize=16);
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/u_Amp_' + U0_name + '.png');
	plt.close();
	
	plt.pcolor(x_grid, y_grid, np.absolute(v_nd), vmin=0., vmax=0.5);
	plt.text(0.4,0.4,r'$v^{\prime}$',fontsize=26,color='w');
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2),fontsize=0);	
	plt.xlabel('x',fontsize=16);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/v_Amp_' + U0_name + '.png');
	plt.close();

	plt.pcolor(x_grid, y_grid, np.absolute(eta_nd), vmin=0., vmax=0.5);
	plt.text(0.4,0.4,r'$\eta^{\prime}$',fontsize=26,color='w');
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2),fontsize=0);
	plt.xlabel('x',fontsize=16);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/eta_Amp_' + U0_name + '.png');
	plt.close();


#====================================================

# solutionPlotsPhase
# Plots of phase 
def solutionPlotsPhase(x_grid,y_grid,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,U0_name,U0_str,N):

	u_nd = extend(u_nd);
	v_nd = extend(v_nd);
	eta_nd = extend(eta_nd);

	plt.pcolor(x_grid, y_grid, np.angle(u_nd), cmap='hsv',vmin=-np.pi, vmax=np.pi);
	plt.text(0.4,-0.4,r'$u^{\prime}$',fontsize=26,color='k');
	plt.text(-0.45,-0.4,U0_str,fontsize=22,color='k');
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2),);	
	plt.xlabel('x',fontsize=16);
	plt.ylabel('y',fontsize=16);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/u_Phase_' + U0_name + '.png');
	plt.close();

	plt.pcolor(x_grid, y_grid, np.angle(v_nd),cmap='rainbow',vmin=-np.pi, vmax=np.pi);
	plt.text(0.4,-0.4,r'$v^{\prime}$',fontsize=26,color='k');
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2),fontsize=0);
	plt.xlabel('x',fontsize=16);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/v_Phase_' + U0_name + '.png');
	plt.close();

	fig, ax = plt.subplots()
	cax = ax.pcolor(x_grid, y_grid, np.angle(eta_nd),vmin=-np.pi,vmax=np.pi);
	cbar = fig.colorbar(cax, ticks=[-np.pi, 0, np.pi])
	cbar.ax.set_yticklabels([r'$-\pi$', r'$0$', r'$\pi$'],fontsize=18)  # vertically oriented colorbar
	plt.text(0.4,-0.4,r'$\eta^{\prime}$',fontsize=26,color='k');
	plt.xticks((-1./2,-1./4,0,1./4,1./2));	
	plt.yticks((-1./2,-1./4,0,1./4,1./2),fontsize=0);
	plt.xlabel('x',fontsize=16);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/eta_Phase_' + U0_name + '.png');
	plt.close();


#====================================================

# pvPlots
# Plots of PV and footprint
def pvPlots(PV_full,PV_prime,x_nd,y_nd):

	N = len(y_nd);

	plt.figure(1,figsize=[13,6]);

	plt.subplot(121);
	plt.contourf(x_nd[0:N],y_nd,PV_full[:,:,1]);
	plt.text(0.05,0.4,'PV FULL',color='k',fontsize=22);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();

	plt.subplot(122);
	plt.contourf(x_nd[0:N],y_nd,PV_prime[:,:,1]);
	plt.text(0.05,0.4,'PV PRIME',color='k',fontsize=22);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();

	plt.tight_layout()
	plt.show()

#====================================================

# pvPlots_save
def pvPlots_save(PV_full,PV_prime,x_nd,y_nd,ts,FORCE,BG,Fpos,N,U0_str,x_grid,y_grid,U0_name):

	PV_full_lim = np.max(abs(PV_full[:,:,ts]));
	PV_prime_lim = np.max(abs(PV_prime[:,:,ts]));

	PV_full = PV_full / PV_full_lim;
	PV_prime = PV_prime / PV_prime_lim;

	PV1_str = 'max=' + str(round(PV_full_lim,2));
	PV2_str = 'max=' + str(round(PV_prime_lim,2));

	plt.figure(1);
	plt.pcolor(x_grid, y_grid, PV_full[:,:,ts], cmap='bwr', vmin=-.5, vmax=.5);
	plt.text(0.4,0.4,r'$q$',color='k',fontsize=26);
	plt.text(-0.45,0.4,U0_str,color='k',fontsize=22);
	#plt.text(-0.45,-0.4,PV1_str,color='k',fontsize=18);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.xlabel('x',fontsize=18);
	plt.ylabel('y',fontsize=18);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/PV_full_' + str(Fpos) + '_' + U0_name + '.png');
	plt.close();

	plt.figure(2);
	plt.pcolor(x_grid, y_grid, PV_prime[:,:,ts], cmap='bwr', vmin=-.5, vmax=.5);
	plt.text(0.4,0.4,r'$q^{\prime}$',color='k',fontsize=26);
	plt.text(-0.45,0.4,U0_str,color='k',fontsize=22);
	#plt.text(-0.45,-0.4,PV2_str,color='k',fontsize=18);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.xlabel('x',fontsize=18);
	plt.ylabel('y',fontsize=18);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/PV_prime_' + str(Fpos) + '_'  + U0_name + '.png');
	plt.close();

#====================================================

# footprintPlots_save
def footprintPlots_save(P,P_xav,x_nd,y_nd,ts,FORCE,BG,Fpos,N,U0_str,x_grid,y_grid,U0_name):

	Plim = np.max(np.absolute(P));
	
	P = P / Plim;
	
	P_xav = P_xav * 1.0e5;

	P_str = 'max=' + str(round(Plim,2));

	plt.figure(1);
	plt.pcolor(x_grid, y_grid, P, cmap='bwr',vmin=-.5,vmax=.5);
	plt.text(0.4,0.4,r'$P$',fontsize=26);
	#plt.text(-0.45,-0.4,P_str,color='k',fontsize=18);
	#plt.text(0.25,0.4,str(Fpos),fontsize=18);		# Comment out this line if text on the plot isn't wanted.
	#plt.text(0.15,0.4,'r0 = '+str(r0/1000) + ' km' ,fontsize=18);	
	#plt.text(0.25,0.4,str(int(period_days))+' days',fontsize=18)
	#plt.text(0.25,0.4,'U0 = ' + str(U*U0_nd[0]),fontsize=18);
	#plt.text(0.25,0.4,r'$\nu$ = ' + str(int(nu)),fontsize=18);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2),fontsize=0);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.xlabel('x',fontsize=18);
	#plt.ylabel('y',fontsize=18);
	plt.colorbar();
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/P_' + str(Fpos) + '_' + U0_name + '.png');
	plt.close();
		
	plt.figure(2);
	plt.plot(P_xav,y_nd,'k-',linewidth=2)
	plt.text(0.8*max(abs(P_xav)),0.40,r'$\langle P\rangle$',fontsize=26)
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.ylim(-0.5,0.5);
	# Change these next two lines depending on which BG flow we're working with.
	#plt.xlim(-7.,7);	
	#plt.xticks((-6.,-4.,-2.,0.,2.,4.,6.));
	plt.xlim(-1.1*np.max(abs(P_xav)),1.1*np.max(abs(P_xav)));
	# ===
	plt.ylabel('y',fontsize=18);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.xlabel('x',color='white',fontsize=18);
	plt.tight_layout()
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/P_xav_' + str(Fpos) + '_' + U0_name + '.png');
	plt.close();
	

#====================================================

# forcingPlots
# Forcing plots 
def footprintPlots(x,y,P,P_xav):
# Function that plots the forcing, and its Fourier representation.

	Plim = np.max(abs(P));

	plt.figure(1,figsize=(15,7))
	plt.subplot(121)
	plt.contourf(x,y,P);
	plt.text(0.0,0.4,'PV FOOTPRINT',fontsize=22);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.xlabel('x',fontsize=16);
	plt.ylabel('y',fontsize=16);
	plt.clim(-Plim,Plim)
	plt.colorbar()

	plt.subplot(122)
	plt.plot(P_xav,y,linewidth=2)
	plt.text(10,0.40,'ZONAL AVERAGE',fontsize=22)
	plt.yticks((-1./2,0,1./2))
	plt.ylim(-0.5,0.5)
	plt.xlim(-1.1*np.max(abs(P_xav)),1.1*np.max(abs(P_xav)));
	plt.tight_layout()
	plt.show();
	


#====================================================

# footprintComponentsPlot
# A function that plots the footprint components.
def footprintComponentsPlot(uq,Uq,uQ,vq,vQ,P,P_uq,P_Uq,P_uQ,P_vq,P_vQ,P_xav,P_uq_xav,P_uQ_xav,P_Uq_xav,P_vq_xav,P_vQ_xav,x_nd,y_nd,N,Nt):

	uq_tav = np.zeros((N,N));
	Uq_tav = np.zeros((N,N));
	uQ_tav = np.zeros((N,N));
	vq_tav = np.zeros((N,N));
	vQ_tav = np.zeros((N,N));
	for j in range(0,N):
		for i in range(0,N):
			uq_tav[i,j] = sum(uq[i,j,:]) / Nt;
			Uq_tav[i,j] = sum(Uq[i,j,:]) / Nt;
			uQ_tav[i,j] = sum(uQ[i,j,:]) / Nt;
			vq_tav[i,j] = sum(vq[i,j,:]) / Nt;
			vQ_tav[i,j] = sum(vQ[i,j,:]) / Nt;

	print('Plotting time-averages of footprint components');

	plt.figure(1);

	plt.subplot(231);
	plt.contourf(uq_tav);
	plt.title('uq');
	plt.colorbar();

	plt.subplot(232);
	plt.contourf(uQ_tav);
	plt.title('uQ');
	plt.colorbar();

	plt.subplot(233);
	plt.contourf(Uq_tav);
	plt.title('Uq');
	plt.colorbar();

	plt.subplot(234);
	plt.contourf(vq_tav);
	plt.title('vq');
	plt.colorbar();

	plt.subplot(235);
	plt.contourf(vQ_tav);
	plt.title('vQ');
	plt.colorbar();

	plt.tight_layout();
	plt.show();

	#==

	print('Plotting contribution to footprint of each component (i.e. differentiated in space)');

	plt.figure(2);

	plt.subplot(231);
	plt.contourf(P_uq);
	plt.title('uq');
	plt.colorbar();

	plt.subplot(232);
	plt.contourf(P_uQ);
	plt.title('uQ');
	plt.colorbar();

	plt.subplot(233);
	plt.contourf(P_Uq);
	plt.title('Uq');
	plt.colorbar();

	plt.subplot(234);
	plt.contourf(P_vq);
	plt.title('vq');
	plt.colorbar();

	plt.subplot(235);
	plt.contourf(P_vQ);
	plt.title('vQ');
	plt.colorbar();

	plt.subplot(236);
	plt.contourf(P);
	plt.title('P');
	plt.colorbar();

	plt.tight_layout();
	plt.show();

	#==

	print('Plotting zonal averages');	
	
	plt.figure(3);

	plt.subplot(231);
	plt.plot(P_uq_xav,y_nd,linewidth=2);
	plt.title('uq');

	plt.subplot(232);
	plt.plot(P_uQ_xav,y_nd,linewidth=2);
	plt.title('uQ');

	plt.subplot(233);
	plt.plot(P_Uq_xav,y_nd,linewidth=2);
	plt.title('Uq');

	plt.subplot(234);
	plt.plot(P_vq_xav,y_nd,linewidth=2);
	plt.title('vq');

	plt.subplot(235);
	plt.plot(P_vQ_xav,y_nd,linewidth=2);
	plt.title('vQ');

	plt.subplot(235);
	plt.plot(P_xav,y_nd,linewidth=2);
	plt.title('P_xav');

	plt.tight_layout()
	plt.show();

#====================================================

# plotFluxes
def plotPrimaryComponents(P_uq,P_vq,P_uq_xav,P_vq_xav,x_nd,y_nd,FORCE,BG,Fpos,N):
# This function plots two footprint components uq and vq, and saves the output.
# These are the two fluxes that we are most interested in.

	print('Now plotting Primary components: contributions from uq and vq');
	
	P_uq_lim = np.max(abs(P_uq));
	P_vq_lim = np.max(abs(P_vq));

	U0_str = r'$U_{0}=0.16$';

	plt.figure(10);
	plt.contourf(x_nd,y_nd,P_uq);
	plt.text(0.4,0.4,r'$P_{u}$',color='k',fontsize=26);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.xlabel('x',fontsize=18);
	plt.ylabel('y',fontsize=18);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.clim(-P_uq_lim,P_uq_lim);
	plt.colorbar();
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/P_uq_' + str(Fpos) + '_'  + str(N) + '.png');
	plt.close();
	
	plt.figure(11);
	plt.contourf(x_nd,y_nd,P_vq);
	plt.text(0.4,0.4,r'$P_{v}$',color='k',fontsize=26);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.xlabel('x',fontsize=18);
	plt.ylabel('y',fontsize=18);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.clim(-P_vq_lim,P_vq_lim);
	plt.colorbar();
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/P_vq_' + str(Fpos) + '_'  + str(N) + '.png');
	plt.close();

	plt.figure(12);
	fig, ax1 = plt.subplots();
	ax1.plot(P_uq_xav,y_nd,'b-',linewidth=2);
	#ax1.set_xlabel(r'$\langle P_{v}\rangle$',color='b',fontsize=24);
	ax1.tick_params('x',colors='b');
	# Make the y-axis label, ticks and tick labels match the line color.
	ax1.set_yticks((-1./2,-1./4,0,1./4,1./2));
	ax1.set_ylabel('y',fontsize=18);
	ax1.text(-0.05,0.4,r'$\langle P_{v}\rangle$',color='r',fontsize=26);
	ax1.text(0.04,-0.4,r'$\langle P_{u}\rangle$',color='b',fontsize=26);
	#ax1.tick_params('y');
	ax2 = ax1.twiny();
	ax2.plot(P_vq_xav,y_nd,'r-',linewidth=2);
	#ax2.set_xlabel(r'$\langle P_{v}\rangle$', color='r',fontsize=24);
	ax2.tick_params('x',colors='r');
	#ax2.tick_params('y', colors='r');
	fig.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) + '/P_uv_' + str(Fpos) + str(N) + '.png')

#====================================================

# bgPlots
# Background state plots
def bgPlots(y_nd,H0_nd,U0_nd,PV_BG):

	plt.figure(2);
	plt.subplot(131);
	plt.plot(H0_nd,y_nd);
	plt.yticks((-1./2,0,1./2));
	plt.title('BG SSH');
	plt.subplot(132);
	plt.plot(U0_nd,y_nd);
	plt.yticks((-1./2,0,1./2));
	plt.title('BG flow U0');
	plt.subplot(133);
	plt.plot(PV_BG,y_nd);
	plt.yticks((-1./2,0,1./2));
	plt.title('BG PV');
	plt.show()

#====================================================

# plot_xshift
def plot_xshift():

	xshift=np.load('/home/mike/Documents/GulfStream/RSW/DATA/1L/Paper1/xshift.npy');	

	N = len(xshift);
	U0 = np.linspace(-0.5,0.5,N);

	plt.plot(U0,xshift,'k',linewidth=2);
	plt.ylim(-.12,.12);
	plt.xlim(-0.5,0.5);
	plt.xticks((-0.5,-.4,-.3,-.2,-.1,0.,0.1,0.2,0.3,0.4,0.5));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.xlabel('U0',fontsize=18)
	plt.ylabel('Footprint zonal shift',fontsize=18);
	plt.show();	


#====================================================

# forcingPlots
def forcingPlots(x_nd,y_nd,F1_nd,F2_nd,F3_nd,Ftilde1_nd,Ftilde2_nd,Ftilde3_nd,N):
# Function that plots the forcing, and its Fourier representation.

	plt.figure(1);

	plt.subplot(331);
	plt.contourf(x_nd,y_nd,F1_nd);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.subplot(332);
	plt.contourf(x_nd,y_nd,F2_nd);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.subplot(333);
	plt.contourf(x_nd,y_nd,F3_nd);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar()

	plt.subplot(334);
	plt.contourf(np.real(Ftilde1_nd));
	plt.colorbar()
	plt.subplot(335);
	plt.contourf(np.real(Ftilde2_nd));
	plt.colorbar()
	plt.subplot(336);
	plt.contourf(np.real(Ftilde3_nd));
	plt.colorbar()

	plt.subplot(337);
	plt.contourf(np.imag(Ftilde1_nd));
	plt.colorbar()
	plt.subplot(338);
	plt.contourf(np.imag(Ftilde2_nd));
	plt.colorbar()
	plt.subplot(339);
	plt.contourf(np.imag(Ftilde3_nd));
	plt.colorbar()

	plt.show();


#====================================================

# extend
def extend(f):
# A function used to replace the extra x-gridpoint on a solution.

	dimx = np.shape(f)[1];
	dimy = np.shape(f)[0];
	if f.size != dimx * dimy:
		dimt = np.shape(f)[2];

		f_new = np.zeros((dimy,dimx+1,dimt),dtype=f.dtype);
		for i in range(0,dimx):
			f_new[:,i,:] = f[:,i,:];
	
		f_new[:,dimx,:] = f[:,0,:];
	
	else:
		f_new = np.zeros((dimy,dimx+1),dtype=f.dtype);
		for i in range(0,dimx):
			f_new[:,i] = f[:,i];
	
		f_new[:,dimx] = f[:,0];

	return f_new


