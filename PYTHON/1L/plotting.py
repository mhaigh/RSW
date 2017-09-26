# plotting.py
# File containing numerous plotting functions.
#====================================================

import numpy as np
import matplotlib.pyplot as plt
from diagnostics import extend

#====================================================

# forcingPlot_save
def forcingPlot_save(x_grid,y_grid,F3_nd,FORCE,BG,Fpos,N):

	aa = 1./24;
	Fmax = np.max(np.max(F3_nd,axis=0));
	Flim = np.max(abs(F3_nd/(1.05*Fmax)));
	#F3[0,0] = - Fmax;

	plt.pcolor(x_grid,y_grid,F3_nd/(1.1*Fmax),cmap='bwr',vmin=-Flim,vmax=Flim);
	plt.xlabel('x',fontsize=22);
	plt.ylabel('y',fontsize=22);
	plt.text(0.4,0.4,r'$F_{3}$',fontsize=26);
	plt.arrow(-aa,2*aa+0.25,2*aa,0,head_width=0.02, head_length=0.02,color='k');
	plt.arrow(2*aa,aa+.25,0,-2*aa,head_width=0.02, head_length=0.02,color='k');
	plt.arrow(aa,-2*aa+.25,-2*aa,0,head_width=0.02, head_length=0.02,color='k');
	plt.arrow(-2*aa,-aa+.25,0,2*aa,head_width=0.02, head_length=0.02,color='k');
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));	
	plt.xlabel('x',fontsize=16);
	plt.ylabel('y',fontsize=16);
	plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();
	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/F_' + str(Fpos) + '_'  + str(N) + '.png');
	plt.close();

# solutionPlots
def solutionPlots(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,div):

	ulim = np.max(abs(u_nd[:,:,ts]));
	vlim = np.max(abs(v_nd[:,:,ts]));
	etalim = np.max(abs(eta_nd[:,:,ts]));
	
	if div:
			
		plt.figure(1,figsize=(22,6.4));

		plt.subplot(131);
		plt.pcolor(x_grid, y_grid, u_nd[:,:,ts], cmap='bwr', vmin=-ulim, vmax=ulim);
		plt.text(0.3,0.45,'u',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.subplot(132);
		plt.pcolor(x_grid, y_grid, v_nd[:,:,ts], cmap='bwr', vmin=-vlim, vmax=vlim);
		plt.text(0.3,0.45,'v',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.colorbar();

		plt.subplot(133);
		plt.pcolor(x_grid, y_grid, eta_nd[:,:,ts], cmap='bwr', vmin=-etalim, vmax=etalim);
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

		# Another good option is 'seismic' (slightly darker shades).

	else:

		plt.figure(1,figsize=(22,6.4));

		plt.subplot(131);
		plt.contourf(x_nd[0:N],y_nd,u_nd[:,:,ts]);
		plt.text(0.3,0.45,'u',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));	
		plt.xlabel('x',fontsize=16);
		plt.ylabel('y',fontsize=16);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-ulim,ulim);
		plt.colorbar();

		plt.subplot(132);
		plt.contourf(x_nd[0:N],y_nd,v_nd[:,:,ts]);
		plt.text(0.3,0.45,'v',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-vlim,vlim);
		plt.colorbar();

		plt.subplot(133);
		plt.contourf(x_nd[0:N],y_nd,eta_nd[:,:,ts]);
		plt.text(0.3,0.45,'eta',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-etalim,etalim);
		plt.colorbar();

		plt.tight_layout();
		#plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/' + str(Fpos) + '_'  + str(N) + '.png');
		plt.show();

		#cmap='coolwarm'

#====================================================

# solutionPlots_save
def solutionPlots_save(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,div):
# Function that saves plots of the solutions (including PV) separately.

	ulim = np.max(abs(u_nd[:,:,ts]));
	vlim = np.max(abs(v_nd[:,:,ts]));
	etalim = np.max(abs(eta_nd[:,:,ts]));

	u_nd = u_nd / ulim;
	v_nd = v_nd / vlim;
	eta_nd = eta_nd / etalim

	U0_str = r'$U_{0}=-0.16$';

	if div:
			
		plt.pcolor(x_grid, y_grid, u_nd[:,:,ts], cmap='bwr', vmin=-1., vmax=1.);
		plt.text(0.4,0.4,r'$u^{\prime}$',fontsize=26);
		plt.text(-0.45,0.4,U0_str,color='k',fontsize=22);
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
		#plt.text(-0.45,0.4,U0_str,color='k',fontsize=22)
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
		#plt.text(-0.45,0.4,U0_str,color='k',fontsize=22)
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

# solutionPlotsDim
# Plots of the dimensional solutions
def solutionPlotsDim(x,y,u,v,eta,ts,L,FORCE,BG,Fpos,N):

	plt.figure(1,figsize=(22,6));

	plt.subplot(131);
	plt.contourf(x,y,u[:,:,ts]);
	plt.text(y[N-N/6],x[N-N/12],'u',fontsize=22);
	plt.xticks((-L/2,0,L/2),['-L/2','0','L/2']);
	plt.yticks((-L/2,0,L/2),['-L/2','0','L/2']);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();

	plt.subplot(132);
	plt.contourf(x,y,v[:,:,ts]);
	plt.text(y[N-N/6],x[N-N/12],'v',fontsize=22);
	plt.xticks((-L/2,0,L/2),['-L/2','0','L/2']);
	plt.yticks((-L/2,0,L/2),['-L/2','0','L/2']);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();

	plt.subplot(133);
	plt.contourf(x,y,eta[:,:,ts]);
	plt.text(y[N-N/6],x[N-N/12],'eta',fontsize=22);
	plt.xticks((-L/2,0,L/2),['-L/2','0','L/2']);
	plt.yticks((-L/2,0,L/2),['-L/2','0','L/2']);
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();

	plt.tight_layout();
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/' + str(Fpos) + '_'  + str(N) + '.png');
	plt.show();


#====================================================

# solutionPlotsAmp
# Plots of amplitude 
def solutionPlotsAmp(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N):

	plt.figure(1);

	plt.subplot(131);
	plt.contourf(x_nd,y_nd,np.absolute(u_nd[:,:,ts]));
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();

	plt.subplot(132);
	plt.contourf(x_nd,y_nd,np.absolute(v_nd[:,:,ts]));
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();

	plt.subplot(133);
	plt.contourf(x_nd,y_nd,np.absolute(eta_nd[:,:,ts]));
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.colorbar();

	plt.tight_layout()
	plt.show()

#====================================================

# solutionPlotsPhase
# Plots of phase 
def solutionPlotsPhase(x_nd,y_nd,u_nd,v_nd,eta_nd,ts,FORCE,BG,Fpos,N):

	plt.figure(1)

	plt.subplot(131);
	plt.contourf(x_nd,y_nd,np.angle(u_nd[:,:,ts]));
	plt.text(0.05,0.4,'u PHASE',color='w',fontsize=12);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');

	plt.subplot(132);
	plt.contourf(x_nd,y_nd,np.angle(v_nd[:,:,ts]));
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');

	plt.subplot(133);
	plt.contourf(x_nd,y_nd,np.angle(eta_nd[:,:,ts]));
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');

	plt.tight_layout()
	plt.show()

#====================================================

# PV_avPlots
# Plots of the zonal averages of PV, both the BG state and the forced state
def PV_avPlots(x_nd,y_nd,PV_prime,PV_BG,PV_full,ts,FORCE,BG,Fpos,N):


	PV_full_av = np.trapz(PV_full[:,:,ts],x_nd[0:N],x_nd[1]-x_nd[0],1);

	plt.figure(1)

	plt.plot(PV_full_av,y_nd);
	plt.plot(PV_BG,y_nd);

	plt.tight_layout()
	plt.show()


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
def pvPlots_save(PV_full,PV_prime,P,P_xav,x_nd,y_nd,ts,FORCE,BG,Fpos,N,x_grid,y_grid,div):

	PV_full_lim = np.max(abs(PV_full[:,:,ts]));
	PV_prime_lim = np.max(abs(PV_prime[:,:,ts]));
	Plim = np.max(abs(P));

	PV_full = PV_full / PV_full_lim;
	PV_prime = PV_prime / PV_prime_lim;
	P = P / Plim; 
	P_xav = P_xav / Plim;

	U0_str = r'$U_{0}=-0.16$';

	if div:

		plt.figure(1);
		plt.pcolor(x_grid, y_grid, PV_full[:,:,ts], cmap='bwr', vmin=-1., vmax=1.);
		plt.text(0.4,0.4,r'$q$',color='k',fontsize=26);
		plt.text(-0.45,0.4,U0_str,color='k',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.xlabel('x',fontsize=18);
		plt.ylabel('y',fontsize=18);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.tight_layout();
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/PV_full_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();

		plt.figure(1);
		plt.pcolor(x_grid, y_grid, PV_prime[:,:,ts], cmap='bwr', vmin=-1., vmax=1.);
		plt.text(0.4,0.4,r'$q^{\prime}$',color='k',fontsize=26);
		plt.text(-0.45,0.4,U0_str,color='k',fontsize=22);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.xlabel('x',fontsize=18);
		plt.ylabel('y',fontsize=18);
		plt.axis([x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()]);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.tight_layout();
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/PV_prime_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();
	
		plt.figure(1);
		plt.pcolor(x_grid, y_grid, P, cmap='bwr', vmin=-1., vmax=1.);
		plt.text(0.4,0.4,r'$P$',fontsize=26);
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
		plt.tight_layout()
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/P_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();
		
		plt.figure(4);
		plt.plot(P_xav,y_nd,'k-',linewidth=2)
		plt.text(0.8*max(abs(P_xav)),0.40,r'$\langle P\rangle$',fontsize=26)
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.ylim(-0.5,0.5);
		plt.ylabel('y',fontsize=18);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.xlim(-1.1*np.max(abs(P_xav)),1.1*np.max(abs(P_xav)));
		plt.xlabel('x',color='white',fontsize=18);
		plt.tight_layout()
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/P_xav_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();

	else:

		plt.figure(1);
		plt.contourf(x_nd,y_nd,PV_full[:,:,ts]);
		plt.text(0.4,0.4,r'$q$',color='k',fontsize=26);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.xlabel('x',fontsize=18);
		plt.ylabel('y',fontsize=18);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-PV_full_lim,PV_full_lim);
		plt.colorbar();
		plt.tight_layout();
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/PV_full_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();

		plt.figure(2);
		plt.contourf(x_nd,y_nd,PV_prime[:,:,ts]);
		plt.text(0.4,0.4,r'$q^{\prime}$',color='k',fontsize=26);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.xlabel('x',fontsize=18);
		plt.ylabel('y',fontsize=18);
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.clim(-PV_prime_lim,PV_prime_lim);
		plt.colorbar();
		plt.tight_layout();
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/PV_prime_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();
	
		plt.figure(3);
		plt.contourf(x_nd,y_nd,P);
		plt.text(0.4,0.4,r'$P$',fontsize=26);
		#plt.text(0.25,0.4,str(Fpos),fontsize=18);		# Comment out this line if text on the plot isn't wanted.
		#plt.text(0.15,0.4,'r0 = '+str(r0/1000) + ' km' ,fontsize=18);	
		#plt.text(0.25,0.4,str(int(period_days))+' days',fontsize=18)
		#plt.text(0.25,0.4,'U0 = ' + str(U*U0_nd[0]),fontsize=18);
		#plt.text(0.25,0.4,r'$\nu$ = ' + str(int(nu)),fontsize=18);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.xlabel('x',fontsize=18);
		plt.ylabel('y',fontsize=18);
		plt.clim(-Plim,Plim);
		plt.colorbar();
		plt.tight_layout()
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/P_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();
		
		plt.figure(4);
		plt.plot(P_xav,y_nd,linewidth=2)
		#plt.text(10,0.40,r'$\langle P\rangle$',fontsize=26)
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.ylim(-0.5,0.5);
		plt.ylabel('y',fontsize=18);
		plt.xlim(-1.1*np.max(abs(P_xav)),1.1*np.max(abs(P_xav)));
		plt.xlabel(r'$\langle P\rangle$',fontsize=26);
		plt.tight_layout()
		plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/P_xav_' + str(Fpos) + '_'  + str(N) + '.png');
		plt.close();

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

# forcingPlots
# Forcing plots 
def footprintPlots(x_nd,y_nd,P,P_xav,Fpos,BG,GAUSS,FORCE,nu,r0,period_days,U0_nd,U,N):
# Function that plots the forcing, and its Fourier representation.

	Plim = np.max(abs(P));

	plt.figure(1,figsize=(15,7))
	plt.subplot(121)
	#plt.contourf(x_nd,y_nd,P,cmap='coolwarm')
	plt.contourf(x_nd[0:N],y_nd,P)
	plt.text(0.0,0.4,'PV FOOTPRINT',fontsize=22);
	#plt.text(0.25,0.4,str(Fpos),fontsize=18);		# Comment out this line if text on the plot isn't wanted.
	#plt.text(0.15,0.4,'r0 = '+str(r0/1000) + ' km' ,fontsize=18);	
	#plt.text(0.25,0.4,str(int(period_days))+' days',fontsize=18)
	#plt.text(0.25,0.4,'U0 = ' + str(U*U0_nd[0]),fontsize=18);
	#plt.text(0.25,0.4,r'$\nu$ = ' + str(int(nu)),fontsize=18);
	plt.xticks((-1./2,-1./4,0,1./4,1./2));
	plt.yticks((-1./2,-1./4,0,1./4,1./2));
	plt.grid(b=True, which='both', color='0.65',linestyle='--');
	plt.xlabel('x',fontsize=16);
	plt.ylabel('y',fontsize=16);
	plt.clim(-Plim,Plim)
	plt.colorbar()
	plt.subplot(122)
	plt.plot(P_xav,y_nd,linewidth=2)
	plt.text(10,0.40,'ZONAL AVERAGE',fontsize=22)
	plt.yticks((-1./2,0,1./2))
	plt.ylim(-0.5,0.5)
	plt.xlim(-1.1*np.max(abs(P_xav)),1.1*np.max(abs(P_xav)));
	plt.tight_layout()
	plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/FOOTPRINT_nu=' + str(nu) + '.png');
	plt.show();
	
	# These if loops are for constantly altering depending on the test being done.
	if BG == 'GAUSSIAN':
		plt.figure(2)
		plt.contourf(x_nd,y_nd,P)
		plt.plot(U*U0_nd/(1.5*Umag)-0.5,y_nd,'k--',linewidth=2);
		plt.text(0.25,0.4,str(GAUSS),fontsize=18);
		#plt.text(0.25,0.4,str(period_days)+' days',fontsize=18)
		#plt.text(0.25,0.4,str(Fpos),fontsize=18);
		#plt.plot(P_xav[:,ts],y_nd,linewidth=2)
		#plt.text(0.25,0.4,'r0 = '+str(r0/1000),fontsize=18);	
		plt.colorbar()
		plt.ylim([-0.5,0.5]);
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.xlabel('x');
		plt.ylabel('y');
		plt.tight_layout()
		#plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/TEST/FOOTPRINT_' + str(GAUSS) + '.png');
		
	if BG == 'UNIFORM':
		plt.figure(2)
		plt.contourf(x_nd,y_nd,P)
		plt.text(0.25,0.4,'U0 = ' + str(U*U0_nd[0]),fontsize=18);
		plt.colorbar()
		plt.xticks((-1./2,-1./4,0,1./4,1./2));
		plt.yticks((-1./2,-1./4,0,1./4,1./2));
		plt.grid(b=True, which='both', color='0.65',linestyle='--');
		plt.tight_layout()
		#plt.savefig('/home/mike/Documents/GulfStream/RSW/IMAGES/1L/' + str(FORCE) + '/' + str(BG) +  '/FOOTPRINT_U0=' + str(U*U0_nd[0]) + '.png');

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
