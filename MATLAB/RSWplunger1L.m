% Code for solving the vhat pde of the inviscid 1-L RSW system.

clear all
close all

FORCE=2;        % 1 for manually tranformed delta, 2 for cosine-distributed forcing
RES=4;          % 1-8 increasing resolution
EDDYFORCING=0;  % 1 for calculation of forcing, 0 to turn it off 
OMEGA=2*pi/(3600*24)*[1/10,1/20,1/30,1/40,1/50,1/60,1/70,1/80,1/90,1/100,1/110,1/120];   % The set of forcing frequencies

Nsel=[18,34,66,130,258,514,1026,2050];  % Resolution, 2^n +2. Choices 7,8 take a while!

w=OMEGA(7);      %1.4544e-6;          % Periodicity of plunger, once every 50 days (e-6)
ts=(1)*2*pi/w;       % Time at which to plot the snapshot (if no EDDYFORCING, otherwies it's overwritten)

N=Nsel(RES);     % Number of gridpoints
N2=N-2;
Ly=3000000;     % Size of the basin
Lx=3000000;
dy=Ly/(N-1);     % Distance between gridpoints
dx=Lx/(N2-1);

y=transpose(linspace(-Ly/2,Ly/2,N));    % Array of all gridpoints in physical space
x=linspace(-Lx/2,Lx/2,N2);

% Array of x-gridpoints in wavenumber space
K=[-N2/2:(N2-1)/2]*2*pi/Lx;
K=fftshift(K,2);      % Shift vector so that zero lies at index 1

% Remove dead points in the y-direction
yd=zeros(N2,1);
for j=1:N2
    yd(j,1)=y(j+1,1);
end

I=complex(0,1);

% Some coefficients
H=4000;     % Ocean depth
g=9.81;     % Gravity
gH=g*H;
beta=2e-11;     % Rate of change of Coriolis
f0=0.83e-4;      % Base value of Coriolis parameter

f=f0+beta.*y;      % Coriolis


%Define the RHS of the ODE - the plunger forcing in k-y space
%==========================================================================
% F is the forcing term on the RHS of the original continuity equation.
% Fhat is the Fourier transform of F.
% RHS is a combination of Fhat and its derivative that's on the RHS of the
% vhat differential equation.
% Each of these forcing terms are initially defined on the smaller domain before being
% extended to include the two dead points in the y-direction, e.g. Ffull.

if FORCE==1;
    Fmag=1e-3;
    F=zeros(N2,N2);        % These dimensions as we want to solve the ODE N2 times, once for each K
    if mod(N2,2)==1;
        F((N2+1)/2,(N2+1)/2)=Fmag;
    else
        %F(N2/2,N2/2)=Fmag/2;
        F(N2/2+1,N2/2)=Fmag;
    end
    Fhat=zeros(N2,N2);
    if mod(N2,2)==1;
        Fhat((N2+1)/2,:)=Fmag/N2;
    else
        %Fhat(N2/2,:)=Fmag/N2;
        Fhat(N2/2+1,:)=Fmag/N2;
    end
    RHS=-Fhat/H;
       
elseif FORCE==2;
    r0=Ly/32;
    F=zeros(N2,N2);
    for i=1:N2
        for j=1:N2
            r=sqrt(x(1,i)^2+yd(j,1)^2);
            if r<r0
                F(j,i)=0.02*cos((pi/2)*r/r0)/r0;
            else
            end            
        end
    end
    %F=F-sum(sum(F))/(N2*N2);
    Fhat=fftshift(F,2);
    %Fhat=fft(Fhat,[],2);  % Need to sort this: F is even so Fhat should be real. 
    Fhat=N2*ifft(Fhat,[],2,'symmetric');
    %Fhat=conj(Fhat);      % Also need to check if this conj is needed in general. 
    RHS1=zeros(N2,N2);
    for i=1:N2
        for j=1:N2
            RHS1(j,i)=-f(j+1,1)*K(1,i)*Fhat(j,i)/(w*H);
        end
    end
    RHS2=ddy(Fhat,dy)/H;
            
    RHS=RHS2+RHS1;    
end

% Extend the forcings to fill the whole y-range
RHSfull=zeros(N,N2);
Ffull=zeros(N,N2);
Fhatfull=zeros(N,N2);
for j=1:N2
    RHSfull(j+1,:)=RHS(j,:);
    Ffull(j+1,:)=F(j,:);
    Fhatfull(j+1,:)=Fhat(j,:);
end
    
% Coefficients of the ODE that don't depend on wavenumber
a=-beta^2/gH;
b=-2*f0*beta/gH;

vd=zeros(N2,N2);     % Initialise the solution, excluding deadpoints
c=zeros(1,N2);

% Here we solve the problem
A=zeros(N2,N2);
i=1;
while(i<=N2)
    k=K(1,i);
    
    c(i)=(w^2-f0^2)/gH-k^2+beta*k/w;   % The coefficient of the PDE that depends on k
        
    % Now define the matrix that characterises the problem
    A(1,1)=-2/dy^2+(a*yd(1,1)^2+b*yd(1,1)+c(i));
    A(1,2)=1/dy^2;
    A(N2,N2)=-2/dy^2+(a*yd(N2,1)^2+b*yd(N2,1)+c(i));
    A(N2,N2-1)=1/dy^2;
    for j=2:N2-1
        A(j,j)=-2/dy^2+(a*yd(j,1)^2+b*yd(j,1)+c(i));
        A(j,j-1)=1/dy^2;
        A(j,j+1)=1/dy^2;
    end
    
    if FORCE==1
        for j=1:N2
            A(j,:)=yd(j,1)*A(j,:);
        end
    end
        
    % Solve the problem for each k
    vd(:,i)=linsolve(A,RHS(:,i));
        
    i=i+1;    
end

% Now extend the solution in the y-direction, to include the deadpoints
vhat=zeros(N,N2);
for j=1:N2
    vhat(j+1,:)=vd(j,:);
end
vhatlim=max(max(abs(vhat)));

% Plot the Fourier-physical representation of meridional velocity
%surf(fftshift(vhat,2),'edgecolor','none'); colormap(jet); title('$$\hat{v}$$','interpreter','latex');...
%    view(0,90); xlabel('$$k_{i}$$','interpreter','latex'); ylabel('$$y_{j}$$','interpreter','latex');...
%    shading interp; caxis([-vhatlim vhatlim]); axis([1,N2,1,N]);...
%    set(gca,'XTick',[1 32 64 96 128]); line([65 65],ylim,'color','k','LineStyle','--');...
%    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/1-Layer/Tests/','vhat1L'],'png');
    
    
%pause

waveinfo=wave(vhat,K,Lx);

vtilde=ifft(vhat,[],2);   % parameter 2 ensures ifft over rows, rather than columns
%vtilde=conj(vtilde);
vtilde=ifftshift(vtilde,2);           % shifts zero-frequency back to the centre
v=real(exp(I*w*ts)*vtilde);

%=================

% Find hhat from vhat - currently using a periodic boundary 
hhat=zeros(N,N2);
wk2=zeros(1,N);

% We require the y-derivative of vhat
vhat_y=ddy(vhat,dy);

% A loop to asign the values of hhat
for i=1:N2
    k=K(i);
    wk2=w^2-g*k^2*H;
    for j=1:N
        hhat(j,i)=1/wk2*(w*H*vhat_y(j,i)+f(j,1)*k*H*vhat(j,i)-w*Fhatfull(j,i));
    end
end
hhat=I*hhat;
htilde=ifft(hhat,[],2);
%htilde=conj(htilde);
htilde=ifftshift(htilde,2);
h=real(exp(I*w*ts)*htilde);

%====================

% Four methods of retrieving u:
% 1 - Find u/utilde from original u momentum equation
% 2,3,4 - Find uhat from one of the three expressions for uhat & use IFFT

uOpt=2;
if uOpt==1
    htilde_x=ddx(htilde,dx);
    utilde=zeros(N,N2);
    for i=1:N2
        for j=1:N
            utilde(j,i)=1/w*(f(j,1)*vtilde(j,i)-g*htilde_x(j,i));
        end
    end
    u=real(utilde*exp(I*w*ts));      % vtilde and htilde have already been ifftshifted.
elseif uOpt==2
    uhat=zeros(N,N2);
    for i=1:N2
        k=K(i);
        for j=1:N
            uhat(j,i)=-I*f(j,1)/w*vhat(j,i)-g*k*hhat(j,i);
        end
    end
    utilde=ifft(uhat,[],2);
    %utilde=conj(utilde);
    utilde=ifftshift(utilde,2);
    u=real(exp(I*w*ts)*utilde);
elseif uOpt==3
    uhat=zeros(N,N2);
    hhat_y=ddy(hhat,dy);
    for i=1:N2
        for j=1:N
            uhat(j,i)=-I*w/f(j,1)*vhat(j,i)-g/f(j,1)*hhat_y(j,i);
        end
    end
    utilde=ifft(uhat,[],2);
    utilde=ifftshift(utilde,2);
    u=real(exp(I*w*ts)*utilde);
elseif uOpt==4
    uhat=zeros(N,N2);      % Want to avoid dividing by zero
    for i=2:N2
        k=K(i);
        for j=1:N
            uhat(j,i)=I/k*vhat_y(j,i)-w/(k*H)*hhat(j,i)-I/(k*H)*Fhatfull(j,i);
        end
    end
    utilde=ifft(uhat,[],2);
    utilde=ifftshift(utilde,2);
    u=real(exp(I*w*ts)*utilde);
end

% EDDY FORCING
%==============================================================================

% Now the calculation of the eddy forcing
% Here we calculate the nonlinear Jacobian at every time-step
if EDDYFORCING==1
    
    t0=0;               % Initial time
    T=2*365*3600*24;      % End time
    nt=365;             % Number of steps to take
    dt=(T-t0)/nt;       % Size of time-step. ()/nt where nt=365 means calculate forcing every day
    
    EFaverage=zeros(N,N2);  % Initialise the average of the eddy forcing
    
    % Start the time loop
    for t=t0:dt:T
        
        u=real(exp(I*w*t)*utilde);
        v=real(exp(I*w*t)*vtilde);
        h=real(exp(I*w*t)*htilde);
            
        % Vorticity
        q=ddx(v,dx)-ddy(u,dy);
        
        % Eddy forcing
        EF=u.*ddx(q,dx)+v.*ddy(q,dy);
        
        % Amend the average eddy forcing
        EFaverage=EFaverage+EF;
              
    end
    
    EFaverage=EFaverage/(nt+1);     % Divide by the number of EFs in the sum
    EF=EF-EFaverage;
    
    EFaveragelim=max(max(abs(EFaverage)));
    EFlim=max(max(abs(EF)));

    figure(1)
    surf(x,y,EFaverage,'edgecolor','none'); view(0,90); caxis([-EFaveragelim EFaveragelim]); colormap(jet); colorbar;

    figure(2)
    surf(x,y,EF,'edgecolor','none'); view(0,90); caxis([-EFlim EFlim]); colormap(jet); colorbar;
    
    ts=t;
    
else
    q=ddx(v,dx)-ddy(u,dy);
end

% Error
error=Error1L(utilde,vtilde,htilde,f,w,ts,g,H,N,N2,dy,dx,Ffull);

% Quantification of variation in y
uS=sum(sum(abs(u(1:floor(N/3),:)))); 
uN=sum(sum(abs(u(floor(2*N/3):N,:))));
vS=sum(sum(abs(v(1:floor(N/3),:)))); 
vN=sum(sum(abs(v(floor(2*N/3):N,:))));
hS=sum(sum(abs(h(1:floor(N/3),:)))); 
hN=sum(sum(abs(h(floor(2*N/3):N,:))));
yVar=[uN/uS,vN/vS,hN/hS]
clear uS uN vS vN hS hN
    
%surf(x(i),y(j),v(j,i));

% PLOTS
%==========================================================================

ulim=max(max(abs(u)));
vlim=max(max(abs(v)));
hlim=max(max(abs(h)));

%figure(3)
surf(x,y,u,'edgecolor','none'); colorbar; colormap(jet); title('u');...
    view(0,90); xlabel('x'); ylabel('y'); shading interp; caxis([-ulim ulim]);
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/1-Layer/','u1L'],'png');

%figure(4)
surf(x,y,v,'edgecolor','none'); colorbar; colormap(jet); title('v');...
    view(0,90); xlabel('x'); ylabel('y'); shading interp; caxis([-vlim vlim]);
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/1-Layer/','v1L'],'png');

%figure(5)
surf(x,y,h,'edgecolor','none'); colorbar; colormap(jet); title('h');...
    view(0,90); xlabel('x'); ylabel('y'); shading interp; caxis([-hlim hlim]);
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/1-Layer/','h1L'],'png');

figure(6)
subplot(1,3,1)
surf(x,y,u,'edgecolor','none'); colorbar; colormap(jet); title('u'); axis square;...
    xlabel('x'); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-ulim ulim]);
subplot(1,3,2)
surf(x,y,v,'edgecolor','none'); colorbar; colormap(jet); title('v'); axis square;...
    xlabel('x'); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-vlim vlim]);
subplot(1,3,3)
surf(x,y,h,'edgecolor','none'); colorbar; colormap(jet); title('h'); axis square;...
    xlabel('x'); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); %caxis([-hlim hlim]);


FORCEplot=0;
if FORCEplot==1
    if FORCE==2;
        figure(7)
        surf(x,y,Ffull,'edgecolor','none'); colormap(jet); colorbar; view(0,90);
        figure(8)
        surf(K,y,Fhatfull,'edgecolor','none'); colormap(jet); colorbar; view(0,90);
    else
    end
else
end
clear FORCEplot

movie=0;
if movie==1 
    nt=100;
    t=linspace(0,2*pi/w,nt);
    h=zeros(N,N2,nt);
    for ti=1:nt
        for j=1:N
            for i=1:N2
                h(j,i,ti)=real(htilde(j,i)*exp(I*w*t(ti)));
            end
        end
    end   
                        
    hlim=max(max(max(abs(h))));
    for ti=1:nt
        figure(9)
        surf(x,y,h(:,:,ti),'edgecolor','none'); shading interp; view(0,90); colorbar; colormap(jet);...
            axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); caxis([-hlim hlim]);
     
        mov(ti)=getframe(9);
        im=frame2im(mov(ti));
        [imind,cm] = rgb2ind(im,256);

        if ti == 1;

            imwrite(imind,cm,'planewave.gif', 'Loopcount',inf);

        else

            imwrite(imind,cm,'planewave.gif','WriteMode','append');

        end
    end
end
clear movie

qlim=max(max(abs(q)));
figure(10)
surf(x,y,q,'edgecolor','none'); colormap(jet); colorbar; view(0,90); caxis([-qlim qlim]);

% Plot of c
%Omega=linspace(2*pi/(3600*24*150),2*pi/(3600*24*30),N2);
%Omega=linspace(0,2,N2);
%K=ifftshift(K);
%for i=1:N2
%    for j=1:N2
%        c1(j,i)=(Omega(j)^2-f0^2)/gH-K(i)^2+beta*K(i)/Omega(j);
%    end
%    for j=1:12
%        c(j,i)=(OMEGA(j)^2-f0^2)/gH-K(i)^2+beta*K(i)/OMEGA(j);
%    end
%end
%figure(7)
%subplot(1,2,1);
%surf(K,OMEGA,c,'edgecolor','none'); view(0,90); axis([K(1) K(N2) Omega(1) Omega(N2)]);...
%    colormap(jet); colorbar; 
%subplot(1,2,2);
%surf(K,Omega,c1,'edgecolor','none'); view(0,90); axis([K(1) K(N2) Omega(1) Omega(N2)]);...
%    colormap(jet); colorbar; shading interp;