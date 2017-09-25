% Code for solving the vhat pde of the viscous 1-L RSW system.

clear all
close all

FORCE=2;        % 1 for manually tranformed delta, 2 for cosine-distributed forcing
RES=4;          % 1-8 increasing resolution
EDDYFORCING=0;  % 1 for calculation of forcing, 0 to turn it off 
OMEGA=2*pi/(3600*24)*[1/10,1/20,1/30,1/40,1/50,1/60,1/70,1/80,1/90,1/100,1/110,1/120];   % The set of forcing frequencies

Nsel=[18,34,66,130,258,514,1026,2050];  % Resolution, 2^n +2. Choices 7,8 take a while!

w=OMEGA(6);      %1.4544e-6;          % Periodicity of plunger, once every 50 days (e-6)
ts=(1)*2*pi/w;       % Time at which to plot the snapshot (if no EDDYFORCING, otherwise it's overwritten)

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

nu=2000;         % Viscosity
  
% Coefficients of the matrix
%=================================================

% Coefficients of the ODE that don't depend on k or y
a2=-w*nu;
b4=I*gH-nu*w;
b5=I*g;

% All other coefficients
a1=zeros(1,N2);
a3=zeros(N2,1);
a4=zeros(1,N2);
a5=zeros(1,N2);
b1=zeros(N2,1);
b2=zeros(1,N2);
b3=zeros(1,N2);

% Coefficients that depend on k
for i=1:N2
    k=K(1,i);
    a1(1,i)=I*(w^2-gH*k^2)+w*nu*k^2;
    a4(1,i)=-gH*k;
    a5(1,i)=-g*k;
    b2(1,i)=-gH*k;
    b3(1,i)=I*w^2+k^2*w*nu;
end

% Coefficients that depend on y
for j=1:N2
    a3(j,1)=-w*f(j+1,1);
    b1(j,1)=w*f(j+1,1);
end

% Normalise the coefficients by their respective step-sizes
a2=a2/dy^2;
a4=a4/(2*dy);
b2=b2/(2*dy);
b4=b4/dy^2;
b5=b5/(2*dy);

%==========================================================================


%Define the RHS of the ODE - the plunger forcing in k-y space
%==========================================================================
% F is the forcing term on the RHS of the original continuity equation.
% Fhat is the Fourier transform of F.
% RHS is a combination of Fhat and its derivative that's on the RHS of the
% vhat differential equation.
% Each of these forcing terms are initially defined on the smaller domain before being
% extended to include the two dead points in the y-direction, e.g. Ffull.

if FORCE==1;
    Fmag=1e-8;
    F=zeros(N2,N2);        % These dimensions as we want to solve the ODE N times, once for each K
    if mod(N2,2)==1;
        F((N2+1)/2,(N2+1)/2)=Fmag;
    else
        %F(N/2,N/2)=Fmag;
        F(N2/2+1,N2/2)=Fmag;
    end
    Fhat=zeros(N2,N2);
    if mod(N2,2)==1;
        Fhat((N2+1)/2,:)=Fmag/N;
    else
        %Fhat(N/2,:)=Fmag/N;
        Fhat(N2/2+1,:)=Fmag/N;
    end
    RHS=zeros(N+N2,N2);
    for j=1:N2
        RHS(N2+j,:)=-g*Fhat(j,1);  % minus sign because derivative of delta introduces it
    end
    
       
elseif FORCE==2;
    r0=Ly/32;
    F=zeros(N2,N2);
    for i=1:N2
        for j=1:N2
            r=sqrt(x(i)^2+yd(j)^2);
            if r<r0
                F(j,i)=1e-6*cos((pi/2)*r/r0);
            else
            end            
        end
    end
    %F=F-sum(sum(F))/(N2*N2);
    Fhat=fftshift(F,2);
    %Fhat=fft(Fhat,[],2);  % Need to sort this: F is even so Fhat should be real. 
    Fhat=N2*ifft(Fhat,[],2,'symmetric');
    Fhat_y=ddy(Fhat,dy);
    RHS=zeros(2*N2,N2);
    for i=1:N2
        for j=1:N2
            RHS(j,i)=a5(1,i)*Fhat(j,i);
            RHS(N2+j,i)=b5*Fhat_y(j,i);
        end
    end
        
end

% Extend the forcings to fill the whole y-range
RHSfull=zeros(2*N,N2);
Ffull=zeros(N,N2);
Fhatfull=zeros(N,N2);
for j=1:N2
    RHSfull(j+1,:)=RHS(j,:);
    RHSfull(N+j+1,:)=RHS(N2+j,:);
    Ffull(j+1,:)=F(j,:);
    Fhatfull(j+1,:)=Fhat(j,:);
end

solution=zeros(2*N2,N2);     % Initialise the solution, excluding deadpoints:
                             % vhat has initial y-dimension=N2, 
                             % uhat has initial y-dimension=N.

% Here we solve the problem for i=N2 wavenumbers
A=zeros(2*N2,2*N2);
i=1;
while(i<=N2)        
    
    % Define the matrix that characterises the problem
    
    % Upper left corner BCs
    A(1,1)=a1(1,i)-2*a2;
    A(1,2)=a2;
    A(N2,N2)=a1(1,i)-2*a2;
    A(N2,N2-1)=a2;
    % Lower left corner BCs
    A(N2+1,1)=b1(1,1);
    A(N2+1,2)=b2(1,i);
    A(2*N2,N2)=b1(N2,1);
    A(2*N2,N2-1)=-b2(1,i);
    % Upper right corner BCs
    A(1,N2+1)=a3(1,1);
    A(1,N2+2)=a4(1,i);
    A(N2,2*N2)=a3(N2,1);
    A(N2,2*N2-1)=-a4(1,i);
    % Lower right corner BCs
    A(N2+1,N2+1)=b3(1,i)-2*b4;
    A(N2+1,N2+2)=b4;
    A(2*N2,2*N2)=b3(1,i)-2*b4;
    A(2*N2,2*N2-1)=b4;
    
    % Assign the remainder of the matrix's values
    for j=2:N2-1
        % Upper left corner
        A(j,j)=a1(1,i)-2*a2;
        A(j,j-1)=a2;
        A(j,j+1)=a2;
        % Lower left corner
        A(N2+j,j)=b1(j,1);
        A(N2+j,j-1)=-b2(1,i);
        A(N2+j,j+1)=b2(1,i);
        % Upper right corner
        A(j,N2+j)=a3(j,1);
        A(j,N2+j-1)=-a4(1,i);
        A(j,N2+j+1)=a4(1,i);
        % Lower right corner
        A(N2+j,N2+j)=b3(1,i)-2*b4;
        A(N2+j,N2+j-1)=b4;
        A(N2+j,N2+j+1)=b4;
    end
    
    if FORCE==1
        for j=1:N2
            A(j,:)=y(j,1)*A(j,:);
            A(N2+j,:)=y(j,1)*A(N2+j,:);
        end
    end
        
    % Solve the problem for each k
    
    solution(:,i)=linsolve(A,RHS(:,i));
        
    i=i+1;    
end

% Extract the solution for uhat and vhat. Include dead points on which we
% have no-slip and no-normal flow
vhat=zeros(N,N2);
uhat=zeros(N,N2);
for j=1:N2
    vhat(j,:)=solution(N2+j,:);
    uhat(j,:)=solution(j,:);
end

vhat_y=ddy(vhat,dy);

hhat=zeros(N,N2);
for i=1:N2
    k=K(1,i);
    for j=1:N
        hhat(j,i)=-H*k*uhat(j,i)/w+I*H*vhat_y(j,i)/w-I*Fhatfull(j,i)/w;
    end
end


vtilde=ifft(vhat,[],2);  % parameter 2 ensures ifft over rows, rather than columns
vtilde=ifftshift(vtilde,2); % shifts zero-frequency back to the centre
v=real(exp(I*w*ts)*vtilde);        

utilde=ifft(uhat,[],2); 
utilde=ifftshift(utilde,2);
u=real(exp(I*w*ts)*utilde);       

htilde=ifft(hhat,[],2); 
htilde=ifftshift(htilde,2);
h=real(exp(I*w*ts)*htilde);

%h(1,:)=zeros(1,N);
%h(N,:)=zeros(1,N);

%error=Error1Lvisc(utilde,vtilde,htilde,nu,f,w,t,g,H,N,N,dy,dx,Ffull);

% PLOTS
%==========================================================================

ulim=max(max(abs(u)));
vlim=max(max(abs(v)));
hlim=max(max(abs(h)));

figure(8)
subplot(1,3,1)
surf(x,y,u,'edgecolor','none'); colorbar; colormap(jet); title('u'); axis square;...
    xlabel('x'); ylabel('y'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-ulim ulim]);
subplot(1,3,2)
surf(x,y,v,'edgecolor','none'); colorbar; colormap(jet); title('v'); axis square;...
    xlabel('x'); ylabel('y'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-vlim vlim]);
subplot(1,3,3)
surf(x,y,h,'edgecolor','none'); colorbar; colormap(jet); title('h'); axis square;...
    xlabel('x'); ylabel('y'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 -hlim hlim]); view(0,90); caxis([-hlim hlim]);

% Plots for v
%figure(1)
%subplot(2,1,1)
%surf(K,y,vhat,'edgecolor','none'); colorbar; colormap(jet); title('Fourier transform of meridional velocity')
%subplot(2,1,2)
%surf(K,y,vhat,'edgecolor','none'); view(0,90);

%figure(2)
%subplot(2,2,[1,2])
%surf(x,y,v,'edgecolor','none'); colorbar; colormap(jet); title('Meridional velocity response to plunger forcing');...
%    xlabel('x'); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 min(min(v)) max(max(v))]);
%subplot(2,2,3);
%surf(x,y,v,'edgecolor','none'); view(90,0); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 min(min(v)) max(max(v))]);
%subplot(2,2,4);
%surf(x,y,v,'edgecolor','none'); title('Top view'); xlabel('x'); ylabel('y'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 min(min(v)) max(max(v))]); view(0,90);

% Plots for h
%figure(3)
%subplot(2,1,1)
%surf(K,y,imag(hhat),'edgecolor','none'); colorbar; colormap(jet); title('Fourier transform of SSH')
%subplot(2,1,2)
%surf(K,y,imag(hhat),'edgecolor','none')
%view(0,90);

%figure(4)
%subplot(2,2,[1,2])
%surf(x,y,h,'edgecolor','none'); colorbar; colormap(jet); title('SSH response to plunger forcing');...
%    xlabel('x'); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 min(min(h)) max(max(h))]);
%subplot(2,2,3);
%surf(x,y,h,'edgecolor','none'); view(90,0); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 min(min(h)) max(max(h))]);
%subplot(2,2,4);
%surf(x,y,h,'edgecolor','none'); title('Top view'); xlabel('x'); ylabel('y'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 min(min(h)) max(max(h))]); view(0,90);

% Plots for u
%figure(5)
%subplot(2,1,1)
%surf(K,y,imag(uhat),'edgecolor','none'); colorbar; colormap(jet); title('Fourier transform of zonal velocity')
%subplot(2,1,2)
%surf(K,y,imag(uhat),'edgecolor','none')
%view(0,90);

%figure(6)
%subplot(2,2,[1,2])
%surf(x,y,u,'edgecolor','none'); colorbar; colormap(jet); title('Zonal velocity response to plunger forcing');...
%    xlabel('x'); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 min(min(u)) max(max(u))]);
%subplot(2,2,3);
%surf(x,y,u,'edgecolor','none'); view(90,0); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 min(min(u)) max(max(u))]);
%subplot(2,2,4);
%surf(x,y,u,'edgecolor','none'); title('Top view'); xlabel('x'); ylabel('y'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 min(min(u)) max(max(u))]);                                                                            view(0,90);