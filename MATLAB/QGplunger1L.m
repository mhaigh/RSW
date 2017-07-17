% Code for solving QG plunger-forced system

clear all
close all

FORCE=1;        % 1 for manually tranformed delta, 2 for cosine-distributed forcing
RES=4;          % 1-8 increasing resolution

Nsel=[18,34,66,130,258,514,1026,2050];  % Resolution, 2^n +2. Choices 7,8 take a while!

w=2*pi/(3600*24*60);      %1.4544e-6;          % Periodicity of plunger, once every 50 days (e-6)
ts=(0.4)*2*pi/w; 

N=Nsel(RES);     % Number of gridpoints
N2=N-2;
Ly=3000000;     % Size of the basin
Lx=3000000;
dy=Ly/(N-1);     % Distance between gridpoints
dx=Lx/(N2-1);

y=linspace(-Ly/2,Ly/2,N);    % Array of all gridpoints in physical space
x=linspace(-Lx/2,Lx/2,N2);

% Array of x-gridpoints in wavenumber space
K=[-N2/2:(N2-1)/2]*2*pi/Lx;
K=fftshift(K,2);      % Shift vector so that zero lies at index 1

% Remove dead points in the y-direction
yd=zeros(N2,1);
for l=1:N2
    yd(l,1)=y(l+1);
end

I=complex(0,1);

beta=2e-11;     % Rate of change of Coriolis

S=0;        % Rossby def radius squared

%Define the RHS of the ODE - the plunger forcing in k-y space
%==========================================================================
% F is the forcing term on the RHS of the original continuity equation.
% Fhat is the Fourier transform of F.
% RHS is a combination of Fhat and its derivative that's on the RHS of the
% vhat differential equation.
% Each of these forcing terms are initially defined on the smaller domain before being
% extended to include the two dead points in the y-direction, e.g. Ffull.

if FORCE==1;
    Fmag=1e-6;
    F=zeros(N2,N2);        % These dimensions as we want to solve the ODE N2 times, once for each K
    if mod(N2,2)==1;
        F((N2+1)/2,(N2+1)/2)=Fmag;
    else
        %F(N2/2,N2/2)=Fmag/2;
        F(N2/2+1,N2/2)=Fmag/2;
    end
    Fhat=zeros(N2,N2);
    if mod(N2,2)==1;
        Fhat((N2+1)/2,:)=Fmag/N2;
    else
        %Fhat(N2/2,:)=Fmag/N2;
        Fhat(N2/2+1,:)=Fmag/N2;
    end
    RHS=Fhat/w;     % For FORCE=1, we solve a real system and then multiply phihat by I.
       
elseif FORCE==2;
    r0=Ly/32;
    F=zeros(N2,N2);
    for i=1:N2
        for j=1:N2
            r=sqrt(x(i)^2+yd(j)^2);
            if r<r0
                F(j,i)=1e-15*cos((pi/2)*r/r0);
            else
            end            
        end
    end
    Fhat=fftshift(F,2);
    Fhat=fft(Fhat,[],2);
        
    RHS=-I*Fhat/w;    
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
   
phid=zeros(N2,N2);     % Initialise the solution, excluding deadpoints

% Here we solve the problem
A=zeros(N2,N2);
i=1;
while(i<=N2)
    k=K(i);
    
    c=k*beta/w-k^2-S;  
        
    % Now define the matrix that characterises the problem
    A(1,1)=c-2/dy^2;;
    A(1,2)=1/dy^2;
    A(N2,N2)=c-2/dy^2;
    A(N2,N2-1)=1/dy^2;
    for j=2:N2-1
        A(j,j)=c-2/dy^2;
        A(j,j-1)=1/dy^2;
        A(j,j+1)=1/dy^2;
    end
    
    if FORCE==1
        for j=1:N2
            A(j,:)=yd(j)*A(j,:);
        end
    end
        
    % Solve the problem for each k
    phid(:,i)=linsolve(A,RHS(:,i));
        
    i=i+1;    
end



% Now extend the solution in the y-direction, to include the deadpoints
phihat=zeros(N,N2);
for j=1:N2
    phihat(j+1,:)=phid(j,:);
end

if FORCE==1
    phihat=I*phihat;
end

phitilde=ifft(phihat,[],2);
phitilde=ifftshift(phitilde,2);
phi=real(phitilde*exp(I*w*ts));

error=Error1LQG(phitilde,beta,w,S,ts,N,N2,dy,dx,Ffull);

u=ddy(phi,dy);
v=-ddx(phi,dx);

% PLOTS
%==========================================================================

ulim=0.6*max(max(abs(u)));
vlim=0.6*max(max(abs(v)));
philim=max(max(abs(phi)));


figure(1)
surf(x,y,phi,'edgecolor','none'); view(0,90); colormap(jet); colorbar; shading interp; %caxis([-philim philim]);

figure(2)
subplot(1,2,1)
surf(x,y,u,'edgecolor','none'); view(0,90); colormap(jet); colorbar; shading interp; caxis([-ulim ulim]);
subplot(1,2,2)
surf(x,y,v,'edgecolor','none'); view(0,90); colormap(jet); colorbar; shading interp; caxis([-vlim vlim]);


