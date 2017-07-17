% Code for solving the pde vhat...

clear all
close all

FORCE=2;        % 1 for manually tranformed delta, 2 for cosine-distributed forcing
RES=4;          % 1-8 increasing resolution

Nsel=[18,34,66,130,258,514,1026,2050];  % Resolution, 2^n +2. Choices 7,8 take a while!

w=1.4544e-4;          % Periodicity of plunger, once every 50 days (e-6)

N=Nsel(RES);     % Number of gridpoints
Ly=3000000;     % Size of the basin
Lx=3000000;
dy=Ly/(N-1);     % Distance between gridpoints
dx=Lx/(N-1);

y=transpose(linspace(-Ly/2,Ly/2,N));    % Array of all gridpoints in physical space
x=linspace(-Lx/2,Lx/2,N);

% Some coefficients
H=4000;     % Ocean depth
g=9.81;     % Gravity
gH=g*H;
beta=2e-11;     % Rate of change of Coriolis
f0=0.83e-4;      % Base value of Coriolis parameter

f=f0+beta.*y;       % Coriolis

%Define the RHS of the ODE - the plunger forcing in k-y space
%==========================================================================
% F is the forcing term on the RHS of the continuity equation.

if FORCE==1;
    Fmag=dy^2/10000;
    F=zeros(N,N);        % These dimensions as we want to solve the ODE N2 times, once for each K
    if mod(N,2)==1;
        F((N+1)/2,(N+1)/2)=Fmag;
    else
        F(N/2,N/2)=Fmag/2;
        F(N/2+1,N/2)=Fmag/2;
    end         
elseif FORCE==2;
    r0=Ly/12;
    F=zeros(N,N);
    for i=1:N
        for j=1:N
            r=sqrt(x(i)^2+y(j)^2);
            if r<r0
                F(j,i)=0.002*cos((pi/2)*r/r0);
            else
            end            
        end
    end   
end

uOld=zeros(N,N);
vOld=zeros(N,N);
hOld=zeros(N,N);

t0=0;
nt=10000000;
T=8*2*pi/w;
dt=(T-t0)/nt;

hlim=0.1;


fvOld=zeros(N,N);
fuOld=zeros(N,N);
t=t0;
for t=t0:nt
    hNew=hOld+dt*(-H*ddx(uOld,dx)-H*ddy(vOld,dy)+F*cos(w*t)-ddx(ddx(hOld,dx),dx)-ddy(ddy(hOld,dy),dy));
    for j=1:N
        fvOld(j,:)=f(j,1)*vOld(j,:);
        fuOld(j,:)=f(j,1)*uOld(j,:);
    end
    uNew=uOld+dt*(fvOld-g*ddx(hOld,dx));
    vNew=vOld+dt*(-fuOld-g*ddy(hOld,dy));
    
    uOld=uNew;
    vOld=vNew;
    hOld=hNew;
    
    figure(1)   
    surf(x,y,hNew,'edgecolor','none'); view(0,90); colormap(jet); colorbar; caxis([-hlim hlim]);
    if mod(t,50)==0
        Fig(t/50+1)=getframe(1);
        im=frame2im(Fig(t/50+1));
        [imind,cm] = rgb2ind(im,256);
        pause(0.001)
    end

    if t == 0;

        imwrite(imind,cm,'SSH.gif', 'Loopcount',inf);

    else

        imwrite(imind,cm,'SSH.gif','WriteMode','append');

    end
    
    
end

    
    
    

