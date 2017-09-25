% Code for solving the pde vhat...

clear all
close all

FORCE=1;        % 1 for manually tranformed delta, 2 for cosine-distributed forcing
RES=5;          % 1-8 increasing resolution

Nsel=[18,34,66,130,258,514,1026,2050];  % Resolution, 2^n +2. Choices 7,8 take a while!

w=1.4544e-4;          % Periodicity of plunger, once every 50 days (e-6)
t=(1.4)*2*pi/w;

wt=w*tan(w*t);

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
    Fmag=dy^2/100;
    F=zeros(N2,N2);        % These dimensions as we want to solve the ODE N2 times, once for each K
    if mod(N2,2)==1;
        F((N2+1)/2,(N2+1)/2)=Fmag;
    else
        F(N2/2,N2/2)=Fmag/2;
        F(N2/2+1,N2/2)=Fmag/2;
    end
    Fhat=zeros(N2,N2);
    if mod(N2,2)==1;
        Fhat((N2+1)/2,:)=Fmag;
    else
        Fhat(N2/2,:)=Fmag/2;
        Fhat(N2/2+1,:)=Fmag/2;
    end
    RHS=-Fhat/H;
       
elseif FORCE==2;
    r0=Ly/64;
    F=zeros(N2,N2);
    for i=1:N2
        for j=1:N2
            r=sqrt(x(i)^2+yd(j)^2);
            if r<r0
                F(j,i)=0.002*cos((pi/2)*r/r0);
            else
            end            
        end
    end
    Fhat=fftshift(F,2);
    Fhat=fft(Fhat,[],2);
    RHS1=zeros(N2,N2);
    for i=1:N2
        for j=1:N2
            RHS1(j,i)=f(1,j+1)*K(1,i)*Fhat(j,i)/(wt*H);
        end
    end
    RHS2=ddy(Fhat,dy);
    for j=1:N2
        RHS2(j,:)=RHS2(j,:)/H;
    end
        
    RHS=dy*dy*(RHS2-RHS1);    
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

% Here we solve the problem
A=zeros(N2,N2);
i=1;
while(i<=N2)
    k=K(i);
    
    c=-(wt^2+f0^2)/gH-k^2+I*beta*k/wt;   % The coefficient of the PDE that depends on k
        
    % Now define the matrix that characterises the problem
    A(1,1)=-2+dy^2*(a*yd(1)^2+b*yd(1)+c);
    A(1,2)=1;
    A(N2,N2)=-2+dy^2*(a*yd(N2)^2+b*yd(N2)+c);
    A(N2,N2-1)=1;
    for j=2:N2-1
        A(j,j)=-2+dy^2*(a*yd(j)^2+b*yd(j)+c);
        A(j,j-1)=1;
        A(j,j+1)=1;
    end
    
    if FORCE==1
        for j=1:N2
            A(j,:)=yd(j)*A(j,:);
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

vtilde=ifft(vhat,[],2);  % parameter 2 ensures ifft over rows, rather than columns
v=real(cos(w*t)*vtilde);
v=ifftshift(v,2);          % shifts zero-frequency back to the centre
%=================

% Find hhat from vhat - currently using a periodic boundary 
hhat=zeros(N,N2);
wk2=zeros(1,N);

% We require the y-derivative of vhat
vhat_y=ddy(vhat,dy);

% A loop to asign the values of hhat
for i=1:N2
    k=K(i);
    wk2=wt^2+g*k^2*H;
    for j=1:N
        hhat(j,i)=1/wk2*(wt*H*vhat_y(j,i)-I*f(1,j)*H*k*vhat(j,i)-wt*Fhatfull(j,i));
    end
end
htilde=ifft(hhat,[],2);
h=real(cos(w*t)*htilde);           % Having trouble with the ifft and taking real part
h=ifftshift(h,2);
%====================

% Four methods of retrieving u:
% 1 - Find u/utilde from original u momentum equation
% 2,3,4 - Find uhat from one of the three expressions for uhat & use IFFT

uOpt=1;
if uOpt==1
    htilde_x=ddx(htilde,dx);
    utilde=zeros(N,N2);
    for i=1:N2
        for j=1:N
            utilde(j,i)=1/w*(f(1,j)*vtilde(j,i)-g*htilde_x(j,i));
        end
    end
    u=-real(utilde.*sin(w*t));
    u=ifftshift(u,2);
elseif uOpt==2
    uhat=zeros(N,N2);
    for i=1:N2
        k=K(i);
        for j=1:N
            uhat(j,i)=-I*f(1,j)/w*vhat(j,i)-g*k*hhat(j,i);
        end
    end
    utilde=ifft(uhat,[],2);
    u=real(exp(I*w*t)*utilde);
    u=ifftshift(u,2);
elseif uOpt==3
    uhat=zeros(N,N2);
    hhat_y=ddy(hhat,dy);
    for i=1:N2
        for j=1:N
            uhat(j,i)=-I*w/f(1,j)*vhat(j,i)-g/f(1,j)*hhat_y(j,i);
        end
    end
    utilde=ifft(uhat,[],2);
    u=real(exp(I*w*t)*utilde);
    u=ifftshift(u,2);
elseif uOpt==4
    uhat=zeros(N,N2);      % Want to avoid dividing by zero
    for i=2:N2
        k=K(i);
        for j=1:N
            uhat(j,i)=I/k*vhat_y(j,i)-w/(k*H)*hhat(j,i)-I/(k*H)*Fhatfull(j,i);
        end
    end
    utilde=ifft(uhat,[],2);
    u=real(exp(I*w*t)*utilde);
    u=ifftshift(u,2);
end

% Vorticity
q=zeros(N,N2);
zeta=ddx(v,dx)-ddy(u,dy);
for i=1:N2
    for j=1:N
        q(j,i)=(zeta(j,i))/(H+h(j,i));
    end
end

% Eddy forcing
uq=u.*q;
vq=v.*q;

E=-(ddx(uq,dx)+ddy(vq,dy));

Elim=0.1*max(max(abs(E)));
figure(11)
surf(x,y,E,'edgecolor','none'); view(0,90); caxis([-Elim Elim]); colormap(jet); colorbar;


% Error
error=Error1L(u,v,h,f,w,t,g,H,N,N2,dy,dx,Ffull);

%surf(x(i),y(j),v(j,i));

% PLOTS
%==========================================================================

ulim=max(max(abs(u)));
vlim=max(max(abs(v)));
hlim=max(max(abs(h)));

figure(8)
subplot(1,3,1)
surf(x,y,u,'edgecolor','none'); colorbar; colormap(jet); title('u'); axis square;...
    xlabel('x'); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-ulim ulim]);
subplot(1,3,2)
surf(x,y,v,'edgecolor','none'); colorbar; colormap(jet); title('v'); axis square;...
    xlabel('x'); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-vlim vlim]);
subplot(1,3,3)
surf(x,y,h,'edgecolor','none'); colorbar; colormap(jet); title('h'); axis square;...
    xlabel('x'); ylabel('y'); zlabel('z'); shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2 -hlim hlim]); view(0,90); caxis([-hlim hlim]);

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

FORCEplot=1;
if FORCEplot==1
    if FORCE==2;
        figure(7)
        surf(x,y,Ffull,'edgecolor','none'); colormap(jet); colorbar; view(0,90);
    else
    end
else
end

movie=0;
if movie==1  
    dt=0.1/w;
    minu=min(min(u));
    maxu=max(max(u));
    for t=0:dt:12*pi/w
        figure(9)
        surf(x,y,u*sin(w*t),'edgecolor','none'); view(0,90); shading interp; colorbar; colormap(jet); legend('$t$');...
            title('Zonal velocity response to plunger forcing'); xlabel('x'); ylabel('y'); zlabel('z'); axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); caxis([-ulim ulim]);
    end
end

qlim=max(max(abs(q)));
figure(10)
surf(x,y,q,'edgecolor','none'); colormap(jet); colorbar; view(0,90); caxis([-qlim qlim]);