% code to plot two Rossby plane waves
close all
clear all
load mri

RES=4;
      
% define discretisations
Nsel=[16,32,64,128,256,512,1024,2048];  % Resolution options,

nx = Nsel(RES); % spatial discretisation (number of grids)
ny = nx;

U=0; % Background flow

Ly=3000000;     % Size of the basin
Lx=3000000;

dy=Ly/ny;       % Distance between grid points
dx=Lx/nx;

x=[-Lx/2:dx:Lx/2];
y=[-Ly/2:dx:Ly/2];

beta=2e-11;
f0=0.83e-4;      % Base value of Coriolis parameter

f=f0+beta.*y;       % Coriolis

dt = 0.5*86400; % Time step, n days
ntot = 10; % Total number of iterations


% define variables

A0 = [0.1,0.1]; % Wave amplitude
lambda = [Lx/2 Ly/3; Lx/2 -Ly/3]; % Wavelength in x and y directions

k=2*pi*[1/lambda(1,1),1/lambda(1,2)];    % Zonal wave numbers
l=2*pi*[1/lambda(2,1),1/lambda(2,2)];    % Meridional wave numbers

phi = [0,0]; % Phase

kmag=zeros(1,2);
for j = 1:2;
    kmag(j) = sqrt(k(j)^2+l(j)^2); % Magnitude of wave vector
end

omega = [U*k-beta*k(1,1)/kmag(1)^2,U*k-beta*k(1,2)/kmag(2)^2]; % Wave frequency

A=zeros(nx,ny,2);
A_tot=zeros(nx,ny,ntot);
A1=zeros(nx,ny,ntot);
A2=zeros(nx,ny,ntot);
% Now simulate two propagating waves
n=1;
while(n<=ntot)
    t = n*dt;
    for i = 1:nx;
        for j = 1:ny;            
            for m = 1:2;
                arg_m = k(m)*x(i) + l(m)*y(j) - omega(m)*t + phi(m);
                A(i,j,m)=A0(m)*sin(arg_m);
                
            end
           
        end
    end
    A1(:,:,n)=A(:,:,1);
    A2(:,:,n)=A(:,:,2);
    
    n=n+1;
end

zlim=A0(1)+A0(2);

for n=1:ntot
    A_tot(:,:,n)=A1(:,:,n)+A2(:,:,n);
end

for n = 1:ntot;
    figure(2)
    surf(1:nx,1:ny,A_tot(:,:,n),'edgecolor','none'); view(0,90); colormap(jet); colorbar;...
        axis([1 nx 1 ny]);
    F(n)=getframe(2);
    im=frame2im(F(n));
    [imind,cm] = rgb2ind(im,256);

    if n == 1;

        imwrite(imind,cm,'planewave.gif', 'Loopcount',inf);

    else

        imwrite(imind,cm,'planewave.gif','WriteMode','append');

    end
end

% Now work with a single time snapshot
s=5;
A1=A1(:,:,s);
A2=A2(:,:,s);

% Define the derivatives of A1 and A2
A1_x=ddx(A1,dx);
A1_y=ddy(A1,dy);
A2_x=ddx(A2,dx);
A2_y=ddy(A2,dy);

A1_t=zeros(nx,ny);
A2_t=zeros(nx,ny);
% Now the time-derivative
for i = 1:nx;
    for j = 1:ny;            
        arg_1 = k(1)*x(i) + l(1)*y(j) - omega(1)*s + phi(1);
        A1_t(i,j)=omega(1)*A0(1)*cos(arg_1);
        arg_2 = k(2)*x(i) + l(2)*y(j) - omega(2)*s + phi(2);
        A2_t(i,j)=omega(2)*A0(2)*cos(arg_2);
    end
end

A1_xt=ddx(A1_t,dx);
A1_yt=ddy(A1_t,dy);
A2_xt=ddx(A2_t,dx);
A2_yt=ddy(A2_t,dy);


% Now find u and v from the streamfunction
u=A1_y+A2_y;
v=-A1_x-A2_x;

figure(2)
subplot(1,2,1)
surf(1:nx,1:ny,u,'edgecolor','none'); view(0,90); colormap(jet); colorbar; axis([1 nx 1 ny]);
subplot(1,2,2)
surf(1:nx,1:ny,v,'edgecolor','none'); view(0,90); colormap(jet); colorbar; axis([1 nx 1 ny]);

% Errors - need the time-derivatives of u and v
u_t=A1_yt+A2_yt;
v_t=-A1_xt-A2_xt;

fv=zeros(nx,ny);
fu=zeros(nx,ny);
for j=1:ny
    fv(:,j)=f(j)*v(:,j);
    fu(:,j)=f(j)*u(:,j);
end


erroru=immse(u_t-fv,zeros(nx,ny));
errrorv=immse(v_t+fu,zeros(nx,ny));



