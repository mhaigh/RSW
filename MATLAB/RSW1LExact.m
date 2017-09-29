% Code for solving the pde vhat...

clear all
close all

RES=6;          % 1-8 increasing resolution

Nsel=[16,32,64,128,256,512,1024,204];  % Resolution, 2^n +2. Choices 7,8 take a while!

I=complex(0,1);
w=1.4544e-4;          % Periodicity of plunger, once every 50 days (e-6)
t=(1.4)*2*pi/w;

N=Nsel(RES);     % Number of gridpoints
Ly=3000000;     % Size of the basin
Lx=3000000;
dy=Ly/(N-1);     % Distance between gridpoints
dx=Lx/(N-1);

y=transpose(linspace(-Ly/2,Ly/2,N));    % Array of all gridpoints in physical space
x=linspace(-Lx/2,Lx/2,N);

% Array of x-gridpoints in wavenumber space
K=[-N/2:(N-1)/2]*2*pi/Lx;
K=fftshift(K,2);      % Shift vector so that zero lies at index 1

% Some coefficients
H=4000;     % Ocean depth
g=9.81;     % Gravity
gH=g*H;
beta=2e-11;     % Rate of change of Coriolis
f0=0.83e-4;      % Base value of Coriolis parameter

f=f0+beta.*y;      % Coriolis

% New coordinates
%====================================
A=(0.25*gH/beta^2)^0.25;


xi=(y+f0/beta)/A;       % New y coordinate

alpha=gH^0.5/(2*beta)*(K.^2-beta.*K/w-(w^2-f0^2)/gH);   % parameter of the ODE

para1=0.5-0.5*alpha;        %para2=0.5;
para3=0.5*xi.^2;

nterms=30; % Number of terms to sum together for the HCF

a=zeros(nterms+1,N);    % depends on alpha, which depends on k
a(1,:)=para1;
b=zeros(nterms+1,1);
b(1,1)=0.5;
n=1;
while(n<=nterms)
    a(n+1,:)=a(n,:).*(para1+n);
    b(n+1,1)=b(n,1)*(n+0.5);
    n=n+1;
end

for n=1:nterms+1
    b(n,1)=b(n,1)*factorial(n-1);
end

para3n=zeros(N,nterms+1);
para3n(:,1)=para3(:,1);
n=1;
while(n<=nterms)
    para3n(:,n+1)=para3n(:,n).*para3(:,1);
    n=n+1;
end

    
vhat=zeros(N,N);
M=zeros(N,N);
n=1;
while(n<=nterms+1)
    for i=1:N
        for j=1:N
            M(j,i)=M(j,i)+a(n,i)*para3n(j,n)/b(n,1);
        end
    end
    n=n+1;
end

for j=1:N
    vhat(j,:)=M(j,:)*exp(-xi(j,1)^2/4);
end


vtilde=ifft(vhat,[],2,'symmetric');  % parameter 2 ensures ifft over rows, rather than columns
v=real(exp(I*w*t)*vtilde);
v=ifftshift(v,2);

figure(2)
surf(K,xi,vhat,'edgecolor','none'); view(0,90);

figure(1)
surf(x,y,v,'edgecolor','none'); view(0,90); colormap(jet); colorbar;
