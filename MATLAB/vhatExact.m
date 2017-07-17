% Code for solving the pde vhat...

clear all
close all

RES=4;          % 1-8 increasing resolution

Nsel=[18,34,66,130,258,514,1026,2050];  % Resolution, 2^n +2. Choices 7,8 take a while!

w=1.4544e-6;          % Periodicity of plunger, once every 50 days (e-6)
t=(1.4)*2*pi/w;

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

% Some coefficients
H=4000;     % Ocean depth
g=9.81;     % Gravity
gH=g*H;
beta=2e-11;     % Rate of change of Coriolis
f0=0.83e-4;      % Base value of Coriolis parameter

f=transpose(f0+beta.*y);      % Coriolis

% New coordinates
%====================================

xi=(0.25*gH/beta^2)^0.25.*y-f0/beta;       % New y coordinate

alpha=gH/(2*beta)*(K.^2-beta.*K/w-(w^2-f0^2)/gH); % parameter of the ODE