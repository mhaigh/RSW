% 1L Wavenumber analysis
% This code solves the vhat ODE for a range of omega and finds information
% on the wavenumbers for each omega

clear all
close all

RES=4;          % 1-8 increasing resolution
N_T=100;
T=linspace(40,120,N_T);     % The set of forcing periods

Nsel=[18,34,66,130,258,514,1026,2050];  % Resolution, 2^n +2. Choices 7,8 take a while!

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

% FORCING. 
% The rest of the forcing is defined within a later loop as it depends on w
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
Fhat=fftshift(F,2);
%Fhat=fft(Fhat,[],2);  % Need to sort this: F is even so Fhat should be real. 
Fhat=N2*ifft(Fhat,[],2,'symmetric');
%Fhat=conj(Fhat);      % Also need to check if this conj is needed in general. 
RHS1=zeros(N2,N2);
   
% Coefficients of the ODE that don't depend on wavenumber
a=-beta^2/gH;
b=-2*f0*beta/gH;

vd=zeros(N2,N2);     % Initialise the solution, excluding deadpoints
c=zeros(1,N2);

waveinfo=zeros(N_T,3);
% Here we solve the problem for a range of omega (w)
for ti=1:N_T
    w=2*pi/(3600*24*T(ti));
    A=zeros(N2,N2);
    i=1;
    
   % Define the RHS of the ODE - the plunger forcing in k-y space

   for ii=1:N2
        for j=1:N2
            RHS1(j,ii)=-f(j+1,1)*K(1,ii)*Fhat(j,ii)/(w*H);
        end
    end
    RHS2=ddy(Fhat,dy)/H;
            
    RHS=RHS2+RHS1;    

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
           
        % Solve the problem for each k
        vd(:,i)=linsolve(A,RHS(:,i));
        
        i=i+1;    
    end
    
    % Now extend the solution in the y-direction, to include the deadpoints
    vhat=zeros(N,N2);
    for j=1:N2
        vhat(j+1,:)=vd(j,:);
    end
     
    waveinfo(ti,:)=wave(vhat,K,Lx);
end

wavesmooth=zeros(N_T,2);
% Can't average the first and last values
wavesmooth(1,1)=waveinfo(1,1);
wavesmooth(1,2)=waveinfo(1,2);
wavesmooth(N_T,1)=waveinfo(N_T,1);
wavesmooth(N_T,2)=waveinfo(N_T,2);
% Average the second and second-last values over 3 samples
wavesmooth(2,1)=(waveinfo(1,1)+waveinfo(2,1)+waveinfo(3,1))/3;
wavesmooth(2,2)=(waveinfo(1,2)+waveinfo(2,2)+waveinfo(3,2))/3;
wavesmooth(N_T-1,1)=(waveinfo(N_T-2,1)+waveinfo(N_T,1)+waveinfo(N_T,1))/3;
wavesmooth(N_T-1,2)=(waveinfo(N_T-2,2)+waveinfo(N_T,2)+waveinfo(N_T,2))/3;
% Average over all other values
for ti=3:N_T-2
    wavesmooth(ti,1)=(waveinfo(ti-2,1)+waveinfo(ti-1,1)+waveinfo(ti,1)+waveinfo(ti+1,1)+waveinfo(ti+2,1))/5;
    wavesmooth(ti,2)=(waveinfo(ti-2,2)+waveinfo(ti-1,2)+waveinfo(ti,2)+waveinfo(ti+1,2)+waveinfo(ti+2,2))/5;
end


figure(1)
plot(T,waveinfo(:,1),'k',T,wavesmooth(:,1),'r','LineWidth',1.2);...
    xlabel('Forcing Period (days)'); ylabel('Dominant Wavenumber');...
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/1-Layer/Tests/','k1'],'png');
figure(2)
plot(T,waveinfo(:,2),'k',T,wavesmooth(:,2),'r','LineWidth',1.2);...
    xlabel('Forcing Period (days)'); ylabel('Second-most Dominant Wavenumber');...
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/1-Layer/Tests/','k2'],'png');

