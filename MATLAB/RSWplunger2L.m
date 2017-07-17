# RSWplunger2L
#=============
# Code for solving the 2-layer plunger-forced governing SW equations.

import numpy as np
import matplotlib.pyplot as plt

import funcs

# All coefficients
# ===========================================================================================================
# ===========================================================================================================

FORCE = 2;        # 1 for manually tranformed delta, 2 for cosine-distributed forcing.
RES = 4;          # 1-8 increasing resolution.

Nsel = [18,34,66,130,258,514,1026,2050];  # Resolution options, 2^n +2.

p = 50									  # Period of the forcing (days).
w = 2 * pi / (3600 * 24 * p);             # Frequency of plunger in (s-1).
ts = (1) * 2 * pi / w;       			  # Plot snapshot time.

N = Nsel[RES];     # Number of gridpoints in y (includes two extra for BCs).
N2 = N-2;		   # Number of gridpoints in x.
TwoN2=2*N2;

Ly=3000000;			# Size of the basin (m).
Lx=3000000;
dy=Ly/(N-1);		# Distance between gridpoints (m).
dx=Lx/(N2-1);

y=np.linspace(-Ly/2,Ly/2,N);     	# Array of all gridpoints.
x=np.linspace(-Lx/2,Lx/2,N2);

K = np.fft.fftfreq(N2,dy)		 # Array of x-gridpoints in wavenumber space.
# Returns an array with 0 as the first wavenumber.
# These differ from the Matlab freqs by a factor or 2pi.

# Remove dead points in the y-direction
yd = np.zeros(N);
for j in range(0,N2):
    yd[j] = y[j+1];

I = complex(0,1);

# Physical parameters
# ===========================================================================================================

U1 = 0.16;      # Background velocity of layer 1 (m s-1)
U2 = 0;      		# Background velocity of layer 2 (m s-1)         

H1 = 250;     # Layer 1 depth (m)
H2 = 3750;    # Layer 2 depth (m)
g = 9.81;     # Acceleration due to gravity (m s-1)

gH1 = g * H1;
gH2 = g * H2;

rho1 = 1000;              # Layer 1 density (kg m-3)
rho2 = 1030;              # Layer 2 density (kg m-3)
rho1prime = rho1 / rho2;                # Reduced density 1
rho2prime = (rho2 - rho1) / rho2;         # Reduced density 2

beta = 2e-11;     # Rate of change of Coriolis (m-1 s-1)
f0 = 0.83e-4;     # Base value of Coriolis parameter (s-1)

f = f0 + beta * y;      # Coriolis paramter (s-1)

# Parameters that depend on k 
# ===========================================================================================================
d1 = np.zeros(N2);
d2 = np.zeros(N2);
d12 = np.zeros(N2);
d22 = np.zeros(N2);
for i in range(0,N2):
    d1[i] = w + U1 * K[i];
    d2[i] = w + U2 * K[i];
    d12[i] = d1[i]**2;
    d22[i] = d2[i]**2;

lambda1 = np.zeros(N2);
zeta1 = np.zeros(N2);
for i in range(0,N2):
    lambda1[i] = d22[i] - gH2 * K[i]**2;
    zeta1[i] = d22[i] - gH2 * K[i]**2 * rho2prime;


# Coefficients appearing directly in the matrix problem
# ===========================================================================================================

a10 = np.zeros((N2,N2));
a12 = np.zeros((N2,N2));
a20 = np.zeros((N2,N2));
a21 = np.zeros((N2,N2));
a22 = np.zeros((N2,N2));
a30 = np.zeros((N2,N2));
a31 = np.zeros((N2,N2));

b10 = np.zeros((N2,N2));
b11 = np.zeros((N2,N2));
b12 = np.zeros((N2,N2));
b20 = np.zeros((N2,N2));
b22 = np.zeros((N2,N2));
b30 = np.zeros((N2,N2));
b31 = np.zeros((N2,N2));

# For wavenumber k=0, the matrix is nearly degenerate, but we can
# substantially simplify the system and find a better-behaved matrix.

# First we define the i=0 rows of the coefficients (wavenumber k=0).
for j in range(0,N2):
    a10[j,0] = f[j+1]**2 - w**2;	# f[j+1] because the array f includes values at the dead endpoints.
    a12[j,0] = - gH1;
    a22[j,0] = - gH2;
    a31[j,0] = - g;
    
    b12[j,0] = - gH1 * rho1prime;
    b20[j,0] = a10[j,0];
    b22[j,0] = - gH2;
    b31[j,0] = - g * rho1prime;

# Now we define the values of the coefficients for all other wavenumbers.
for i in range(1,N2):
    # Define some values that depend only on k, so they aren't computed
    # multiple times.
    
    k = K[i];
    k2 = K[i]**2;
    k3 = K[i]**3;
   
    # The calculation of the coefficients is broken down to maximise efficiency.

    a10i = (d1[i] * gH1 * k2 - gH1 * k * beta) * zeta1[i];   # For a10
    a12i = - d1[i] * gH1 * zeta1[i];                  		 # For a12
    a22i = - d12[i] * d2[i] * gH2;                    		 # For a22
    a31i = - d1[i] * g * zeta1[i];                    		 # For a31
    
    b11i=d2(1,i)*gH1*k2*rho1prime*(U1-U2)*zeta1(1,i);   % For b11
    b12i=-d1(1,i)*d22(1,i)*gH1*rho1prime*zeta1(1,i);   % For b12
    
    b20i1=d12(1,i)*d2(1,i)*lambda1(1,i);                % All for b20
    b20i2=-d12(1,i)*d2(1,i)*lambda1(1,i)*zeta1(1,i);
    b20i3=d2(1,i)*gH1*k2*zeta1(1,i)^2;
    b20i4=d12(1,i)*d2(1,i)*gH2*k2*rho1prime;
    b20i5=-d12(1,i)*d22(1,i)*gH2*k*beta;
    b20i6=d12(1,i)*gH2^2*k3*beta*rho2prime;
    
    b22i=-d2(1,i)*gH2*(d12(1,i)-gH1*k2*rho2prime)*zeta1(1,i); % For b22
    b30i=d1(1,i)*d2(1,i)*g*k*rho1prime*zeta1(1,i);             % For b30
    b31i=-d1(1,i)*d22(1,i)*g*rho1prime*zeta1(1,i);            % For b31
    
    for j=1:N2
        % Coefficients of the first PDE
        % a10
        a10(j,i)=(d1(1,i)*f(j+1,1)^2-d1(1,i)^3)*lambda1(1,i)+a10i;
        % a12
        a12(j,i)=a12i;
        % a20
        a20(j,i)=d1(1,i)*gH2*k*(f(j+1,1)^2*k-d1(1,i)*beta);
        % a21
        a21(j,i)=d1(1,i)*f(j+1,1)*gH2*k2*(U2-U1);
        % a22
        a22(j,i)=a22i;
        % a30
        a30(j,i)=f(j+1,1)*g*k*zeta1(1,i);
        % a31
        a31(j,i)=a31i;
        
        % Coefficients of the second PDE
        % b10
        b10(j,i)=d2(1,i)*gH1*k*rho1prime*(f(j+1,1)^2*k-d2(1,i)*beta)*zeta1(1,i);
        % b11
        b11(j,i)=f(j+1,1)*b11i;
        % b12
        b12(j,i)=b12i;
        % b20
        b20(j,i)=f(j+1,1)^2*b20i1+b20i2...
            +(gH1*gH2*k3*beta*rho2prime-d2(1,i)*f(j+1,1)^2*gH1*k2)*zeta1(1,i)...
            +b20i3+f(j+1,1)^2*b20i4+b20i5+b20i6;
        % b22
        b22(j,i)=b22i;
        % b30
        b30(j,i)=f(j+1,1)*b30i;
        % b31
        b31(j,i)=b31i;
    end
end

% To save time in the matrix problem, divide coefficients by the
% corresponding space-step size: dy^2 for a second-order derivative or 2*dy
% for a first-order derivative.

a12=a12/dy^2;
a21=a21/(2*dy);
a22=a22/dy^2;
a31=a31/(2*dy);

b11=b11/(2*dy);
b12=b12/dy^2;
b22=b22/dy^2;
b31=b31/(2*dy);

% All coefficients are now defined within the confines of the y-loop above
% To tidy up, we can clear some of the intermidiate values that were
% defined in order to define the coeffs of the PDEs

clear a10i a12i a22i a31i b11i b12i b20i1 b20i2 b20i3 b20i4 b20i5 b20i6 b22i b30i b31i

# ===========================================================================================================
# ===========================================================================================================

%Define the RHS of the ODE - the plunger forcing in k-y space
# ===========================================================================================================
% F is the forcing term on the RHS of the original continuity equation.
% Fhat is the Fourier transform of F.
% RHS is a combination of Fhat and its derivative that's on the RHS of the
% vhat differential equation.
% Each of these forcing terms are initially defined on the smaller domain before being
% extended to include the two dead points in the y-direction, e.g. Ffull.

if FORCE==1;
    Fmag=1;
    F=zeros(TwoN2,N2);        % These dimensions as we want to solve the ODEs N2 times, once for each K
    if mod(N2,2)==1;
        F((N2+1)/2,(N2+1)/2)=Fmag;
    else
        %F(N2/2,N2/2)=Fmag;
        F(N2/2+1,N2/2)=Fmag;
    end
    Fhat=zeros(N2,N2);
    if mod(N2,2)==1;
        Fhat((N2+1)/2,:)=Fmag/N2;
    else
        %Fhat(N2/2,:)=Fmag/N2;
        Fhat(N2/2+1,:)=Fmag/N2;
    end
    RHS=zeros(TwoN2,N2);
    for i=1:N2
        for j=1:N2
            RHS(j,i)=-a31(j,i)*Fhat(j,i);       % - signs because its the derivative
            RHS(j+N2,i)=-b31(j,i)*Fhat(j,i);    % of a delta function
        end
    end
    
elseif FORCE==2;
    r0=Ly/32;
    F=zeros(N2,N2);
    for i=1:N2
        for j=1:N2
            r=sqrt(x(i)^2+yd(j)^2);
            if r<r0
                F(j,i)=1e-5*cos((pi/2)*r/r0);       % Magnitude of the forcing arbritrary for now
            else
            end            
        end
    end
    Fhat=fftshift(F,2);
    %Fhat=fft(Fhat,[],2);
    Fhat=N2*ifft(Fhat,[],2,'symmetric');
    Fhat_y=ddy(Fhat,dy);
    RHS=zeros(TwoN2,N2);
    for i=1:N2
        for j=1:N2
            RHS(j,i)=a30(j,i)*Fhat(j,i)+a31(j,i)*Fhat_y(j,i);
            RHS(j+N2,i)=b30(j,i)*Fhat(j,i)+b31(j,i)*Fhat_y(j,i);
        end
    end
    clear Fhat_y
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

# ===========================================================================================================

% Now we solve the problem
# ===========================================================================================================

vd=zeros(TwoN2,N2);     % Initialise the solution, excluding deadpoints

A=zeros(TwoN2,TwoN2);         % Initialise the matrix that characterises the problem

% There is a 
% Northern BCs, equation 1
A(1,1)=a10(1,1)-2*a12(1,1);
A(1,2)=a12(1,1);
A(1,N2+1)=-2*a22(1,1);
A(1,N2+2)=a22(1,1);
    
% Southern BCs, equation 1
A(N2,N2-1)=a12(N2,1);
A(N2,N2)=a10(N2,1)-2*a12(N2,1);
A(N2,TwoN2-1)=a22(N2,1);
A(N2,TwoN2)=-2*a22(N2,1);
    
% Northern BCs, equation 2
A(N2+1,1)=-2*b12(1,1);
A(N2+1,2)=b12(1,1);
A(N2+1,N2+1)=b20(1,1)-2*b22(1,1);
A(N2+1,N2+2)=b22(1,1);

% Southern BCs, equation 2
A(TwoN2,N2-1)=b12(N2,1);
A(TwoN2,N2)=-2*b12(N2,1);
A(TwoN2,TwoN2-1)=b22(2,1);
A(TwoN2,TwoN2)=b20(N2,1)-2*b22(N2,1);

% Assign all of the inner values of the matrix
for j=2:N2-1;
    % Upper left corner
    A(j,j-1)=a12(j,1);
    A(j,j)=a10(j,i)-2*a12(j,1);
    A(j,j+1)=a12(j,1);
     
    % Upper right corner
    A(j,j+N2-1)=a22(j,1);
    A(j,j+N2)=-2*a22(j,1);
    A(j,j+N2+1)=a22(j,1);
    
    % Lower left corner
    A(j+N2,j-1)=b12(j,1);
    A(j+N2,j)=-2*b12(j,1);
    A(j+N2,j+1)=b12(j,1);
     
    % Lower right corner
    A(j+N2,j+N2-1)=b22(j,1);
    A(j+N2,j+N2)=b20(j,1)-2*b22(j,1);
    A(j+N2,j+N2+1)=b22(j,1);
end
    
% If we have a delta-function forcing, we have to multiply through by y
% to avoid delta(y)/y.
if FORCE==1
    for j=1:N2
        A(j,:)=yd(j,1)*A(j,:);
        A(j+N2,:)=yd(j,1)*A(j+N2,:);
    end
end
           
vd(:,1)=linsolve(A,RHS(:,1));       % First solve the k=0 problem.

i=2;
while(i<=N2)
    % Northern BCs, equation 1
    A(1,1)=a10(1,i)-2*a12(1,i);
    A(1,2)=a12(1,i);
    A(1,N2+1)=a20(1,i)-2*a22(1,i);
    A(1,N2+2)=a22(1,i)+a21(1,i);
    
    % Southern BCs, equation 1
    A(N2,N2-1)=a12(N2,i);
    A(N2,N2)=a10(N2,i)-2*a12(N2,i);
    A(N2,TwoN2-1)=a22(N2,i)-a21(N2,i);
    A(N2,TwoN2)=a20(N2,i)-2*a22(N2,i);
    
    % Northern BCs, equation 2
    A(N2+1,1)=b10(1,i)-2*b12(1,i);
    A(N2+1,2)=b12(1,i)+b11(1,i);
    A(N2+1,N2+1)=b20(1,i)-2*b22(1,i);
    A(N2+1,N2+2)=b22(1,i);
    
    % Southern BCs, equation 2
    A(TwoN2,N2-1)=b12(N2,i)-b11(N2,i);
    A(TwoN2,N2)=b10(N2,i)-2*b12(N2,i);
    A(TwoN2,TwoN2-1)=b22(2,i);
    A(TwoN2,TwoN2)=b20(N2,i)-2*b22(N2,i);
    
    % Assign all of the inner values of the matrix
    for j=2:N2-1;
        % Upper left corner
        A(j,j-1)=a12(j,i);
        A(j,j)=a10(j,i)-2*a12(j,i);
        A(j,j+1)=a12(j,i);
        
        % Upper right corner
        A(j,j+N2-1)=a22(j,i)-a21(j,i);
        A(j,j+N2)=a20(j,i)-2*a22(j,i);
        A(j,j+N2+1)=a22(j,i)+a21(j,i);
        
        % Lower left corner
        A(j+N2,j-1)=b12(j,i)-b11(j,i);
        A(j+N2,j)=b10(j,i)-2*b12(j,i);
        A(j+N2,j+1)=b12(j,i)+b11(j,i);
        
        % Lower right corner
        A(j+N2,j+N2-1)=b22(j,i);
        A(j+N2,j+N2)=b20(j,i)-2*b22(j,i);
        A(j+N2,j+N2+1)=b22(j,i);
    end
    
    % If we have a delta-function forcing, we have to multiply through by y
    % to avoid delta(y)/y.
    if FORCE==1
        for j=1:N2
            A(j,:)=yd(j,1)*A(j,:);
            A(j+N2,:)=yd(j,1)*A(j+N2,:);
        end
    end
           
    vd(:,i)=linsolve(A,RHS(:,i));
        
    i=i+1;
end

% Extract v1hat and v2hat from the solution vd
v1hat=zeros(N,N2);
v2hat=zeros(N,N2);
for j=1:N2
    v1hat(j+1,:)=vd(j,:);
    v2hat(j+1,:)=vd(j+N2,:);
end

%%

% We will need the derivative of v1hat and v2hat in order to retrieve all
% other solutions of the 2-L SW system
v1hat_y=ddy(v1hat,dy);
v1hat_yy=ddy(v1hat_y,dy);
v2hat_y=ddy(v2hat,dy);
v2hat_yy=ddy(v2hat,dy);

% Calculate u1hat from v1hat and v2hat and their derivatives
% (Some algebraic manipulations have been made in order to increase
% efficiency when calculating these terms)
u1hat=zeros(N,N2);
for i=1:N2
    k=K(1,i);
    kd=gH1*k^2*zeta1(1,i)-d12(1,i)*lambda1(1,i);
    for j=1:N
        u1hat(j,i)=(d1(1,i)*f(j,1)*lambda1(1,i)/kd)*v1hat(j,i)...
            +(gH1*k*zeta1(1,i)/kd)*v1hat_y(j,i)...
            +(d1(1,i)*f(j,1)*gH2*k^2/kd)*v2hat(j,i)...
            +(d1(1,i)*d2(1,i)*gH2*k/kd)*v2hat_y(j,i);...
            -g*k*zeta1(1,i)/kd*Fhatfull(j,i);
    end
end
u1hat=I*u1hat;

% Now calculate u2hat
u2hat=zeros(N,N2);
for i=1:N2
    k=K(1,i);
    for j=1:N
        u2hat(j,i)=1/zeta1(1,i)*(d1(1,i)*d2(1,i)*rho1prime*u1hat(j,i)...
            +I*d2(1,i)*f(j,1)*rho1prime*v1hat(j,i)...
            -I*d2(1,i)*f(j,1)*v2hat(j,i)...
            -I*gH2*k*rho2prime*v2hat_y(j,i));
    end
end
    
% From the velocities in each layer, we can retrieve the interface heights
eta0hat1=zeros(N,N2);       % We split eta0hat into two parts and use IFFT on each part;
eta0hat2=zeros(N,N2);       % otherwise the code has to subtract infinity from infinity
eta1hat=zeros(N,N2);
for i=1:N2
    k=K(1,i);
    for j=1:N
        eta0hat1(j,i)=-1/(g*k)*d1(1,i)*u1hat(j,i);
        eta0hat2(j,i)=I*f(j,1)/(g*k)*v1hat(j,i);
        eta1hat(j,i)=1/d2(1,i)*(-H2*k*u2hat(j,i)+I*H2*v2hat_y(j,i));
    end
end
eta0hat1(:,1)=0;        % We are evaluating 0/0 in the above loop, which returns NaN.
eta0hat2(:,1)=0;        % Reset the values here so that the the algorithm produces values
%eta0hat1(N,1)=0;         % for the whole domain.
%eta0hat2(N,1)=0;

% Inverse Fourier transforms
# ===========================================================================================================

% Find u1 from its Fourier transform
u1tilde=ifft(u1hat,[],2);  % parameter 2 ensures ifft over rows, rather than columns
u1tilde=ifftshift(u1tilde,2);
u1=real(exp(I*w*ts)*u1tilde);

% Find u2 from its Fourier transform
u2tilde=ifft(u2hat,[],2);
u2tilde=ifftshift(u2tilde,2);
u2=real(exp(I*w*ts)*u2tilde);

% Find v1 from its Fourier transform
v1tilde=ifft(v1hat,[],2);
v1tilde=ifftshift(v1tilde,2);
v1=real(exp(I*w*ts)*v1tilde);

% Find v2 from its Fourier transform
v2tilde=ifft(v2hat,[],2);
v2tilde=ifftshift(v2tilde,2);
v2=real(exp(I*w*ts)*v2tilde);

% Find eta01 from its Fourier transform
eta0tilde=ifft(eta0hat1,[],2)+ifft(eta0hat2,[],2); 
eta0tilde=ifftshift(eta0tilde,2);
eta0=real(exp(I*w*ts)*eta0tilde);

% Find eta1 from its Fourier transform
eta1tilde=ifft(eta1hat,[],2);
eta1tilde=ifftshift(eta1tilde,2);
eta1=real(exp(I*w*ts)*eta1tilde);

% Layer thicknesses
h2=eta1;
h1=eta0-eta1;

%%

% EDDY FORCING & VORTICITY
# ===========================================================================================================

% Now the calculation of the eddy forcing
% Here we calculate the nonlinear Jacobian at every time-step
if EDDYFORCING==1
    
    t0=0;               % Initial time
    T=1*2*pi/w;      % End time
    nt=1000;             % Number of steps to take
    dt=(T-t0)/nt;       % Size of time-step. ()/nt where nt=365 means calculate forcing every day
    
    EF1average=zeros(N,N2);  % Initialise the average of the eddy forcing
    EF2average=zeros(N,N2);
    
    % Start the time loop
    for t=t0:dt:T
        
        u1=real(exp(I*w*t)*u1tilde);
        v1=real(exp(I*w*t)*v1tilde);
                
        u2=real(exp(I*w*t)*u2tilde);
        v2=real(exp(I*w*t)*v2tilde);
        
        eta0=real(exp(I*w*t)*eta0tilde);
        eta1=real(exp(I*w*t)*eta1tilde);
        
        h2=eta1;
        h1=eta0-eta1;
        
        % Vorticity
        q1=ddx(v1,dx)-ddy(u1,dy);
        q2=ddx(v2,dx)-ddy(u2,dy);
        
        % Eddy forcing
        EF1=u1.*ddx(q1,dx)+v1.*ddy(q1,dy);
        EF2=u2.*ddx(q2,dx)+v2.*ddy(q2,dy);
        
        % Amend the average eddy forcing
        EF1average=EF1average+EF1;
        EF2average=EF2average+EF2;
                
    end
    
    EF1average=EF1average/(nt+1);     % Divide by the number of EFs in the sum
    EF2average=EF2average/(nt+1);
    
    EF1zonal=sum(EF1average,2)/N2;    % Calcualate the zonal average of the eddy forcings
    EF2zonal=sum(EF2average,2)/N2;
    
    EF1=EF1-EF1average;
    EF2=EF2-EF2average;
    
    EF1averagelim=max(max(abs(EF1average)));    % Limits to be used for the plots
    EF1lim=max(max(abs(EF1)));
    EF2averagelim=max(max(abs(EF2average)));    % Limits to be used for the plots
    EF2lim=max(max(abs(EF2)));

    surf(x,y,EF1average,'edgecolor','none'); view(0,90); caxis([-EF1averagelim EF1averagelim]); colormap(jet);...
        colorbar; title('$$P_{1}$$','interpreter','latex'); shading interp;
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','EF1av_2L'],'png');

    surf(x,y,EF1,'edgecolor','none'); view(0,90); caxis([-EF1lim EF1lim]); colormap(jet);...
        colorbar; title('$$EF_{1}$$','interpreter','latex'); shading interp;
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','EF1_2L'],'png');
        
    surf(x,y,EF2average,'edgecolor','none'); view(0,90); caxis([-EF2averagelim EF2averagelim]); colormap(jet);...
        colorbar; title('$$P_{2}$$','interpreter','latex'); shading interp;
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','EF2av_2L'],'png');

    surf(x,y,EF2,'edgecolor','none'); view(0,90); caxis([-EF2lim EF2lim]); colormap(jet);...
        colorbar; title('$$EF_{2}$$','interpreter','latex'); shading interp;
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','EF2_2L'],'png');
        
    plot(y,EF1zonal,'linewidth',1.5); view([-90,90]); title('EF1 zonal av'); 
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','EF1zonal_2L'],'png');
    
    plot(y,EF2zonal,'linewidth',1.5); view([-90,90]); title('EF2 zonal av');
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','EF2zonal_2L'],'png');
    
    ts=t;
    
else
end

    
    
    
% For EDDYFORCING=0, still compute vorticity and solutions' amplitude
        
q1=ddx(v1,dx)-ddy(u1,dy);       % Vorticity
q2=ddx(v2,dx)-ddy(u2,dy);
    
u1amp=abs(exp(I*w*ts)*u1tilde);     % All amplitudes
u2amp=abs(exp(I*w*ts)*u2tilde);
v1amp=abs(exp(I*w*ts)*v1tilde);
v2amp=abs(exp(I*w*ts)*v2tilde);
eta0amp=abs(exp(I*w*ts)*eta0tilde);
eta1amp=abs(exp(I*w*ts)*eta1tilde);

h2amp=eta1amp;
h1amp=eta0amp-eta1amp;

% Error
error=Error2L(u1tilde,v1tilde,eta0tilde,u2tilde,v2tilde,eta1tilde,...
    f,w,ts,g,H1,H2,rho1prime,rho2prime,U1,U2,N,N2,dy,dx,Ffull);


%% PLOTS 
# ===========================================================================================================

u1lim=max(max(abs(u1)));
u2lim=max(max(abs(u2)));

v1lim=max(max(abs(v1)));
v2lim=max(max(abs(v2)));

h1lim=max(max(abs(h1)));
h2lim=max(max(abs(h2)));

U1str=num2str(100*U1);

% Plots for u
%figure(1)
surf(x,y,u1,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-u1lim u1lim]);...
    title('u1');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','u1_2L_U1=',U1str,],'png');

%figure(2)    
surf(x,y,u2,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-u2lim u2lim]);...
    title('u2');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','u2_2L_U1=',U1str,],'png');

% Plots for v
%figure(3)
surf(x,y,v1,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-v1lim v1lim]);...
    title('v1');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','v1_2L_U1=',U1str,],'png');

%figure(4)
surf(x,y,v2,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-v2lim v2lim]);...
    title('v2');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','v2_2L_U1=',U1str,],'png');

% Plots for h
%figure(5)
surf(x,y,h1,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-h1lim h1lim]);...
    title('h1');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','h1_2L_U1=',U1str,],'png');

%figure(6)
surf(x,y,h2,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-h2lim h2lim]);...
    title('h2');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','h2_2L_U1=',U1str,],'png');
  
    
% Amplitude Plots

% Plots for u
%figure(1)
surf(x,y,u1amp,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-u1lim u1lim]);...
    title('u1');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','u1amp_2L'],'png');

%figure(2)    
surf(x,y,u2amp,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-u2lim u2lim]);...
    title('u2');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','u2amp_2L'],'png');

% Plots for v
%figure(3)
surf(x,y,v1amp,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-v1lim v1lim]);...
    title('v1');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','v1amp_2L'],'png');

%figure(4)
surf(x,y,v2amp,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-v2lim v2lim]);...
    title('v2');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','v2amp_2L'],'png');

% Plots for h
%figure(5)
surf(x,y,h1amp,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-h1lim h1lim]);...
    title('h1');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','h1amp_2L'],'png');

%figure(6)
surf(x,y,h2amp,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-h2lim h2lim]);...
    title('h2');
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/','h2amp_2L'],'png');

% Final combined plots
figure(8)
subplot(3,2,1)
surf(x,y,u1,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-u1lim u1lim]);...
    title('u1');
subplot(3,2,2)
surf(x,y,u2,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-u2lim u2lim]);...
    title('u2');
subplot(3,2,3)
surf(x,y,v1,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-v1lim v1lim]);...
    title('v1');
subplot(3,2,4)
surf(x,y,v2,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-v2lim v2lim]);...
    title('v2')
subplot(3,2,5)
surf(x,y,h1,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-h1lim h1lim]);...
    title('h1');
subplot(3,2,6)
surf(x,y,h2,'edgecolor','none'); view(0,90); colorbar; colormap(jet);...
    shading interp; axis([-Lx/2 Lx/2 -Ly/2 Ly/2]); view(0,90); caxis([-h2lim h2lim]);...
    title('h2')
saveas(gcf,'2layer.png');

if EDDYFORCING==1
    figure(9)
    subplot(3,2,1)
    surf(x,y,EF1average,'edgecolor','none'); view(0,90); caxis([-EF1averagelim EF1averagelim]); colormap(jet);...
        colorbar; title('$$P_{1}$$','interpreter','LaTex'); shading interp;
    subplot(3,2,2)
    surf(x,y,EF2average,'edgecolor','none'); view(0,90); caxis([-EF2averagelim EF2averagelim]); colormap(jet);...
        colorbar; title('$$P_{2}$$','interpreter','latex'); shading interp;
    subplot(3,2,3)
    surf(x,y,EF1,'edgecolor','none'); view(0,90); caxis([-EF1lim EF1lim]); colormap(jet);...
        colorbar; title('$$EF_{1}$$','interpreter','latex'); shading interp; 
    subplot(3,2,4)
    surf(x,y,EF2,'edgecolor','none'); view(0,90); caxis([-EF2lim EF2lim]); colormap(jet);...
        colorbar; title('$$EF_{2}$$','interpreter','latex'); shading interp;        
    subplot(3,2,5)
    plot(y,EF1zonal,'linewidth',1.5); view([-90,90]); title('EF1 zonal av'); 
    subplot(3,2,6)    
    plot(y,EF2zonal,'linewidth',1.5); view([-90,90]); title('EF2 zonal av');
end

