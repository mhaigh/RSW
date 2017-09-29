% 2L Wavenumber analysis
% This code solves the v1hat, v2hat ODE for a range of omega and finds information
% on the wavenumbers for each omega

%clear all
close all

test=1; % 0 for plotting wavenumber vs frequency, 1 for vs U1

RES=4;          % 1-8 increasing resolution

if test==0
    N_T=120;
    T=linspace(40,120,N_T);     % The set of forcing periods
    U1=-0.08;                   % Set the BG flow
    waveinfo1=zeros(N_T,3);
    waveinfo2=zeros(N_T,3);
end
if test==1
    Ti=4;
    N_T=50; 
    U1range=linspace(-0.16,0.32,N_T);
    Trange=[50,60,70,80];       % Here the period is kept fixed, four choices
    T=Trange(Ti);     
    w=2*pi/(3600*24*T);
    if Ti==1;
        waveinfo1=zeros(N_T,3,4);
        waveinfo2=zeros(N_T,3,4);
    end
end

Nsel=[18,34,66,130,258,514,1026,2050];  % Resolution options, 2^n +2. Choices 7,8 take a while!

N=Nsel(RES);     % Number of gridpoints
N2=N-2;
TwoN2=2*N2;
Ly=3000000;     % Size of the basin
Lx=3000000;
dy=Ly/(N-1);     % Distance between gridpoints
dx=Lx/(N2-1);

y=linspace(-Ly/2,Ly/2,N);    % Array of all gridpoints
x=linspace(-Lx/2,Lx/2,N2);

% Array of x-gridpoints in wavenumber space
K=[-N2/2:(N2-1)/2]*2*pi/Lx;
K=fftshift(K,2);      % Shift vector so that zeros lies at index 1

% Remove dead points in the y-direction
yd=zeros(N2,1);
for j=1:N2
    yd(j,1)=y(1,j+1);
end

I=complex(0,1);

% Coefficients
%=========================================================================
%=========================================================================

% Some coefficients
% U1 is defined at the top of the code
U2=0;      % Background velocity of layer 2          

H1=250;     % Layer 1 depth
H2=3750;     % Layer 2 depth
g=9.81;      % Gravity

gH1=g*H1;
gH2=g*H2;

rho1=1000;              % Layer 1 density km m-3
rho2=1030;              % Layer 2 density
rho1prime=rho1/rho2;                % Reduced density 1
rho2prime=(rho2-rho1)/rho2;         % Reduced density 2

beta=2e-11;     % Rate of change of Coriolis
f0=0.83e-4;      % Base value of Coriolis parameter

f=transpose(f0+beta.*y);       % Coriolis

for ti=1:N_T
    if test==0
       w=2*pi/(3600*24*T(ti));
    end
    if test==1
        U1=U1range(ti);
    end
    
    % Coefficients that depend on k
    d1=zeros(1,N2);
    d2=zeros(1,N2);
    d12=zeros(1,N2);
    d22=zeros(1,N2);
    for i=1:N2
        d1(1,i)=w+U1*K(1,i);
        d2(1,i)=w+U2*K(1,i);
        d12(1,i)=d1(1,i)^2;
        d22(1,i)=d2(1,i)^2;
    end

    lambda=zeros(1,N2);
    zeta1=zeros(1,N2);
    for i=1:N2
        lambda(1,i)=d22(1,i)-gH2*K(1,i)^2;
        zeta1(1,i)=d22(1,i)-gH2*K(1,i)^2*rho2prime;
    end

    % Coefficients appearing directly in the matrix problem
    %=========================================================================

    a10=zeros(N2,N2);
    a12=zeros(N2,N2);
    a20=zeros(N2,N2);
    a21=zeros(N2,N2);
    a22=zeros(N2,N2);
    a30=zeros(N2,N2);
    a31=zeros(N2,N2);

    b10=zeros(N2,N2);
    b11=zeros(N2,N2);
    b12=zeros(N2,N2);
    b20=zeros(N2,N2);
    b22=zeros(N2,N2);
    b30=zeros(N2,N2);
    b31=zeros(N2,N2);

    % For wavenumber k=0, the matrix is nearly degenerate, but we can
    % substantially simplify the system and find a better-behaved matrix.

    % First we define the i=1 columns of the coefficients.
    for j=1:N2
        a10(j,1)=f(j+1,1)^2-w^2;
        a12(j,1)=-gH1;
        a22(j,1)=-gH2;
        a31(j,1)=-g;
    
        b12(j,1)=-gH1*rho1prime;
        b20(j,1)=a10(j,1);
        b22(j,1)=-gH2;
        b31(j,1)=-g*rho1prime;
    end

    % Now we define the values of the coefficients for all other wavenumbers
    for i=2:N2
        % Define some values that depend only on k, so they aren't computed
        % multiple times.
    
        k=K(1,i);
        k2=K(1,i).^2;
        k3=K(1,i).^3;
   
    
        a10i=(d1(1,i)*gH1*k2-gH1*k*beta)*zeta1(1,i);   % For a10
        a12i=-d1(1,i)*gH1*zeta1(1,i);                 % For a12
        a22i=-d12(1,i)*d2(1,i)*gH2;                   % For a22
        a31i=-d1(1,i)*g*zeta1(1,i);                   % For a31
    
        b11i=d2(1,i)*gH1*k2*rho1prime*(U1-U2)*zeta1(1,i);   % For b11
        b12i=-d1(1,i)*d22(1,i)*gH1*rho1prime*zeta1(1,i);   % For b12
    
        b20i1=d12(1,i)*d2(1,i)*lambda(1,i);                % All for b20
        b20i2=-d12(1,i)*d2(1,i)*lambda(1,i)*zeta1(1,i);
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
            a10(j,i)=(d1(1,i)*f(j+1,1)^2-d1(1,i)^3)*lambda(1,i)+a10i;
            % a12
            a12(j,i)=a12i;
            % a20
            a20(j,i)=d1(1,i)*gH2*k*(f(j+1,1)^2*k-d1(1,i)*beta);
            % a21
            a21(j,i)=d1(1,i)*f(j,1)*gH2*k2*(U2-U1);
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

    %===========================================================================
    %===========================================================================

    %Define the RHS of the ODE - the plunger forcing in k-y space
    %==========================================================================
    % F is the forcing term on the RHS of the original continuity equation.
    % Fhat is the Fourier transform of F.
    % RHS is a combination of Fhat and its derivative that's on the RHS of the
    % vhat differential equation.
    % Each of these forcing terms are initially defined on the smaller domain before being
    % extended to include the two dead points in the y-direction, e.g. Ffull.

    r0=Ly/32;
    F=zeros(N2,N2);
    for i=1:N2
        for j=1:N2
            r=sqrt(x(i)^2+yd(j)^2);
            if r<r0
                F(j,i)=1e-5*cos((pi/2)*r/r0);
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
    
    %============================================================================================

    % Now we solve the problem
    %============================================================================================

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
    
    if test==0
        waveinfo1(ti,:)=wave(v1hat,K,Lx);   % Upper layer
        waveinfo2(ti,:)=wave(v2hat,K,Lx);   % Lower layer
    end
    if test==1
        waveinfo1(ti,:,Ti)=wave(v1hat,K,Lx);   % Upper layer
        waveinfo2(ti,:,Ti)=wave(v2hat,K,Lx);   % Lower layer
    end
    
    
end

%%

wavesmooth1=zeros(N_T,2);
% Can't average the first and last values
wavesmooth1(1,1)=waveinfo1(1,1);
wavesmooth1(1,2)=waveinfo1(1,2);
wavesmooth1(N_T,1)=waveinfo1(N_T,1);
wavesmooth1(N_T,2)=waveinfo1(N_T,2);
% Average the second and second-last values over 3 samples
wavesmooth1(2,1)=(waveinfo1(1,1)+waveinfo1(2,1)+waveinfo1(3,1))/3;
wavesmooth1(2,2)=(waveinfo1(1,2)+waveinfo1(2,2)+waveinfo1(3,2))/3;
wavesmooth1(N_T-1,1)=(waveinfo1(N_T-2,1)+waveinfo1(N_T,1)+waveinfo1(N_T,1))/3;
wavesmooth1(N_T-1,2)=(waveinfo1(N_T-2,2)+waveinfo1(N_T,2)+waveinfo1(N_T,2))/3;
% Average over all other values
for ti=3:N_T-2
    wavesmooth1(ti,1)=(waveinfo1(ti-2,1)+waveinfo1(ti-1,1)+waveinfo1(ti,1)+waveinfo1(ti+1,1)+waveinfo1(ti+2,1))/5;
    wavesmooth1(ti,2)=(waveinfo1(ti-2,2)+waveinfo1(ti-1,2)+waveinfo1(ti,2)+waveinfo1(ti+1,2)+waveinfo1(ti+2,2))/5;
end

wavesmooth2=zeros(N_T,2);
% Can't average the first and last values
wavesmooth2(1,1)=waveinfo2(1,1);
wavesmooth2(1,2)=waveinfo2(1,2);
wavesmooth2(N_T,1)=waveinfo2(N_T,1);
wavesmooth2(N_T,2)=waveinfo2(N_T,2);
% Average the second and second-last values over 3 samples
wavesmooth2(2,1)=(waveinfo2(1,1)+waveinfo2(2,1)+waveinfo2(3,1))/3;
wavesmooth2(2,2)=(waveinfo2(1,2)+waveinfo2(2,2)+waveinfo2(3,2))/3;
wavesmooth2(N_T-1,1)=(waveinfo2(N_T-2,1)+waveinfo2(N_T-1,1)+waveinfo2(N_T,1))/3;
wavesmooth2(N_T-1,2)=(waveinfo2(N_T-2,2)+waveinfo2(N_T-1,1)+waveinfo2(N_T,2))/3;
% Average over all other values
for ti=3:N_T-2
    wavesmooth2(ti,1)=(waveinfo2(ti-2,1)+waveinfo2(ti-1,1)+waveinfo2(ti,1)+waveinfo2(ti+1,1)+waveinfo2(ti+2,1))/5;
    wavesmooth2(ti,2)=(waveinfo2(ti-2,2)+waveinfo2(ti-1,2)+waveinfo2(ti,2)+waveinfo2(ti+1,2)+waveinfo2(ti+2,2))/5;
end

%%


if test==0
    % Upper Layer
    figure(1)
    plot(T,waveinfo1(:,1),'k',T,wavesmooth1(:,1),'r','LineWidth',1.2);...
        xlabel('Forcing Period (days)'); ylabel('Dominant Wavenumber - Upper Layer');...
        text(50,max(waveinfo1(:,1))-0.5,['U_{1}=',num2str(U1),],'FontSize',14);...
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/Tests/Freq/','k1_upper'],'png');
    figure(2)
    plot(T,waveinfo1(:,2),'k',T,wavesmooth1(:,2),'r','LineWidth',1.2);...
        xlabel('Forcing Period (days)'); ylabel('Second-most Dominant Wavenumber - Upper Layer');...
        text(50,max(waveinfo1(:,2))-0.5,['U_{1}=',num2str(U1),],'FontSize',14);...
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/Tests/Freq/','k2_upper'],'png');

    % Lower layer
    figure(3)
    plot(T,waveinfo2(:,1),'k',T,wavesmooth2(:,1),'r','LineWidth',1.2);...
        xlabel('Forcing Period (days)'); ylabel('Dominant Wavenumber - Lower Layer');...
        text(50,max(waveinfo2(:,1))-1.5,['U_{1}=',num2str(U1),],'FontSize',14);...
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/Tests/Freq/','k1_lower'],'png');
    figure(4)
    plot(T,waveinfo2(:,2),'k',T,wavesmooth2(:,2),'r','LineWidth',1.2);...
        xlabel('Forcing Period (days)'); ylabel('Second-most Dominant Wavenumber - Lower Layer');...
        text(50,max(waveinfo2(:,2))-1.5,['U_{1}=',num2str(U1),],'FontSize',14);...
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/Tests/Freq/','k2_lower'],'png')
end

if test==1
    figure(1)
    plot(U1range,waveinfo1(:,1,1),'g',U1range,waveinfo1(:,1,2),'b',U1range,waveinfo1(:,1,3),'r',...
        U1range,waveinfo1(:,1,4),'k','LineWidth',1.2);...
        legend('\omega=1/50','\omega=1/60','\omega=1/70','\omega=1/80');...
        xlabel('$U_{1} \,$ (m s$^{-1}$)','Interpreter','LaTex');...
        ylabel({'Dominant Wavenumber'; 'Upper Layer'}); axis([-0.16 0.32 -50 50]);...
        line([-0.16 0.32], [0 0],'color','k', 'LineStyle','--'); line([0 0],ylim,'color','k','LineStyle','--');...
        set(gca, 'XTick', [-0.16 -0.08 0 0.08 0.16 0.24 0.32]);...
        set(gcf,'units','points','position',[0,0,800,200]);...
        set(gcf,'PaperPositionMode','auto');...
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/Tests/Freq/','k1_upper_U1'],'png');
    figure(2)
    plot(U1range,waveinfo2(:,1,1),'g',U1range,waveinfo2(:,1,2),'b',U1range,waveinfo2(:,1,3),'r',...
        U1range,waveinfo2(:,1,4),'k','LineWidth',1.2);...
        legend('\omega=1/50','\omega=1/60','\omega=1/70','\omega=1/80');...
        xlabel('$U_{1} \,$ (m s$^{-1}$)','Interpreter','LaTex');...
        ylabel({'Dominant Wavenumber'; 'Lower Layer'}); axis([-0.16 0.32 -50 50]);...
        line([-0.16 0.32], [0 0],'color','k', 'LineStyle','--'); line([0 0],ylim,'color','k','LineStyle','--');...
        set(gca, 'XTick', [-0.16 -0.08 0 0.08 0.16 0.24 0.32]);...
        set(gcf,'units','points','position',[0,0,800,200]);...
        set(gcf,'PaperPositionMode','auto');...
        saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/2-Layer/Tests/Freq/','k1_lower_U1'],'png');
end
