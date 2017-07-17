% 1L Wavenumber analysis
% The code is set up to test how the magnitude of the solution changes
% as with r0. It calculates the magnitude for a range of omega (w)

%clear all
close all

wi=1;           % Need to cycle this parameter through 1-4

N_r0=100;

RES=4;          % 1-8 increasing resolution
OMEGA=2*pi/(3600*24)*[1/50,1/60,1/70,1/80,1/90,1/100,1/110,1/120];   % Forcing frequencies

Nsel=[18,34,66,130,258,514,1026,2050];  % Resolution, 2^n +2. Choices 7,8 take a while!

w=OMEGA(wi);
ts=(1)*2*pi/w;       % Time at which to plot the snapshot (if no EDDYFORCING, otherwies it's overwritten)

N=Nsel(RES);     % Number of gridpoints
N2=N-2;
Ly=3000000;     % Size of the basin
Lx=3000000;
dy=Ly/(N-1);     % Distance between gridpoints
dx=Lx/(N2-1);

r0min=Lx/82;
r0max=Lx/12;
r0range=linspace(r0min,r0max,N_r0);

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


%Define the RHS of the ODE - the plunger forcing in k-y space
%==========================================================================
% F is the forcing term on the RHS of the original continuity equation.
% Fhat is the Fourier transform of F.
% RHS is a combination of Fhat and its derivative that's on the RHS of the
% vhat differential equation.
% Each of these forcing terms are initially defined on the smaller domain before being
% extended to include the two dead points in the y-direction, e.g. Ffull.

for ri=1:N_r0
    r0=r0range(ri);
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
    %F=F-sum(sum(F))/(N2*N2);
    Fhat=fftshift(F,2);
    %Fhat=fft(Fhat,[],2);  % Need to sort this: F is even so Fhat should be real. 
    Fhat=N2*ifft(Fhat,[],2,'symmetric');
    %Fhat=conj(Fhat);      % Also need to check if this conj is needed in general. 
    RHS1=zeros(N2,N2);
    for i=1:N2
        for j=1:N2
            RHS1(j,i)=-f(j+1,1)*K(1,i)*Fhat(j,i)/(w*H);
        end
    end
    RHS2=ddy(Fhat,dy)/H;
        
    RHS=RHS2+RHS1;    

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
    c=zeros(1,N2);

    % Here we solve the problem
    A=zeros(N2,N2);
    i=1;
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
    vhatlim=max(max(abs(vhat)));

    vtilde=ifft(vhat,[],2);   % parameter 2 ensures ifft over rows, rather than columns
    %vtilde=conj(vtilde);
    vtilde=ifftshift(vtilde,2);           % shifts zero-frequency back to the centre
    v=real(exp(I*w*ts)*vtilde);
    
    vmag(ri,wi)=sum(sum(abs(v)))/N2^2;

end

%%

plot(r0range,vmag(:,1),'g',r0range,vmag(:,2),'b',r0range,vmag(:,3),'r',r0range,vmag(:,4),'k','LineWidth',1.2);...
    legend('\omega=1/50','\omega=1/60','\omega=1/70','\omega=1/80','Location','Northwest');
    xlabel('r_{0} (m)'); ylabel('v abs average');...
    axis([r0min r0max 0.0 1.1*max(max(vmag))]);...
    %set(gca, 'XTick', [Lx/64 Lx/16],'XTickLabels'...
    %,{'L_{x}/64=46.875','L_{x}/16=187.5'});...
    %set(gca,'XTickLabelRotation',45);...
    line([Lx/16 Lx/16],ylim,'color','k','LineStyle','--'); text(Lx/16-4000,0.48e-4,'L_{x}/16','rotation',90);...
    line([Lx/20 Lx/20],ylim,'color','k','LineStyle','--'); text(Lx/20-4000,0.48e-4,'L_{x}/20','rotation',90);...
    line([Lx/32 Lx/32],ylim,'color','k','LineStyle','--'); text(Lx/32-4000,0.42e-4,'L_{x}/32','rotation',90);...
    line([Lx/64 Lx/64],ylim,'color','k','LineStyle','--'); text(Lx/64-4000,0.48e-4,'L_{x}/64','rotation',90);...
    set(gcf,'units','points','position',[0,0,660,220]);...
    set(gcf,'PaperPositionMode','auto');...
    saveas(gcf,['~/Documents/GulfStream/Interim_Report/Images/1-Layer/Tests/r0/','vav_r0_1L'],'png');