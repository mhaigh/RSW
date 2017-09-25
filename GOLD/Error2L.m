function f=Error2L(u1tilde,v1tilde,eta0tilde,u2tilde,v2tilde,eta1tilde,f,w,t,g,H1,H2,rho1prime,rho2prime,U1,U2,N,N2,dy,dx,Ffull)
% A function to find the error of the RSW code

I=complex(0,1);

u1=real(exp(I*w*t)*u1tilde);
u2=real(exp(I*w*t)*u2tilde);
v1=real(exp(I*w*t)*v1tilde);
v2=real(exp(I*w*t)*v2tilde);
eta0=real(exp(I*w*t)*eta0tilde);
eta1=real(exp(I*w*t)*eta1tilde);

h1=eta0-eta1;

% u1 equation
%=========================================

u11=real(I*w*exp(I*w*t)*u1tilde)+U1*ddx(u1,dx);
v11=zeros(N,N2);
for j=1:N
    v11(j,:)=-f(j,1)*v1(j,:);
end
eta01=g*ddx(eta0,dx);

E1=u11+v11+eta01;

% v1 equation
%=========================================

v11=real(I*w*exp(I*w*t)*v1tilde)+U1*ddx(v1,dx);
for j=1:N
    u11(j,:)=f(j,1)*u1(j,:);
end
eta01=g*ddy(eta0,dy);

E2=u11+v11+eta01;

% h1 equation
%=========================================

h11=real(I*w*exp(I*w*t)*(eta0tilde-eta1tilde))+U1*ddx(h1,dx);
u11=H1*ddx(u1,dx);
v11=H1*ddy(v1,dy);

E3=h11+u11+v11-Ffull*cos(w*t);

% u2 equation
%=========================================

u21=real(I*w*exp(I*w*t)*u2tilde)+U2*ddx(u2,dx);
v21=zeros(N,N2);
for j=1:N
    v21(j,:)=-f(j,1)*v2(j,:);
end
eta11=g*rho1prime*ddx(eta0,dx)+g*rho2prime*ddx(eta1,dx);

E4=u21+v21+eta11;

% v2 equation
%=========================================

v21=real(I*w*exp(I*w*t)*v2tilde)+U2*ddx(v2,dx);
for j=1:N
    u21(j,:)=f(j,1)*u2(j,:);
end
eta11=g*rho1prime*ddy(eta0,dy)+g*rho2prime*ddy(eta1,dy);

E5=u21+v21+eta11;

% h2 equation
%=========================================

eta11=real(I*w*exp(I*w*t)*eta1tilde)+U1*ddx(h1,dx);
u21=H2*ddx(u2,dx);
v21=H2*ddy(v2,dy);

E6=eta11+u21+v21;

%==========================================

e1=sqrt(immse(E1,zeros(N,N2)));
e2=immse(E2,zeros(N,N2));
e3=immse(E3,zeros(N,N2));
e4=immse(E4,zeros(N,N2));
e5=immse(E5,zeros(N,N2));
e6=immse(E6,zeros(N,N2));


f=[e1,e2,e3,e4,e5,e6];

end
