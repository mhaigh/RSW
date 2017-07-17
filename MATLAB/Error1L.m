function f=Error1L(utilde,vtilde,htilde,f,w,t,g,H,N,N2,dy,dx,Ffull)
% A function to find the error of the RSW code

I=complex(0,1);

u=real(utilde*exp(I*w*t));
v=real(vtilde*exp(I*w*t));
h=real(htilde*exp(I*w*t));


% U equation
%============================

u1=real(I*w*utilde*exp(I*w*t));
v1=zeros(N,N2);
for j=1:N
    v1(j,:)=-f(j,1)*v(j,:);
end
h1=g*ddx(h,dx);

E1=u1+v1+h1;

% V equation
%===============================

v2=real(I*w*vtilde*exp(I*w*t));
u2=zeros(N,N2);
for j=1:N
    u2(j,:)=f(j,1)*u(j,:);
end
h2=g*ddy(h,dy);

E2=u2+v2+h2;

% H equation
%===============================

h3=real(I*w*htilde*exp(I*w*t));
u3=H*ddx(u,dx);
v3=H*ddy(v,dy);

E3=u3+v3+h3-Ffull*cos(w*t);

%================================

e1=immse(E1,zeros(N,N2));
e2=immse(E2,zeros(N,N2));
e3=immse(E3,zeros(N,N2));

%f=zeros(N,N2,3);
%f(:,:,1)=E1;
%f(:,:,2)=E2;
%f(:,:,3)=E3;

f=[e1,e2,e3];

end