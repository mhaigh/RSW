function f=Error1Lvisc(utilde,vtilde,htilde,nu,f,w,t,g,H,N,N2,dy,dx,Ffull)
% A function to find the error of the RSW code

I=complex(0,1);

u=real(utilde*cos(w*t));
v=real(vtilde*cos(w*t));
h=real(htilde*sin(w*t));

% U equation
%=======================================================================

u1=real(-w*utilde*sin(w*t));
v1=zeros(N,N2);
for j=1:N
    v1(j,:)=-f(j,1)*v(j,:);
end
h1=g.*ddx(h,dx);
visc1=nu*ddx(ddx(u,dx),dx)+nu*ddy(ddy(u,dy,'periodic'),dy,'periodic');

E1=u1+v1+h1-visc1;

% V equation
%=======================================================================

v2=real(-w*vtilde*sin(w*t));
u2=zeros(N,N2);
for j=1:N
    u2(j,:)=f(j,1)*u(j,:);
end
h2=g.*ddy(h,dy,'periodic');

E2=u2+v2+h2;

% H equation
%=======================================================================

h3=real(-w*htilde.*sin(w*t));
u3=H.*ddx(u,dx);
v3=H.*ddy(v,dy,'periodic');

E3=u3+v3+h3-Ffull*cos(w*t);

%======================================================================

e1=immse(E1,zeros(N,N2));
e2=immse(E2,zeros(N,N2));
e3=immse(E3,zeros(N,N2));

f=[e1,e2,e3];

end