function f=Error1LQG(phitilde,beta,w,S,ts,N,N2,dy,dx,Ffull)
% A function to find the error of the RSW code

I=complex(0,1);

phi=real(exp(I*w*ts)*phitilde);

phi_t=real(I*w*exp(I*w*ts)*phitilde);

% Calculate the del^2 of time-derivative of phi
del2phi_t=ddx(ddx(phi_t,dx),dx)+ddy(ddy(phi_t,dy),dy);

% The full QG equation
Total=del2phi_t-S*phi_t+beta*ddx(phi,dx)-real(Ffull*exp(I*w*ts));

e=immse(Total,zeros(N,N2));

f=e;

end
