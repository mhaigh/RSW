function f=ddx(u,dx)
% A function to find the x-derivative of a vector u

xlen=size(u,2);
ylen=size(u,1);

twodx=2*dx;

ddxu=zeros(ylen,xlen);
for j=1:ylen
    for i=2:xlen-1
        ddxu(j,i)=(u(j,i+1)-u(j,i-1))/twodx;
    end
    ddxu(j,1)=(u(j,2)-u(j,xlen))/twodx;
    ddxu(j,xlen)=(u(j,1)-u(j,xlen-1))/twodx;
end    

f=ddxu;
