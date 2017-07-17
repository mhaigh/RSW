function f=ddy(u,dy,option)
% A function to find the x-derivative of a vector u

if nargin<3
    option=1;
end

xlen=size(u,2);
ylen=size(u,1);

twody=2*dy;

ddyu=zeros(ylen,xlen);
if option=='periodic'
    for i=1:xlen
        for j=2:ylen-1
            ddyu(j,i)=(u(j+1,i)-u(j-1,i))/twody;
        end
        ddyu(1,i)=(u(2,i)-u(ylen,i))/twody;
        ddyu(ylen,i)=(u(1,i)-u(ylen-1,i))/twody;
    end
else
    for i=1:xlen
        for j=2:ylen-1
            ddyu(j,i)=(u(j+1,i)-u(j-1,i))/twody;
        end
        ddyu(1,i)=(u(2,i)-u(1,i))/dy;
        ddyu(ylen,i)=(u(ylen,i)-u(ylen-1,i))/dy;
    end
end

    

f=ddyu;