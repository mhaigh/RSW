function F=wave(vhat,K,Lx)
% Calculates the dominant wavenumber of a vector in Fourier-physical space

N=size(vhat,2);

vabs=abs(vhat);
vsum=sum(vabs,1);

vsorted=sort(vsum,'descend');

k1mag=vsorted(1);
k2mag=vsorted(2);

if k2mag==k1mag
    k2mag=vsorted(3);
end

kratio=k1mag/k2mag;
    

for i=1:N
    if k1mag==vsum(i)
        k1=K(i);
        i;
    end
    if k2mag==vsum(i)
        k2=K(i);
        j=i;
        j;
    end
end

k1=k1*Lx/(2*pi);
k2=k2*Lx/(2*pi);

F=[k1,k2,kratio];



    

