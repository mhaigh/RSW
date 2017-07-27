% readForcing

clear all

loc = '~/cluster/gold4/';

cd(loc);

%T = ncread('timestats.nc','Time');
files = dir('forc*');

%nt = size(T,1);
nf = size(files);
nf = nf(1);


i = nf-1;

b_new = ncread(files(i).name,'buoy');
b = transpose(squeeze(b_new(:,:,1))); 

surf(b,'edgecolor','none'); view(0,90); colorbar; colormap(jet);



movie = 1;
if movie == 1
    
    for ii=1:nt
        ii
        figure(1)
        surf(PV1(:,:,ii),'edgecolor','none'); shading interp; view(0,90); colorbar; colormap(jet); axis image;...
            caxis([qlim1 qlim2]);
        
        mov(ii)=getframe(1);
        im=frame2im(mov(ii));
        [imind,cm] = rgb2ind(im,256);

        if ii == 1

            imwrite(imind,cm,'/home/mike/Documents/GulfStream/RSW/IMAGES/GOLD_PV.gif', 'Loopcount',inf);

        else

            imwrite(imind,cm,'/home/mike/Documents/GulfStream/RSW/IMAGES/GOLD_PV.gif','WriteMode','append');
        end

    end
end