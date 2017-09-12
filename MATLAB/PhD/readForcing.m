% readForcing

clear all

loc = '~/cluster/gold3/';

cd(loc);

%T = ncread('timestats.nc','Time');
files = dir('forc*');

%nt = size(T,1);
nf = size(files);
nf = nf(1);

i = nf-1;

wind = 1;
if wind == 1
    taux = ncread(files(i).name,'taux');
    taux = transpose(squeeze(taux));
    
    surf(taux,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;
    
    
end

pause

buoy = 0;
if buoy == 1
    b_new = ncread(files(i).name,'buoy');
    b4 = transpose(squeeze(b_new(:,:,1))); 
    b_new = ncread(files(i-1).name,'buoy');
    b3 = transpose(squeeze(b_new(:,:,1)));
    b_new = ncread(files(i-2).name,'buoy');
    b2 = transpose(squeeze(b_new(:,:,1)));
    b_new = ncread(files(1).name,'buoy');
    b1 = transpose(squeeze(b_new(:,:,1)));
    
    pause

    subplot(221);
    surf(b1,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;
    subplot(222);
    surf(b2,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;
    subplot(223);
    surf(b3,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;
    subplot(224);
    surf(b4,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;

    pause
    
    subplot(131);
    plot(b(:,1));
    subplot(132);
    plot(b(:,100));
    subplot(133);
    plot(b(:,200));
    
end

