% readForcing

clear all

loc = '~/cluster/gold5/';

cd(loc);

%T = ncread('timestats.nc','Time');
files = dir('forc*');

%nt = size(T,1);
nf = size(files);
nf = nf(1);

x = ncread(files(1).name,'xq');
y = ncread(files(1).name,'yh');

nx = size(x,1);
ny = size(y,1);

i = 1%nf-2;

wind = 1;
if wind == 1
    taux = ncread(files(i).name,'taux');
    taux = transpose(squeeze(taux(:,:,1)));
    
    surf(taux,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;
    
    pause
        
    tauy = ncread(files(i).name,'tauy');
    tauy = transpose(squeeze(tauy(:,:,1)));
    
    surf(tauy,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;
    
    pause
    
    N_vec = 16;
    taux_vec = zeros(N_vec,N_vec);
    tauy_vec = zeros(N_vec,N_vec);
    i_set = 1:nx/N_vec:nx;
    j_set = 1:ny/N_vec:ny;
    for i=1:N_vec
        for j=1:N_vec
            taux_vec(j,i) = taux(j_set(j),i_set(i));
            tauy_vec(j,i) = tauy(j_set(j),i_set(i));
        end
    end
    [x_vec,y_vec] = meshgrid(1:nx/N_vec:nx,1:ny/N_vec:ny);
    quiver(x_vec,y_vec,taux_vec,tauy_vec);
    
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

