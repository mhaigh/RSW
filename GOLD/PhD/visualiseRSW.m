% VisualiseRSW
clear all
close all

path = '/home/mike/Documents/GulfStream/Code/PYTHON/1L/';
cd(path);

% Read the necessary data
u_full = ncread('RSW1L.nc','u');
v = ncread('RSW1L.nc','v');
eta = ncread('RSW1L.nc','eta');
x = ncread('RSW1L.nc','x');
y = ncread('RSW1L.nc','y');

% Extract a time-snapshot of the solution
ts = 10;
u=transpose(squeeze(u_full(ts,:,:)));
v=transpose(squeeze(v(ts,:,:)));
eta=transpose(squeeze(eta(ts,:,:)));

% Define the colormaps for each plot
Nb = 201;
[umap, ulim] = makeMap(u,Nb);
[vmap, vlim] = makeMap(v,Nb);
[etamap, etalim] = makeMap(eta,Nb);

path = '/home/mike/Documents/GulfStream/Code/IMAGES/1L/';
cd(path);

surf(x,y,u,'edgecolor','none'); view(0,90); colormap(umap); colorbar(); caxis([-ulim, ulim]); shading interp; axis square;...
    axis([-1/2 1/2 -1/2 1/2]);  set(gca,'LooseInset',get(gca,'TightInset'));
pause
surf(x,y,v,'edgecolor','none'); view(0,90); colormap(vmap); colorbar(); caxis([-vlim, vlim]); shading interp; axis square;...
   axis([-1/2 1/2 -1/2 1/2]);  set(gca,'LooseInset',get(gca,'TightInset'));
pause
surf(x,y,eta,'edgecolor','none'); view(0,90); colormap(etamap); colorbar(); caxis([-etalim, etalim]); shading interp; axis square;...
   axis([-1/2 1/2 -1/2 1/2]); set(gca,'LooseInset',get(gca,'TightInset'));
pause
 
 
% movie
movie = 0;
if movie == 1
[umap, ulim] = makeMap(u_full,Nb);
nt = size(u_full,1);    
    for ii=1:nt
        ii
        u = transpose(squeeze(u_full(ii,:,:)));
        figure(1)
        surf(u,'edgecolor','none'); shading interp; view(0,90); colorbar; colormap(umap); axis image;...
            caxis([-ulim ulim]);
        
        mov(ii)=getframe(1);
        im=frame2im(mov(ii));
        [imind,cm] = rgb2ind(im,256);

        if ii == 1

            %imwrite(imind,cm,'/home/mike/Documents/GulfStream/Code/IMAGES/GOLD_PV.gif', 'Loopcount',inf);

        else

            %imwrite(imind,cm,'/home/mike/Documents/GulfStream/Code/IMAGES/GOLD_PV.gif','WriteMode','append');
        end

    end
end

