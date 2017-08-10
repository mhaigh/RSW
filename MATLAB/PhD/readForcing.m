% readForcing

clear all

loc = '~/cluster/gold/';

cd(loc);

%T = ncread('timestats.nc','Time');
files = dir('forc*');

%nt = size(T,1);
nf = size(files);
nf = nf(1);


i = nf-1;

b_new = ncread(files(i).name,'buoy');
b4 = transpose(squeeze(b_new(:,:,1))); 
b_new = ncread(files(i-1).name,'buoy');
b3 = transpose(squeeze(b_new(:,:,1)));
b_new = ncread(files(i-2).name,'buoy');
b2 = transpose(squeeze(b_new(:,:,1)));
b_new = ncread(files(1).name,'buoy');
b1 = transpose(squeeze(b_new(:,:,1)));

%%
subplot(221);
surf(b1,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;
subplot(222);
surf(b2,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;
subplot(223);
surf(b3,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;
subplot(224);
surf(b4,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;

pause

%%

subplot(131);
plot(b(:,1));
subplot(132);
plot(b(:,100));
subplot(133);
plot(b(:,200));

pause

tau_new = ncread(files(i).name,'taux');
tau = transpose(squeeze(tau_new(:,:,1)));

tau_max = max(max(tau));
tau_min = min(min(tau));

surf(tau,'edgecolor','none'); view(0,90); colorbar; colormap(jet); title('wind'); caxis([tau_min tau_max]);

