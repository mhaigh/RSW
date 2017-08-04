% readForcing

clear all

loc = '~/cluster/gold/';

cd(loc);

%T = ncread('timestats.nc','Time');
files = dir('forc*');

%nt = size(T,1);
nf = size(files);
nf = nf(1);


i = nf-2;

b_new = ncread(files(i).name,'buoy')
b = transpose(squeeze(b_new(:,:,1))); 

%%
surf(b,'edgecolor','none'); view(0,90); colorbar; colormap(jet); axis image;

pause

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

