% visualise
clear all

loc = '~/cluster/gold6/';

cd(loc);

%T = ncread('timestats.nc','Time');
files = dir('prog__*');
x = ncread(files(1).name,'xq');
y = ncread(files(1).name,'yq');

%nt = size(T,1);
nx = size(x,1);
ny = size(y,1);
nf = size(files);
nf = nf(1);

f0 = 0.44e-4;
beta = 2e-11;
f = f0 + beta * y;

RUN = 1;    % Choose to define PV from a completed RUN (1) or a RUN in progress (2)
if RUN == 1
    PV_av = zeros(nx,ny);
    count = 0;
    %for i = 3:3
    for i = nf-20:nf-1
        disp(i);
        PVnew = ncread(files(i).name,'PV');
        unew = ncread(files(i).name,'u');
        PVnew = PVnew(:,:,1,:);
        unew = unew(:,:,1,:);
        nn = size(PVnew,4);    
        for ti = 1:nn
            PV1(:,:,count+ti)=transpose(PVnew(:,:,ti));
            u1(:,:,count+ti)=transpose(unew(:,:,ti));
        end
        count = count + nn;
    end
    nt = size(PV1,3);
end

% Use these lines to define PV from output in progress
if RUN == 2
    PV = ncread(strcat(loc,'prog__0005_002.nc'),'PV');
    PV1 = PV(:,:,1,:);
    h = ncread(strcat(loc,'prog__0005_002.nc'),'h');
    h1 = h(:,:,1,:);
    nt = 73;
    for ti = 1:nt
        PV1(:,:,ti) = transpose(PV1(:,:,ti));
        h1(:,:,ti) = transpose(h1(:,:,ti));
    end
end

%%

% Calcualte the PV average, from some time t > 0
PV_av = zeros(ny,nx);
n0 = 1;
for ti = n0:nt
    PV_av = PV_av + PV1(:,:,ti);
end
PV_av = PV_av / (nt - n0 + 1);

%%

for j = 1:nx
    disp(j);
    for i = 1:ny
        for ti = 1:nt
            PV1(j,i,ti) = PV1(j,i,ti) - f(j)/4000;
        end
    end
end

%%

qlim1 = min(min(PV1(:,:,nt)));
qlim2 = max(max(PV1(:,:,nt)));

surf(PV1(:,:,nt),'edgecolor','none'); view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
    caxis([qlim1 qlim2]);...
    saveas(gcf,['~/Documents/GulfStream/GOLD/Images/','PV_snapshot'],'png');
pause

surf(u1(:,:,nt),'edgecolor','none'); view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
        saveas(gcf,['~/Documents/GulfStream/GOLD/Images/','u_snapshot'],'png');
pause

qlim1 = min(min(PV_av));
qlim2 = max(max(PV_av));

surf(PV_av(:,:),'edgecolor','none');view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
    caxis([qlim1 qlim2]);...
    saveas(gcf,['~/Documents/GulfStream/GOLD/Images/','PV_av'],'png');
pause

qlim1 = min(min(min(PV1)));
qlim2 = max(max(max(PV1)));

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

pause 

%%
% 
h = ncread(files(nf-1).name,'h');

h1 = h(:,:,1,:);
h2 = h(:,:,2,:);
h3 = h(:,:,3,:);

nt = size(h1,4);
h1lim1 = min(min(min(h1)));
h1lim2 = max(max(max(h1)));
h2lim1 = min(min(min(h2))); 
h2lim2 = max(max(max(h2)));
h3lim1 = min(min(min(h3)));
h3lim2 = max(max(max(h3)));


surf(transpose(h1(:,:,1,nt)),'edgecolor','none'); view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
    saveas(gcf,['~/Documents/GulfStream/GOLD/Images/','h_snapshot'],'png');
pause

for ii=1:nt
    ii
    figure(1)
    surf(h1(:,:,ii),'edgecolor','none'); shading interp; view(0,90); colorbar; colormap(jet); axis image;...
        caxis([h1lim1 h1lim2]);
    mov(ii)=getframe(1);
    im=frame2im(mov(ii));
    [imind,cm] = rgb2ind(im,256);
end

%%
% ENERGY

% files = dir('energy__*');
% nf = size(files);
% count = 0;
% %for i = 20:nf
% for i = 1:nf-1
%     disp(i);
%     KEnew = ncread(files(i).name,'KE');
%     KEnew = KEnew(:,:,1,1);
%     nn = size(KEnew,4);    
%     for ti = 1:nn
%         KE1(count+ti+1) =+ sum(sum(KEnew));
%     end
%     count = count + nn;
% end
% nt = size(KE1,3);

TE = load('energy_d_mass');
TE = TE(:,2);

plot(TE); xlabel('days'); ylabel('Total Energy');...
    saveas(gcf,['~/Documents/GulfStream/GOLD/Images/','TE'],'png');
    

