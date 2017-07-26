% massConservation

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

mass1 = [0];
mass2 = [0];
mass3 = [0];
count = 0;
%for i = 3:3
for i = nf-5:nf-1
    disp(i);
    h_new = ncread(files(i).name,'h');
    h1_new = h_new(:,:,1,:);
    h2_new = h_new(:,:,2,:);
    h3_new = h_new(:,:,3,:);
    nn = size(h1_new,4);    
    for ti = 1:nn
        h1(:,:,count+ti) = transpose(h1_new(:,:,ti));
        h2(:,:,count+ti) = transpose(h2_new(:,:,ti));
        h3(:,:,count+ti) = transpose(h3_new(:,:,ti));
        mass1(count+ti) = sum(sum(h1(:,:,count+ti)));
        mass2(count+ti) = sum(sum(h2(:,:,count+ti)));
        mass3(count+ti) = sum(sum(h3(:,:,count+ti)));
    end
    count = count + nn;
end
nt = size(h1,3);

subplot(141); plot(mass1);
subplot(142); plot(mass2);
subplot(143); plot(mass3);
subplot(144); plot(mass1+mass2+mass3);



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


%%

mass = load('mass_time');
mass = mass(:,2);

plot(mass); xlabel('days'); ylabel('Mass');...
    saveas(gcf,['~/Documents/GulfStream/GOLD/Images/','mass'],'png');




