% massConservation

clear all

loc = '~/cluster/gold3/';

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
for i = 1:1
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

pause



%%

mass = load('mass_time');
mass = mass(:,2);

plot(mass); xlabel('days'); ylabel('Mass');...
    saveas(gcf,['~/Documents/GulfStream/GOLD/Images/','mass'],'png');




