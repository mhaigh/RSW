function [umap, ulim] = makeMap(u,Nb)

map_tmp = zeros(10,3);
map_tmp(10,:) = [144,12,15];
map_tmp(9,:) = [178,24,43];
map_tmp(8,:) = [214,96,77];
map_tmp(7,:) = [244,165,130];
map_tmp(6,:) = [253,219,199];
map_tmp(5,:) = [209,229,240];
map_tmp(4,:) = [146,197,222];
map_tmp(3,:) = [67,147,195];
map_tmp(2,:) = [33,102,172];
map_tmp(1,:) = [10,62,138];
map_tmp = map_tmp / 255;

ulim = max(max(max(abs(u))));

umap = zeros(Nb,3);
for i = 1:3
    umap(:,i) = interp1(linspace(-ulim,ulim,10),map_tmp(:,i),linspace(-ulim,ulim,Nb));
end