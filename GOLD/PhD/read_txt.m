cd '~/Documents/GulfStream/Code/PYTHON/TESTS/';

fileID = fopen('results.txt','r');

formatSpec = '%f';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);

plot(A(2,:));
