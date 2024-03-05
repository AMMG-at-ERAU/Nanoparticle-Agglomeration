clear all
clc

filename = 'CNTdata.xlsx';
num = xlsread(filename);
A = num(:, [1 2]);
N = size(A,1);

c = nchoosek(1:size(A, 1),2);

B = zeros( length(c), size(A,2));

for ii = 1: size(c,1)
   [r1, r2] = deal( c(ii,1), c(ii,2) );
   
   B(ii,:) = A(r1,:) - A(r2,:);
   D(ii,:) = sqrt(B(ii,1)^2 + B(ii,2)^2);
end

% filename = 'CNTdistance.xlsx';
% writematrix(D, filename)