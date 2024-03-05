clear all 
clc 
close all

filename = 'CNT_center1.xlsx';
filename_1 = 'endpoints_1_allocation.xlsx';

num = xlsread(filename);
j = num(:, :);
N = size(j,1);
x = j(:,1);
y = j(:,2);
z = j(:,3);

num_1 = xlsread(filename_1);
l = num_1(:, :);
M = size(l,1);
x_1 = l(:,1);
y_1 = l(:,2);
z_1 = l(:,3);
x_2 = l(:,4);
y_2 = l(:,5);
z_2 = l(:,6);


ans(:,1) = [x];
ans(:,2) = [y];
ans(:,3) = [z];
ans(:,5) = [x_1];
ans(:,6) = [y_1];
ans(:,7) = [z_1];
ans(:,8) = [x_2];
ans(:,9) = [y_2];
ans(:,10) = [z_2];

for i=1:N;
ans(i,13) = sqrt(x(i)^2+y(i)^2);
end

K=ans;
writematrix(K, filename)