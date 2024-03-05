clear all 
clc

% Load the Excel spreadsheet into MATLAB
filename = 'CNT_center25.csv';
num = readmatrix(filename);



K=num;
filename = 'CNT_center25.xlsx';
writematrix(K, filename)