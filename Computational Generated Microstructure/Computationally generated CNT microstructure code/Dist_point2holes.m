function [dst] = Dist_point2holes(R1,P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Distance from a tunnelling point to all holes center
% R1 [x y z;...;x y z] : cooerdinate of holes centers
% P [x y z] : coordinates of tunneling point
format long
n = size(R1,1); % number holes
u = ones(n,1)*P - R1; %nx3 matrix
dst = sqrt(sum(u.*u,2));
end

