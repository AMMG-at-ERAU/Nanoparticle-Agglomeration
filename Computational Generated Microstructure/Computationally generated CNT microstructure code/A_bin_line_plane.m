function [I]= A_bin_line_plane(n1,V0,P0,P1)
%plane_line_intersect computes the intersection of a plane and a segment(or
%a straight line)
% Inputs: 
%       n1: normal vector to all the Planes in that direction, 1x3 matrix 
%       V0: any point that belongs to the Planes, kx3 matrix 
%       P0: end point 1 of the segment P0P1
%       P1:  end point 2 of the segment P0P1
%
%Outputs:
%      I    is the point of interection 
k = size(V0,1); %number of planes
%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%If you have any comments or face any problems, please feel free to leave
%your comments and i will try to reply to you as fast as possible.
n = ones(k,1) * n1;
%I = zeros(k,3);
u = ones(k,1) * (P1-P0);
w = (ones(k,1) * P0) - V0;
D = sum(n.*u,2); %column vector
N = -sum(n.*w,2);%column vector

%compute the intersection parameter
sI = N ./ D;    %column vector
I = (ones(k,1) * P0) + (sI * ones(1,3)).*u;  
end