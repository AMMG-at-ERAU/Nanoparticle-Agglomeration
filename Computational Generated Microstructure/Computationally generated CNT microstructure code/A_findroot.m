function [ R ] = A_findroot( n, ptr )
%Find the main root of any leafs
%   Detailed explanation goes here
if ptr(n) < 0
    R = n;
else
    R = A_findroot(ptr(n),ptr);
end
end

