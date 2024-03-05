function [ r, ptr] = find_root( n, ptr )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
p=0;
stack = zeros(1,60000);
r = n;
while ptr(r)>=0
    p = p+1;
    stack(p) = r;
    r = ptr(r);
end
ptr(stack(1:p)) = r;
end

