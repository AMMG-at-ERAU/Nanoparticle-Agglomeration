function [ r, ptr] = A_find_root2( n, ptr )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
r = n;
s = n;
while ptr(r)>=0
    ptr(s) = ptr(r);
    s = r;
    r = ptr(r);
end

end
