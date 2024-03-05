function [ R ] = find_leafs( n, ptr )
%Find the leafs of any main roots without recursive
%   Detailed explanation goes here

R = zeros(1,length(ptr));
a = (ptr ==n);
p = 0;
while ~isempty(find(a,1))
    b = find(a);
    c = length(b);
    R(p+1:p+c) = b;
    p = p + c;
    a(:)=0;
    for i = 1:c
        R1 = (ptr == b(i));
        a = a + R1;
    end
end
R = R(R>0);
end

