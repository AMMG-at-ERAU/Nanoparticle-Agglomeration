function [ R ] = A_findleafs( n, ptr )
%Find the leafs of any main roots
%   Detailed explanation goes here
R = [];
a = (ptr ==n);
if ~isempty(find(a,1))
    b = find(a);
    R = [R b];
    for i = 1:length(b)
        R1 = A_findleafs(b(i),ptr);
        R = [R R1];
    end
end

end

