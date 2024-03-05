function poursave_bis( a,b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if exist(a,'file')==0
save(a,'-struct','b');
end
end

