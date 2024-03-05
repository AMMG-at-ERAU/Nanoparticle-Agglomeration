function poursave_contact( dCNT, ja, t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T.dCNT = dCNT;
T.ja = ja;
save(t,'-struct','T');
end