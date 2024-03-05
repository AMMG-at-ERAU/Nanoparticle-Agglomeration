function poursave_CNT( dCNT, ja, Ptrue1, Length1, RVE, t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T.dCNT = dCNT;
T.ja = ja;
T.Ptrue1 = Ptrue1;
T.Length1 = Length1;
T.RVE = RVE;
save(t,'-struct','T');
end

