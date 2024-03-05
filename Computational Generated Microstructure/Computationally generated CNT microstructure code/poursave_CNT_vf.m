function poursave_CNT_vf( Ptrue1, Center1, Length1, t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T.Ptrue1 = Ptrue1;
T.Center1 = Center1;
T.Length1 = Length1;
save(t,'-struct','T');
end