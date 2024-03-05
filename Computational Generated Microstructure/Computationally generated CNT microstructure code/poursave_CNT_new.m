function poursave_CNT_new( dCNT, ja, Ptrue, Center, Length, UU, Ptrue1, Center1, Length1, store, RVE, t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T.dCNT = dCNT;
T.ja = ja;
T.Ptrue = Ptrue;
T.Center = Center;
T.Length = Length;
T.Ptrue1 = Ptrue1;
T.Center1 = Center1;
T.Length1 = Length1;
T.UU = UU;
T.store = store;
T.RVE = RVE;
save(t,'-struct','T');
end

