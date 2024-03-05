function poursave_new( V, C, V1, C1, R1, P, K, d, j, s, u, t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T.Vector = V;
T.Center = C;
T.Vector1 = V1;
T.Center1 = C1;
T.Resistor1 = R1;
T.PA = P;
T.K1 = K;
T.dCNT = d;
T.ja = j;
T.StorE = s;
T.RVE = u;
save(t,'-struct','T');
end

