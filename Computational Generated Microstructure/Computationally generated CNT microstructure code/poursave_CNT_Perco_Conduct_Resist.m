function poursave_CNT_Perco_Conduct_Resist( Peroclation, conductivity, resistivity, cond, resist, t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T.Peroclation = Peroclation;
T.conductivity = conductivity;
T.resistivity = resistivity;
T.cond = cond;
T.resist = resist;
save(t,'-struct','T');
end