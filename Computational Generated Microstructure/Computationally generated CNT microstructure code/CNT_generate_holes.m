function [Ptrue, center, L, hole] = CNT_generate_holes(RVE,Ptrue1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Vpr = 0;
hole = [];
% Define area where hole can be generated
mi = 0.5*RVE.hole_D;
max = RVE.size(1) - 0.5*RVE.hole_D;
maz = RVE.size(3) - 0.5*RVE.hole_D;
while Vpr<RVE.Vf_holes*RVE.size(1)*RVE.size(3)
    % Generate circle
    r = rand(2,1);
    cx = mi + (max - mi)*r(1,:);
    cy = 0;
    cz = mi + (maz - mi)*r(2,:);
    hole = [hole;[cx,cy,cz]];
    Pnew = [];
    for i = 1 : 0.5*size(Ptrue1,1)
        P0 = Ptrue1(2*i-1:2*i,:);
        P1 = CNT_after_hole(RVE,P0,[cx,cy,cz]);
        Pnew = [Pnew;P1];
    end
    Ptrue1 = Pnew;
    
    Vpr = Vpr + 0.25*pi*RVE.hole_D^2;
end

Ptrue = Ptrue1';
n = size(Ptrue,2);
L = Ptrue(:,1:2:n-1)-Ptrue(:,2:2:n);
L = sqrt(sum(L.*L,1));
center = 0.5*(Ptrue(:,1:2:n-1)+Ptrue(:,2:2:n));
Ptrue = Ptrue';
center = center';

end

