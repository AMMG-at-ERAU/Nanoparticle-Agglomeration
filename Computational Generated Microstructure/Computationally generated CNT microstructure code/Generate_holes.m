function [hole] = Generate_holes(RVE, Vf_holes, Radius_holes)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
format long
Vpr = 0;
hole = [];
Vf = Vf_holes;
Radius = Radius_holes;
% Define area where hole can be generated
mi = Radius;
max = RVE.size(1) - Radius;
maz = RVE.size(3) - Radius;

% First hole
    r = rand(2,1);
    cx = mi + (max - mi)*r(1,:);
    cy = 0;
    cz = mi + (maz - mi)*r(2,:);
    hole = [hole;[cx,cy,cz]];
    Vpr = Vpr + (pi*Radius^2)/(RVE.size(1)*RVE.size(3));

while Vpr<Vf
    k = 1;
    % Generate circle
    r = rand(2,1);
    cx = mi + (max - mi)*r(1,:);
    cy = 0;
    cz = mi + (maz - mi)*r(2,:);
    
    dist1 = Dist_point2holes(hole,[cx,cy,cz]); %distance with all holes center
    dist1 = (dist1<2*Radius);
    dist1 = sum(dist1);
    if (dist1>=1)
        k = 0;
        cx = [];
        cy = [];
        cz = [];
    end
    hole = [hole;[cx,cy,cz]];
    Vpr = Vpr + (pi*k*Radius^2)/(RVE.size(1)*RVE.size(3));
end


end

