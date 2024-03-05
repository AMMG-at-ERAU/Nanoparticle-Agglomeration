function [  Ptrue, Center, Length, Angle, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_angle( RVE)
%Generate agglomerate CNTs
%   Detailed explanation goes here
%   1) Generate all the non agglomerated, Vf = 1-Vf_agg
%   2) Randomly choose the Na CNTs, from the non agglomerated, that will be
%   used as agglomeration center/inclusions (a cnt can be chosen multiple time)
%   3) Generate the CNT in agglomerate, such that their centre is in a
%   certain raduis aroung the agglommerate center
%   PPi1 has all the CNTs without PBC as column
%   PPi has all the CNTs with PBC as column

Vpr=0;
Ptrue1=[];
Center1 = [];
Length1 = [];
Ptrue = [];
Center = [];
Angle = [];
Length = [];
Lx=RVE.size(1);
Lz=RVE.size(3);
agg_rad = RVE.agg_rad; %radius of each agglomerate center
agg_angle = RVE.agg_angle*pi/180; %convert from degree to radians
while Vpr<RVE.Vf*Lx*Lz*(1-RVE.Vf_agg) %generation of non agglomerated CNTs
    r = rand(4,1);
    % Center points (cx,cy,cz) is random
    cx = Lx*r(1,:);
    cy = 0;
    cz = Lz*r(2,:);
    phi = 2*pi*r(3,:);
    li = RVE.Wa*(-log(1-r(4,:))).^(1/RVE.Wb);
    Vp = RVE.D*li;
    % Starting and ending points
    x1 = cx-0.5*li.*(cos(phi));
    y1 = 0;
    z1 = cz-0.5*li.*(sin(phi));
    x2 = cx+0.5*li.*(cos(phi));
    y2 = 0;
    z2 = cz+0.5*li.*(sin(phi));
    pt.P = [x1 y1 z1;x2 y2 z2]; %two end points of CNT
    pt.center = [cx cy cz]; %Center of CNT
    pt.length = li; %length of CNT
    [P, center1, li1] = PBC_2D(pt,RVE);
    Ptrue = [Ptrue pt.P']; %% Matrix with all the CNT generated in the box without BCs
    Center = [Center [cx;cy;cz]]; %Matrix with all the center of CNT generated in the box without BCs
    Angle = [Angle phi];%Matrix with all the angle of CNT generated in the box without BC
    Length = [Length li];%Matrix with all the length of CNT generated in the box without BC
    Ptrue1=[Ptrue1 P];   %% Matrix with all the CNT generated in the box with BCs
    Center1 = [Center1 center1]; %Matrix with all the center of CNT generated in the box with BCs
    Length1 = [Length1 li1];%Matrix with all the length of CNT generated in the box with BC
    Vpr=Vpr+Vp;
end
N = size(Center,2); %Total number of non agglomerate CNT
Na = ceil(RVE.Vf_agg*N/(1-RVE.Vf_agg));
Nrand = datasample(1:N,3*Na); %To sample the indices of CNT. Change Na to 2*Na if error
r1 = rand(1,3*Na);
r2 = rand(1,3*Na);
x_a = sqrt(r1).*cos(2*pi*r2);
y_a = sqrt(r1).*sin(2*pi*r2);
i = 1;
while Vpr<RVE.Vf*Lx*Lz %generation of agglomerated CNTs
    % Center points (cx,cy,cz) is random
    cx = Center(1,Nrand(i)) + agg_rad*x_a(i);
    cy = 0;
    cz = Center(3,Nrand(i)) + agg_rad*y_a(i);
    if cx>Lx || cx<0 % if the center is outside of RVE
        cx = (cx>Lx)*(Lx-0.01) + (cx<=Lx)*0.01; 
    end
    if cz>Lz || cz<0 % if the center is outside of RVE
        cz = (cz>Lz)*(Lz-0.01) + (cz<=Lz)*0.01; 
    end
    r = rand(2,1);
    mi = Angle(1,Nrand(i)) - agg_angle;
    ma = Angle(1,Nrand(i)) + agg_angle;
    phi = mi + (ma - mi)*r(1,:);
    %phi = 2*pi*r(1,:);
    li = RVE.Wa*(-log(1-r(2,:))).^(1/RVE.Wb);
    Vp = RVE.D*li;
    x1 = cx-0.5*li.*(cos(phi));
    y1 = 0;
    z1 = cz-0.5*li.*(sin(phi));
    x2 = cx+0.5*li.*(cos(phi));
    y2 = 0;
    z2 = cz+0.5*li.*(sin(phi));
    % Apply Periodic BC
    pt.P = [x1 y1 z1;x2 y2 z2]; %two end points of CNT
    pt.center = [cx cy cz]; %Center of CNT
    pt.length = li; %length of CNT
    [P, center1,li1] = PBC_2D(pt,RVE);
    Ptrue = [Ptrue pt.P']; %% Matrix with all the CNT generated in the box without BCs
    Center = [Center [cx;cy;cz]]; %Matrix with all the center of CNT generated in the box without BCs
    Angle = [Angle phi];%Matrix with all the angle of CNT generated in the box without BC
    Length = [Length li];%Matrix with all the length of CNT generated in the box without BC
    Ptrue1=[Ptrue1 P];   %% Matrix with all the CNT generated in the box with BCs
    Center1 = [Center1 center1]; %Matrix with all the center of CNT generated in the box with BCs
    Length1 = [Length1 li1];%Matrix with all the length of CNT generated in the box with BC
    Vpr=Vpr+Vp;
    i = i+1;
end

end