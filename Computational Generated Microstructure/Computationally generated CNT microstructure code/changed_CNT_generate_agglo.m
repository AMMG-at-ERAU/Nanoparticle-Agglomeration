function [ Angle1, Center1, Length1, PPi1, PPi ] = changed_CNT_generate_agglo( Wa, Wb, D, Lx, Lz, Vf, agg_rad, Vf_agg)
%Generate agglomerate CNTs
%   Detailed explanation goes here
%   1) Generate all the non agglomerated, Vf = 1-Vf_agg
%   2) Randomly choose the Na CNTs, from the non agglomerated, that will be
%   used as agglomeration center/inclusions (a cnt can be chosen multiple time)
%   3) Generate the CNT in agglomerate, such that their centre is in a
%   certain raduis aroung the agglommerate center
%   PPi1 has all the CNTs without PBC as column
%   PPi has all the CNTs with PBC as column
RVE.size = [Lx,0,Lz]; 
RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
Vpr=0;
PPi=[];
PPi1 = [];
Center1 = [];
Angle1 = [];
Length1 = [];
while Vpr<Vf*Lx*Lz*(1-Vf_agg) %generation of non agglomerated CNTs
    r = rand(4,1);
    % Center points (cx,cy,cz) is random
    cx = Lx*r(1,:);
    cy = 0;
    cz = Lz*r(2,:);
    phi = 2*pi*r(3,:);
    li = Wa*(-log(1-r(4,:))).^(1/Wb);
    Vp = D*li;
    % Starting and ending points
    x1 = cx-0.5*li.*(cos(phi));
    y1 = 0;
    z1 = cz-0.5*li.*(sin(phi));
    x2 = cx+0.5*li.*(cos(phi));
    y2 = 0;
    z2 = cz+0.5*li.*(sin(phi));
    pt= [x1 y1 z1;x2 y2 z2]; %two end points of CNT
    P = PBC_2D(pt,RVE);
    PPi1 = [PPi1 pt']; %% Matrix with all the CNT generated in the box without BCs
    Center1 = [Center1 [cx;cy;cz]]; %Matrix with all the center of CNT generated in the box without BCs
    Angle1 = [Angle1 phi];%Matrix with all the angle of CNT generated in the box without BC
    Length1 = [Length1 li];%Matrix with all the length of CNT generated in the box without BC
    PPi=[PPi P];   %% Matrix with all the CNT generated in the box with BCs
    Vpr=Vpr+Vp;
end
N = size(Center1,2); %Total number of non agglomerate CNT
Na = ceil(Vf_agg*N/(1-Vf_agg));
Nrand = datasample(1:N,3*Na); %To sample the indices of CNT. Change Na to 2*Na if error
r1 = rand(1,3*Na);
r2 = rand(1,3*Na);
x_a = sqrt(r1).*cos(2*pi*r2);
y_a = sqrt(r1).*sin(2*pi*r2);
i = 1;
while Vpr<Vf*Lx*Lz %generation of agglomerated CNTs
    % Center points (cx,cy,cz) is random
    cx = Center1(1,Nrand(i)) + agg_rad*x_a(i);
    cy = 0;
    cz = Center1(3,Nrand(i)) + agg_rad*y_a(i);
    if cx>Lx || cx<0 % if the center is outside of RVE
        cx = (cx>Lx)*(Lx-0.01) + (cx<=Lx)*0.01; 
    end
    if cz>Lz || cz<0 % if the center is outside of RVE
        cz = (cz>Lz)*(Lz-0.01) + (cz<=Lz)*0.01; 
    end
    r = rand(2,1);
    phi = 2*pi*r(1,:);
    li = Wa*(-log(1-r(2,:))).^(1/Wb);
    Vp = D*li;
    x1 = cx-0.5*li.*(cos(phi));
    y1 = 0;
    z1 = cz-0.5*li.*(sin(phi));
    x2 = cx+0.5*li.*(cos(phi));
    y2 = 0;
    z2 = cz+0.5*li.*(sin(phi));
    % Apply Periodic BC
    pt= [x1 y1 z1;x2 y2 z2]; %two end points of CNT
    P = PBC_2D(pt,RVE);
    PPi1 = [PPi1 pt']; %% Matrix with all the CNT generated in the box without BCs
    Center1 = [Center1 [cx;cy;cz]]; %Matrix with all the center of CNT generated in the box without BCs
    Angle1 = [Angle1 phi];%Matrix with all the angle of CNT generated in the box without BC
    Length1 = [Length1 li];%Matrix with all the length of CNT generated in the box without BC
    PPi=[PPi P];   %% Matrix with all the CNT generated in the box with BCs
    Vpr=Vpr+Vp;
    i = i+1;
end
end
