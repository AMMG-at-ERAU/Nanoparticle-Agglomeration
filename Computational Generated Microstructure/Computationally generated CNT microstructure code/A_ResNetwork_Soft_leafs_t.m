function [ G, perco ] = A_ResNetwork_Soft_leafs_t(dCNT, ja, PP, RVE)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here
% ///////////////////////// This is working ///////////////////////////
% G = A_ResNetwork_t(dCNT, PP, RVE)
%
%   RVE.sigma      %intrinsic conductivity S/m
%   dCNT - (nxn) tunelling distance matrix, [upper tri matrix]
%   ja is a N by N matrix with each row representing the different junction
%   on the CNT of that row (their distance from the starting point of the
%   CNT)
%   PP - 3 columns matrix with CNT end points coordinates
% 
% Returns:
%   G - Conductance Network Matrix
% ///////////////// This uses my own function, Faster /////////////////////
choice = 1e6;
if RVE.Dir=='z'
    za = PP(1:2:end,3);
    zb = PP(2:2:end,3);
    ind1 = find((abs(za)<=1e-9)|(abs(zb)<=1e-9));                              % z cross the 0 Boundary (bottom)
    ind2 = find(abs(za-RVE.size(3))<=1e-9|abs(zb-RVE.size(3))<=1e-9);          % z cross the Lz edge, connected to GND
elseif RVE.Dir=='y'
    ya = PP(1:2:end,2);
    yb = PP(2:2:end,2);
    ind1 = find((ya==0)|(yb==0));         % y cross the 1st edge, connected to ISRC
    ind2 = find((ya==RVE.size(2))|(yb==RVE.size(2)));            % y cross the 2nd edge, connected to GND
else
    xa = PP(1:2:end,1);
    xb = PP(2:2:end,1);
    ind1 = find((xa==0)|(xb==0));         % x cross the 1st edge, connected to ISRC
    ind2 = find((xa==RVE.size(1))|(xb==RVE.size(1)));            % x cross the 2nd edge, connected to GND
end
%indE = [ind1; ind2];                                    % Edge CNTs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Deal with tunnelling segments
%% pick out Active CNTs
%actCNT = zeros(1,length(dCNT));                       % Active CNTs: those CNTs which are integrated in percolation cluster
actCNT = zeros(1,choice);
DCNT = dCNT + dCNT';                                    % distance matrix %Have symeetric matrix
DCNT1 = DCNT;
ptr = A_clusters(dCNT);
%%
root1 = zeros(1,length(ind1));
root2 = zeros(1,length(ind2));
for i=1:length(ind1)
    [root1(i),ptr] = A_find_root(ind1(i),ptr); %find roots of CNTs on that edge
end
for i=1:length(ind2)
    [root2(i),ptr] = A_find_root(ind2(i),ptr); %find roots of CNTs on that second edge
end
root = intersect(root1,root2); %find roots of clusters that has CNT from both egdes
if ~isempty(root)
    disp('// CNT network formed...');
else
    G = 0;
    perco = 0;
    disp('// Warning: No percolation cluster exist!');
    return
end
%%
p = 0;
for i = 1:length(root)
    R = A_find_leafs(root(i),ptr); %find all leafs of that root
%     actCNT(p+1:p+length(R)) = R;
%     p = length(actCNT);
    actCNT(p+1:p+length(R)) = R;
    p = p + length(R);

end
perco = p;
%actCNT = actCNT(actCNT>0);
actCNT = actCNT(1:p);
actCNT = sort([root actCNT]);
%indtemp = ismember(1:length(dCNT), actCNT);             % those active CNTs
%DCNT1 = DCNT(indtemp,indtemp)
indI = actCNT(ismember(actCNT,ind1));                   % conductive ISRC CNTs indices
indG = actCNT(ismember(actCNT,ind2));                   % conductive GND CNTs indices
indIG = [indI indG];                                    % ISRC & GND CNTs


%% Delete CNTs that are dead ends
me = length(DCNT1);
DCNT1(indIG,me+1)=1;
DCNT1(me+1,indIG)=1;
[~ , C] = biconnected_components(DCNT1);
clus = full(unique(C(me+1,:)));
clus = clus(clus~=0);
%row =[];
row = zeros(choice,1);
p = 0;
for i = 1:length(clus)
[row1,col1] = find(C==clus(i));
row(p+1:p+2*length(col1)) = [row1;col1];
p = p + 2*length(col1);
%row = [row;row1;col1];
end
%row = unique(row);
row = unique(row(1:p));
innCNT = (row(~ismember(row,[indIG me+1])))';

conCNT = [indI innCNT indG];                            % conductive CNTs
conN = length(conCNT);                                  % conductive CNT index

%% conductance matrix formation

%Value of conductance at electrodes
LCNT = 0.5;
%Conduc = (0.25*pi*RVE.sigma*RVE.D^2)*1e-6/LCNT
                                                        % Gnode structure:
                                                        % |- Inner nodes -|- ISRC node -|- GND node -|
[p,q,d] = find(DCNT(conCNT, conCNT));                   % find effective junctions
g = 1./A_ResTun1( RVE, d );                             % find tunnelling conductance
GMat = sparse(p,q,g,conN+2,conN+2);                     % conductance matrix between CNTs
p = ismember(conCNT,indI);
q = ismember(conCNT,indG);
GMat(p,1+conN) = 9.7e-07;                                     % The CNTs connected to IsrcE are connected to ISRC node by RK
GMat(1+conN,p) = 9.7e-07;
GMat(q,2+conN) = 9.7e-07;                                     % The CNTs connected to GndE are connected to GndE node by RK in THIS program
GMat(2+conN,q) = 9.7e-07;

GMat = -GMat + diag(sum(GMat));                         % conductance adjacent matrix
G = GMat(1:1+conN, 1:1+conN);                           % dimension adjust
end

