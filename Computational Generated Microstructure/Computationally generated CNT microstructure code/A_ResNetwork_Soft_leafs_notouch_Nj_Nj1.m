function [Nj, Nj1, Nj2, Lj, Lj1, Lj2, G, perco ] = A_ResNetwork_Soft_leafs_notouch_Nj_Nj1(dCNT, ja, PP, RVE)
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
Nj = length(find(dCNT)) %non zero distances (number of tunneling contact)
Lj = length(dCNT);
DCNT = dCNT + dCNT';                                    % distance matrix %Have symeetric matrix
DCNT1 = DCNT;
Nj1 = 0;
Lj1 = 0;
Nj2 = 0;
Lj2 = 0;
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
Lj1 = conN;
Lj2 = length(innCNT);
%% Calculation of Rt and Rint

Nj1 = 0.5*length(find(DCNT(conCNT, conCNT)));
Nj2 = 0.5*length(find(DCNT(innCNT,innCNT)));
conJun = DCNT(conCNT, conCNT);                          % find effective junctions
cJ = dCNT(conCNT, conCNT);
conP = ja;                                              
conP = conP(conCNT, conCNT);
IndN = sparse(conN, conN);                              % node index in conJun
node = 1;
for i=1:conN                                            % node index defintaion
    [~,q,d] = find(conJun(i,i+1:conN));
    for j=1:length(d)
        IndN(i,i+q(j)) = node;
        IndN(i+q(j),i) = node + 1;
        node = node + 2;
    end
end
Gnode = sparse(node+1,node+1);                          % Gnode matrix:
                                                        % |- Inner nodes-|- ISRC node -|- GND node -|
for i=1:conN                                            % Inner nodes, Gjunction
    [~,q,d] = find(conJun(i,i+1:conN));
    for j=1:length(d)
        Gnode(IndN(i,i+q(j)),IndN(i+q(j),i)) = 1./A_ResTun1( RVE, d(j) );  % find tunnelling conductance
    end
end

for i=1:conN                                            % Inner nodes, Gcnt
    [~,q,d] = find(conJun(i,:));
    tempP = conP(i,:);
    [sortp,indp] = sort(tempP(tempP~=0));
    for j=2:length(d)
        LCNT = abs(sortp(j-1)-sortp(j));
        if LCNT<1e-5;%1e-6 %close to zero
            LCNT = 1e-5;%0.1;                                 % LCNT should be > 0
        end
        Gnode(IndN(i,q(indp(j-1))),IndN(i,q(indp(j)))) = (0.25*pi*RVE.sigma*RVE.D^2)*1e-6/LCNT;  % cnt    
    end
    if ismember(conCNT(i),indI)                         % ISRC node
        P1 = PP(2*conCNT(i)-1:2*conCNT(i),:);           % find the CNT
        P2 = P1(abs(P1(:,3))<=1e-9,:);                          % find end point on ISRC
        da = sqrt(sum((P2-P1(1,:)).*(P2-P1(1,:))));     % distance between point on ISRC and starting point of CNT
        [E,indp] = sort(abs(da-(conP(i,q)-1e-5)));       % find distance (substract 0.1 added in distance between segment)
        LCNT = E(1) + (E(1)<1e-6)*1e-5;                  % if length close to zero, add 0.1 to its value
        Gnode(IndN(i,q(indp(1))),node) = (0.25*pi*RVE.sigma*RVE.D^2)*1e-6/LCNT;
    end
    if ismember(conCNT(i),indG)                         % GND node
        P1 = PP(2*conCNT(i)-1:2*conCNT(i),:);           % find the CNT
        P2 = P1(abs(P1(:,3)-RVE.size(3))<=1e-9,:);                % find end point on GND
        da = sqrt(sum((P2-P1(1,:)).*(P2-P1(1,:))));     % distance between point on GND and starting point of CNT
        [E,indp] = sort(abs(da-(conP(i,q)-1e-5)));         % find distance (substract 0.1 added in distance between segment)
        LCNT = E(1) + (E(1)<1e-6)*1e-5;                  % if length close to zero, add 0.1 to its value
        Gnode(IndN(i,q(indp(1))),node+1) = (0.25*pi*RVE.sigma*RVE.D^2)*1e-6/LCNT;
    end
end
%% resistance matrix formation
Gnode = Gnode + Gnode';                                 % symmetrical

% %E
% [a1,b1,c1]=find(abs(cJ)>9);
% %cJ(abs(cJ)>9)
% nInd=[];
% for i= 1:length(c1)
%     Gnode(IndN(a1(i),b1(i)),IndN(b1(i),a1(i)))=0; %delete the tunelling resistance for d<=dvdw to d=0
%     Gnode(IndN(b1(i),a1(i)),IndN(a1(i),b1(i)))=0; %delete the tunelling resistance for d<=dvdw (since there is symmetry)
%     Gnode(IndN(a1(i),b1(i)),:)=Gnode(IndN(a1(i),b1(i)),:)+Gnode(IndN(b1(i),a1(i)),:); %merge the two nodes
%     Gnode(:,IndN(a1(i),b1(i)))=Gnode(:,IndN(a1(i),b1(i)))+Gnode(:,IndN(b1(i),a1(i)));
%     nInd=[nInd IndN(b1(i),a1(i))];
% end
% 
%     Gnode(nInd,:)=[]; %delete one of the nodes
%     Gnode(:,nInd)=[];

Gnode = diag(sum(Gnode)) - Gnode;                       % conductance adjacent matrix
G = Gnode(1:node,1:node);                              % dimension adjust
end

