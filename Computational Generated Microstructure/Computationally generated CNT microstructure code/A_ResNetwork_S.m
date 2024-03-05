function [ G, perco ] = A_ResNetwork_S(dCNT, ja, PP, RVE)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here
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
actCNT = zeros(1,length(dCNT)+1);                       % Active CNTs: those CNTs which are integrated in percolation cluster
ptr = A_cluster(dCNT);
root1 = zeros(1,length(ind1),'int8');
root2 = zeros(1,length(ind2),'int8')
for i=1:length(ind1)
    root1(i) = A_findroot(i,ptr); %find roots of CNTs on that edge
end
for i=1:length(ind2)
    root2(i) = A_findroot(i,ptr); %find roots of CNTs on that second edge
end
root = intersect(A,B);
p = 0;
for i = 1:length(root)
    R = A_findleafs(root(i),ptr); %find all leafs of that root
    actCNT(p+1:p+length(R)) = R;
    p = length(actCNT);
end
DCNT = dCNT + dCNT';                                    % distance matrix %Have symeetric matrix
dGraph = sparse(DCNT~=0);                               % define undirected graph
p = 0;
for i=1:length(ind1)                                    % pick out Active CNTs
    order = graphtraverse(dGraph, ind1(i));             % find CNTs that connected to the 1st edge(ISRC)
    if find(ismember(order, ind2),1)                    % check if CNTs in 'order' connect to the 2nd edge(GND)
        actCNT(p+1:p+length(order)) = order;            % feed actCNT
        actCNT = unique(actCNT);                        % find unique elements
        p = length(actCNT);
    end
end
perco = p;
if perco>0
    disp('// CNT network formed...');
else
    G = 0;
    disp('// Warning: No percolation cluster exist!');
    return
end

actCNT = actCNT(2:length(actCNT));                      % exclude 0
indtemp = ismember(1:length(dCNT), actCNT);             % those active CNTs
dGraph = dGraph(indtemp,indtemp);                       % exclude those unactive CNTs
indI = actCNT(ismember(actCNT,ind1));                   % conductive ISRC CNTs indices
indG = actCNT(ismember(actCNT,ind2));                   % conductive GND CNTs indices
indIG = [indI indG];                                    % ISRC & GND CNTs
innCNT = actCNT(~ismember(actCNT,indIG));               % inner CNTs (not the one touching electrodes)

tic
%% Delete CNTs that are dead ends
q = length(innCNT)+1;
while q>length(innCNT)                                  % repeat if exsit unonductive CNTs
    q = length(innCNT);                                 % # of conductive CNTs
    i = 1;
    while i<=length(innCNT)                             % pick out unonductive CNTs
        indtemp = find(ismember(actCNT,innCNT(i)));     % index of innCNT in actCNT
        if isempty(find(dGraph(indtemp,:), 1))          % if innCNT(i) has no connection %not possible
            continue;                                   % continue
        end
        order = graphtraverse(dGraph,indtemp,'depth',1);% CNTs that connected to innCNT(i), inCNT(i) is included
        p = 0;
        if length(order)<3                              % check if innCNT(i) is isolated, so it is connected only to only one other CNT
            dGraph(order(1),:) = 0;                     % if YES, clear innCNT(i)
            dGraph(:,order(1)) = 0;
            innCNT(i) = [];                             % clear innCNT(i)
        else
            Gtemp = dGraph;
            Gtemp(order(1),:) = 0;                      % cut the innCNT(i) connection
            Gtemp(:,order(1)) = 0;
            for j=2:length(order)                       % innCNT(i) connected to the percolation cluster
                tree = graphtraverse(Gtemp, order(j));  % check if order(j) connect to the percolation cluster (is it connected to anny other of the percolated cluster?)
                if ~isempty(ismember(tree,indIG))
                    p = p + 1;
                end
            end
            if p<2
                dGraph(order(1),:) = 0;                 % if p<2, clear inn(i)
                dGraph(:,order(1)) = 0;
                innCNT(i) = [];                         % clear innCNT(i)
            end
        end
        i = i + 1;
    end
end
innCNT;
toc
conCNT = [indI innCNT indG];                            % conductive CNTs
conN = length(conCNT);                                  % conductive CNT index

%% Calculation of Rt and Rint

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
        if LCNT<1e-6 %close to zero
            LCNT = 0.1;                                 % LCNT should be > 0
        end
        Gnode(IndN(i,q(indp(j-1))),IndN(i,q(indp(j)))) = (0.25*pi*RVE.sigma*RVE.D^2)*1e-6/LCNT;  % cnt    
    end
    if ismember(conCNT(i),indI)                         % ISRC node
        P1 = PP(2*conCNT(i)-1:2*conCNT(i),:);           % find the CNT
        P2 = P1(abs(P1(:,3))<=1e-9,:);                          % find end point on ISRC
        da = sqrt(sum((P2-P1(1,:)).*(P2-P1(1,:))));     % distance between point on ISRC and starting point of CNT
        [E,indp] = sort(abs(da-(conP(i,q)-0.1)));       % find distance (substract 0.1 added in distance between segment)
        LCNT = E(1) + (E(1)<1e-6)*0.1;                  % if length close to zero, add 0.1 to its value
        Gnode(IndN(i,q(indp(1))),node) = (0.25*pi*RVE.sigma*RVE.D^2)*1e-6/LCNT;
    end
    if ismember(conCNT(i),indG)                         % GND node
        P1 = PP(2*conCNT(i)-1:2*conCNT(i),:);           % find the CNT
        P2 = P1(abs(P1(:,3)-RVE.size(3))<=1e-9,:);                % find end point on GND
        da = sqrt(sum((P2-P1(1,:)).*(P2-P1(1,:))));     % distance between point on GND and starting point of CNT
        [E,indp] = sort(abs(da-(conP(i,q)-0.1)));         % find distance (substract 0.1 added in distance between segment)
        LCNT = E(1) + (E(1)<1e-6)*0.1;                  % if length close to zero, add 0.1 to its value
        Gnode(IndN(i,q(indp(1))),node+1) = (0.25*pi*RVE.sigma*RVE.D^2)*1e-6/LCNT;
    end
end
%% resistance matrix formation
Gnode = Gnode + Gnode';                                 % symmetrical
[a1,b1,c1]=find(abs(cJ)>9);
%cJ(abs(cJ)>9)
nInd=[];
for i= 1:length(c1)
    Gnode(IndN(a1(i),b1(i)),IndN(b1(i),a1(i)))=0; %delete the tunelling resistance for d<=dvdw to d=0
    Gnode(IndN(b1(i),a1(i)),IndN(a1(i),b1(i)))=0; %delete the tunelling resistance for d<=dvdw (since there is symmetry)
    Gnode(IndN(a1(i),b1(i)),:)=Gnode(IndN(a1(i),b1(i)),:)+Gnode(IndN(b1(i),a1(i)),:); %merge the two nodes
    Gnode(:,IndN(a1(i),b1(i)))=Gnode(:,IndN(a1(i),b1(i)))+Gnode(:,IndN(b1(i),a1(i)));
    nInd=[nInd IndN(b1(i),a1(i))];
end

    Gnode(nInd,:)=[]; %delete one of the nodes
    Gnode(:,nInd)=[];

Gnode = diag(sum(Gnode)) - Gnode;                       % conductance adjacent matrix
G = Gnode(1:length(Gnode)-1,1:length(Gnode)-1);                              % dimension adjust
end

