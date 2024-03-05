clc
clear all
format long
RVE.size = [5 5 5];
RVE.Dir = 'z'; 
X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0;0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0;0,0;RVE.size(1),RVE.size(1);0,0;RVE.size(1),RVE.size(1)];
        Y=[0,RVE.size(2);RVE.size(2),RVE.size(2);RVE.size(2),0;0,0;0,RVE.size(2);RVE.size(2),RVE.size(2);RVE.size(2),0;0,0;0,0;0,0;RVE.size(2),RVE.size(2);RVE.size(2),RVE.size(2)];
        Z=[0,0;0,0;0,0;0,0;RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);0,RVE.size(3);0,RVE.size(3);0,RVE.size(3);0,RVE.size(3)];


PP = [1 1 3 4 0.9 2.1 3 3 4 4 2.9 4.1 2 2 4 5 2 3 2.5 2.5 2.1 3.4 3 3 2.9 3.9 2 2 1 1 3.75 3.75 3.5 4.5 1.9 3 4 4 0.9 2.1 4.7 4.7 4.1 4.9
    2.5 2.5 2.6 2.6 2.6 2.6 2.5 2.5 2.5 2.5 2.6 2.6 2.5 2.5 2.6 2.6 2.6 2.6 2.5 2.5 2.6 2.6 2.5 2.5 2.6 2.6 2.5 2.5 2.5 2.5 2.5 2.5 2.6 2.6 2.6 2.6 2.5 2.5 2.6 2.6 2.5 2.5 2.6 2.6
    5 3.7 4 5 4 4 4.1 2.9 4 3 3.5 3.5 4.1 2.9 3 3 3 3 3.1 1.9 2 2 2.2 0.9 1.5 1.5 1.1 0 1 0 2 0.8 1 1 1 1 1.1 0 0.5 0.5 3.2 2 2.3 2.3]
    Pcenter= 0.5*(PP(:,1:2:end)+PP(:,2:2:end))
    Plength=(sum((PP(:,1:2:end)-PP(:,2:2:end)).^2,1)).^0.5
    N = size(Plength,2)
    RVE.dcut = 0.10;
    RVE.dvdw = 0.01
    RVE.D = 0
   
     for i=1:N
            Xx(i,:)=PP(1,2*i-1:2*i);
            Yy(i,:)=PP(2,2*i-1:2*i);
            Zz(i,:)=PP(3,2*i-1:2*i);
        end
        figure
        plot3(X',Y',Z','r')
        hold on
        plot3(Xx',Yy',Zz','k')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
PP=PP';
dCNT = sparse(N,N);
Pcenter = Pcenter'
ja = dCNT;

%% compute distances betwen CNTs
for i=1:N
    Lspace = RVE.dcut + 0.5*Plength(i) + (0.5*Plength(1:i-1))';     %neighbour range as a colomn vector
    temp = (ones(i-1,1)*Pcenter(i,:)>(Pcenter(1:i-1,:)-Lspace*ones(1,3))) & (ones(i-1,1)*Pcenter(i,:)<(Pcenter(1:i-1,:)+Lspace*ones(1,3))); 
    temp = prod(double(temp),2);                                % find adjcent CNTs
    temp = temp > 0;
    if ~isempty(find(temp, 1))
        for j=1:length(find(temp))
        a = find(temp);
        [dt, dti, dta] = A_DistBetween2Segment_Untitled1(PP(2*i-1,:),PP(2*i,:),PP(2*a(j)-1,:),PP(2*a(j),:),RVE.dvdw,RVE.D);
        dCNT(a(j),i) = (dt<=RVE.dcut+1e-6) .* dt;
        ja(a(j),i) = (dt<=RVE.dcut+1e-6) .* dta;   %distance between junction on CNT a(j) and starting point of CNT a(j)
        ja(i,a(j)) = (dt<=RVE.dcut+1e-6) .* dti;   %distance between junction on CNT i and starting point of CNT i
        end
    end
end

%%


if RVE.Dir=='z'
    za = PP(1:2:end,3);
    zb = PP(2:2:end,3);
    ind1 = find((za==0)|(zb==0));                              % z cross the 0 Boundary (bottom)
    ind2 = find((za==RVE.size(3))|(zb==RVE.size(3)));          % z cross the Lz edge, connected to GND
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
DCNT = dCNT + dCNT';                                    % distance matrix %Have symeetric matrix
dGraph = sparse(DCNT~=0);                               % define undirected graph
me = length(DCNT)
p = 0;
DCNT([ind1;ind2],me+1)=1;
DCNT(me+1,[ind1;ind2])=1;
[a C] = biconnected_components(DCNT)
%biconncomp(DCNT)
% for i=1:length(ind1)                                    % pick out Active CNTs
%     order = graphtraverse(dGraph, ind1(i));             % find CNTs that connected to the 1st edge(ISRC)
%     if find(ismember(order, ind2),1)                    % check if CNTs in 'order' connect to the 2nd edge(GND)
%         actCNT(p+1:p+length(order)) = order;            % feed actCNT
%         actCNT = unique(actCNT);                        % find unique elements
%         p = length(actCNT);
%     end
% end
% perco = p;
% if perco>0
%     disp('// CNT network formed...');
% else
%     G = 0;
%     disp('// Warning: No percolation cluster exist!');
%     return
% end
% 
% actCNT = actCNT(2:length(actCNT));                      % exclude 0
% indtemp = ismember(1:length(dCNT), actCNT);             % those active CNTs
% dGraph = dGraph(indtemp,indtemp);                       % exclude those unactive CNTs
% indI = actCNT(ismember(actCNT,ind1));                   % conductive ISRC CNTs indices
% indG = actCNT(ismember(actCNT,ind2));                   % conductive GND CNTs indices
% indIG = [indI indG];                                    % ISRC & GND CNTs
% innCNT = actCNT(~ismember(actCNT,indIG));               % inner CNTs (not the one touching electrodes)
% 
% 
% %% Delete CNTs that are dead ends
% q = length(innCNT)+1;
% while q>length(innCNT)                                  % repeat if exsit unonductive CNTs
%     q = length(innCNT);                                 % # of conductive CNTs
%     i = 1;
%     while i<=length(innCNT)                             % pick out unonductive CNTs
%         indtemp = find(ismember(actCNT,innCNT(i)));     % index of innCNT in actCNT
%         if isempty(find(dGraph(indtemp,:), 1))          % if innCNT(i) has no connection %not possible
%             continue;                                   % continue
%         end
%         order = graphtraverse(dGraph,indtemp,'depth',1);% CNTs that connected to innCNT(i), inCNT(i) is included
%         p = 0;
%         if length(order)<3                              % check if innCNT(i) is isolated, so it is connected only to only one other CNT
%             dGraph(order(1),:) = 0;                     % if YES, clear innCNT(i)
%             dGraph(:,order(1)) = 0;
%             innCNT(i) = [];                             % clear innCNT(i)
%         else
%             Gtemp = dGraph;
%             Gtemp(order(1),:) = 0;                      % cut the innCNT(i) connection
%             Gtemp(:,order(1)) = 0;
%             for j=2:length(order)                       % innCNT(i) connected to the percolation cluster
%                 tree = graphtraverse(Gtemp, order(j));  % check if order(j) connect to the percolation cluster (is it connected to anny other of the percolated cluster?)
%                 if ~isempty(ismember(tree,indIG))
%                     p = p + 1;
%                 end
%             end
%             if p<2
%                 dGraph(order(1),:) = 0;                 % if p<2, clear inn(i)
%                 dGraph(:,order(1)) = 0;
%                 innCNT(i) = [];                         % clear innCNT(i)
%             end
%         end
%         i = i + 1;
%     end
% end
% conCNT = [indI innCNT indG];                            % conductive CNTs
% conN = length(conCNT);                                  % conductive CNT index
% 
% %% Calculation of Rt and Rint
% 
% conJun = DCNT(conCNT, conCNT);                          % find effective junctions
% conP = ja;                                              
% conP = conP(conCNT, conCNT);
% IndN = sparse(conN, conN);                              % node index in conJun
% node = 1;
% for i=1:conN                                            % node index defintaion
%     [~,q,d] = find(conJun(i,i+1:conN));
%     for j=1:length(d)
%         IndN(i,i+q(j)) = node;
%         IndN(i+q(j),i) = node + 1;
%         node = node + 2;
%     end
% end
% Gnode = sparse(node+1,node+1);                          % Gnode matrix:
%                                                         % |- Inner nodes-|- ISRC node -|- GND node -|
% for i=1:conN                                            % Inner nodes, Gjunction
%     [~,q,d] = find(conJun(i,i+1:conN));
%     for j=1:length(d)
%         Gnode(IndN(i,i+q(j)),IndN(i+q(j),i)) = d(j);  % find tunnelling conductance
%     end
% end
% 
% for i=1:conN                                            % Inner nodes, Gcnt
%     [~,q,d] = find(conJun(i,:));
%     tempP = conP(i,:);
%     [sortp,indp] = sort(tempP(tempP~=0));
%     for j=2:length(d)
%         LCNT = abs(sortp(j-1)-sortp(j));
%         if LCNT<1e-6 %close to zero
%             LCNT = 0.1;                                 % LCNT should be > 0
%         end
%         Gnode(IndN(i,q(indp(j-1))),IndN(i,q(indp(j)))) = LCNT;  % cnt    
%     end
%     if ismember(conCNT(i),indI)                         % ISRC node
%         conCNT(i)
%         P1 = PP(2*conCNT(i)-1:2*conCNT(i),:);           % find the CNT
%         P2 = P1(P1(:,3)==0,:)                          % find end point on ISRC
%         da = sqrt(sum((P2-P1(1,:)).*(P2-P1(1,:))))     % distance between point on ISRC and starting point of CNT
%         conP(i,q)
%         [E,indp] = sort(abs(da-(conP(i,q)-0.5)))          % find distance (substract 0.1 added in distance between segment)
%         LCNT = E(1) + (E(1)<1e-6)*0.1                  % if length close to zero, add 0.1 to its value
%         Gnode(IndN(i,q(indp(1))),node) = LCNT;
%     end
%     if ismember(conCNT(i),indG)                         % GND node
%         conCNT(i)
%         P1 = PP(2*conCNT(i)-1:2*conCNT(i),:);           % find the CNT
%         P2 = P1(P1(:,3)==RVE.size(3),:)                % find end point on GND
%         da = sqrt(sum((P2-P1(1,:)).*(P2-P1(1,:))))     % distance between point on GND and starting point of CNT
%         conP(i,q)
%         [E,indp] = sort(abs(da-(conP(i,q)-0.5)));          % find distance (substract 0.1 added in distance between segment)
%         LCNT = E(1) + (E(1)<1e-6)*0.1                  % if length close to zero, add 0.1 to its value
%         Gnode(IndN(i,q(indp(1))),node+1) = LCNT;
%     end
% end
% %% resistance matrix formation
% Gnode1 = full(Gnode);
% Gnode = Gnode + Gnode';                                 % symmetrical
% Gnode = diag(sum(Gnode)) - Gnode;                       % conductance adjacent matrix
% G = Gnode(1:node,1:node);      