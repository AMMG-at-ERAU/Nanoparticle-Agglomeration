function [ PP1, Pcenter1, Plength1, Pangle1 ] = A_HardCore( PP, Pcenter, Plength, Pangle, RVE)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
% PP has all the CNT ends points as each row with 3 columns
% Pcenter is a 3 columns matrix with all centers
% Pangle is a 2 columns matrix with all angles
% Plength is a row vector with the CNT lengths
% dCNT is a N by N matrix with the tunneling disances
% ja is a N by N matrix with each row representing the different junction
% on the CNT of that row (their distance from the starting point of the
% CNT)
% N is the number of CNTs
N = length(Plength);
dCNT = sparse(N,N);
ja = dCNT;
i = 1;                                                          % useful CNT's index
iL = 1;                                                         % CNT Library index
Overlap = 0;
Vpr=0;
%% compute distances betwen CNTs
while Vpr<RVE.Vf*prod(RVE.size)
    Pcenter1(i,:) = Pcenter(iL,:);
    Plength1(1,i) = Plength(1,iL);
    PP1(2*i-1:2*i,:) = PP(2*iL-1:2*iL,:);
    Pangle1(i,:) = Pangle(iL,:);
    k = 1;
    Lspace = RVE.dcut + 0.5*Plength1(i) + (0.5*Plength1(1:i-1))';     %neighbour range as a colomn vector
    temp = (ones(i-1,1)*Pcenter1(i,:)>(Pcenter1(1:i-1,:)-Lspace*ones(1,3))) & (ones(i-1,1)*Pcenter1(i,:)<(Pcenter1(1:i-1,:)+Lspace*ones(1,3))); 
    temp = prod(double(temp),2);                                % find adjcent CNTs
    temp = temp > 0;
    if ~isempty(find(temp, 1))      
        for j=1:length(find(temp))
        a = find(temp);
        [dt, ~, ~] = A_DistBetween2Segment(PP1(2*i-1,:),PP1(2*i,:),PP1(2*a(j)-1,:),PP1(2*a(j),:),RVE.dvdw,RVE.D);
        if dt<RVE.D+RVE.dvdw
           i = i-1;
           Overlap = Overlap + 1; 
           k = 0;
           break;
        end
       end
    end
    Vp=0.25*k*(pi*RVE.D^2)*Plength1(i);
    Vpr=Vpr+Vp;
    i = i + 1;
    iL = iL + 1;
end
Vpr
Overlap
end

