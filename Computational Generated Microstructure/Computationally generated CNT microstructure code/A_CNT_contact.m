function [ dCNT, ja ] =A_CNT_contact( PP, Pcenter, Plength, RVE)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
% ///////////////////////// This is working, FASTER ///////////////////////////
% PP has all the CNT ends points as each row with 3 columns
% Pcenter is a 3 columns matrix with all centers
% Plength is a row vector with the CNT lengths
% dCNT is a N by N matrix with the tunneling disances
% ja is a N by N matrix with each row representing the different junction
% on the CNT of that row (their distance from the starting point of the
% CNT)
% N is the number of CNTs
% //////////////////////// Here A_DistBetween2Segment_HC1 is used. It computes the distance for each
% adjacents CNT to the current CNT /////////////////////////
tic
N = length(Plength);
dCNT = sparse(N,N);
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
        [dt, dti, dta, ~] = A_DistBetween2Segment_HC1(PP(2*i-1,:),PP(2*i,:),PP(2*a(j)-1,:),PP(2*a(j),:),RVE.dvdw,RVE.D); 
        dCNT(a(j),i) = (dt<=RVE.dcut+1e-6) .* dt;
        ja(a(j),i) = (dt<=RVE.dcut+1e-6) .* dta;   %distance between junction on CNT a(j) and starting point of CNT a(j)
        ja(i,a(j)) = (dt<=RVE.dcut+1e-6) .* dti;   %distance between junction on CNT i and starting point of CNT i
        end
    end
end
toc
end

