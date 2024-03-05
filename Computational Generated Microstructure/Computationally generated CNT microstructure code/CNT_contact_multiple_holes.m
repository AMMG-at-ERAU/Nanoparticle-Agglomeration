function [ dCNT, ja ] =CNT_contact_multiple_holes( PP, Pcenter, Plength, R, Radius_holes, RVE)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
% ///////////////////////// This is working 2nd FAST///////////////////////////
% PP has all the CNT ends points as each row with 3 columns
% Pcenter is a 3 columns matrix with all centers
% Plength is a row vector with the CNT lengths
% dCNT is a N by N matrix with the tunneling disances
% ja is a N by N matrix with each row representing the different junction
% on the CNT of that row (their distance from the starting point of the
% CNT)
% N is the number of CNTs
% //////////////////////// Here A_DistSeg2Seg is used. It computes simultaneously the distance for all
% adjacents CNT to the current CNT /////////////////////////
N = length(Plength);
dCNT = sparse(N,N);
rr = RVE.rr(1:N,1:N); %tunelling distance
ja = dCNT;
rad_hole = Radius_holes;

%% compute distances betwen CNTs
for i=1:N
    Lspace = 2*RVE.D + RVE.dcut + 0.5*Plength(i) + (0.5*Plength(1:i-1))';     %neighbour range as a colomn vector
    temp = (ones(i-1,1)*Pcenter(i,:)>=(Pcenter(1:i-1,:)-Lspace*ones(1,3))) & (ones(i-1,1)*Pcenter(i,:)<=(Pcenter(1:i-1,:)+Lspace*ones(1,3))); 
    temp = prod(double(temp),2);                                % find adjcent CNTs
    temp = temp > 0;
    if ~isempty(find(temp, 1))
        a = find(temp);
        rtun = rr(a,i);
        [dt, dti, dta, PTi, PTa] = DistSeg2Seg_holes(PP(2*i-1,:),PP(2*i,:),PP(2*a-1,:),PP(2*a,:),rtun,RVE.D); 

        for ii = 1:length(a)
            dist1 = Dist_point2holes(R,PTi(ii,:)); %distance with all holes center
            dist2 = Dist_point2holes(R,PTa(ii,:));
            dist1 = (dist1<=rad_hole); % have values 1 everywhere it is less than radius
            dist2 = (dist2<=rad_hole);
            dist = sum(dist1)+sum(dist2);
            if (dist>=1)
                dt(ii) = 2*RVE.dcut + 2e-6; %Make it bigger then tunnelling distance so it acts as if not in percolated network
            end
        end
        dCNT(a,i) = (dt<=RVE.dcut+1e-6) .* dt;
        ja(a,i) = (dt<=RVE.dcut+1e-6) .* dta;   %distance between junction on CNT a(j) and starting point of CNT a(j)
        ja(i,a) = (dt<=RVE.dcut+1e-6) .* dti;   %distance between junction on CNT i and starting point of CNT i
    end
end

end

