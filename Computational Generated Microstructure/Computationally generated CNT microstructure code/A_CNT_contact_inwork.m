function [ dCNT, ja ] =A_CNT_contact_multiple( PP, Pcenter, Plength, RVE)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
% PP has all the CNT ends points as each row with 3 columns
% Pcenter is a 3 columns matrix with all centers
% Plength is a row vector with the CNT lengths
% dCNT is a N by N matrix with the tunneling disances
% ja is a N by N matrix with each row representing the different junction
% on the CNT of that row (their distance from the starting point of the
% CNT)
% N is the number of CNTs
N = length(Plength);
dCNT = sparse(N,N);
ja = dCNT;
ptr = zeros(1,N,'int8');

%% compute distances betwen CNTs
for i=1:N
    Lspace = RVE.dcut + 0.5*Plength(i) + (0.5*Plength(1:i-1))';     %neighbour range as a colomn vector
    temp = (ones(i-1,1)*Pcenter(i,:)>(Pcenter(1:i-1,:)-Lspace*ones(1,3))) & (ones(i-1,1)*Pcenter(i,:)<(Pcenter(1:i-1,:)+Lspace*ones(1,3))); 
    temp = prod(double(temp),2);                                % find adjcent CNTs
    temp = temp > 0;
    if ~isempty(find(temp, 1))
        a = find(temp);
        [dt, dti, dta] = A_DistSeg2Seg(PP(2*i-1,:),PP(2*i,:),PP(2*a-1,:),PP(2*a,:),RVE.dvdw,RVE.D); 
        b = dt<=RVE.dcut+1e-6;
        dCNT(a,i) = b .* dt;
        ja(a,i) = b .* dta;   %distance between junction on CNT a(j) and starting point of CNT a(j)
        ja(i,a) = b .* dti;   %distance between junction on CNT i and starting point of CNT i
        
        if ~isempty(find(b, 1))
            c = a(b);
            for j=1:length(c)
                r2 = A_findroot(c(j), ptr); %find root of intersect CNT (follow A fast Monte Carlo algorithm for site or bond percolation)
                if r2 ~= i % CNT not in the three of CNT "i"
                    if ptr(i) > ptr(r2)
                        ptr(r2) = ptr(r2) + ptr(i); %merge number of CNT under the two threes
                        ptr(i) = r2; %root of i changes to index of the intersect
                    else
                        ptr(i) = ptr(i) + ptr(r2);
                        ptr(r2) = i;
                    end  
                end
            end  
        end
    end
end
end

