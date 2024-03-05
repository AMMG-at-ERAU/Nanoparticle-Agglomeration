function [ dCNT, ja ] =A_CNT_contact_multiple_Grids_small( PP, Plength, RVE)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
% ///////////////////////// This is working but it is NOT FAST///////////////////////////
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
tic
N = length(Plength);
dCNT = sparse(N,N);
ja = dCNT;

trio = [1 2 3];
Nbin = RVE.nBin
ctobin = Nbin./RVE.size; %[ctobinx, ctobiny, ctobinz]
n = [1,0,0;0,1,0;0,0,1];
Grid = cell(Nbin,Nbin,Nbin); %Create cell of array with Nbin in each dimension

%% compute distances betwen CNTs
for i=1:N
    list = cell(1,1);
    Ptrue2 = (ones(2,1)*ctobin).*PP(2*i-1:2*i,:);%to scale everything for each bin with length 1 unit, to avoid rounding errors
    %Spatial patitioning
    Min = min(Ptrue2);%[xmin ymin zmin] 
    Max = max(Ptrue2);%[xmax ymax zmax]
        %find the bins this filler may overlap
        Minn = floor(Min)+1;
        Minn = Minn + (Minn < 1);
        Maxn = floor(Max)+1;
        Maxn = Maxn - (Maxn > Nbin);
        
        B = floor(Ptrue2)+1;
        B = B + (B < 1);
        B = B - (B > Nbin);
    %%%%%% Find bins filler intersects     
        for gi = trio
            n1 = n(gi,:);
            V0 = (Minn(gi):Maxn(gi)-1)' * n1;
            % V0 = [1 0 0;2 0 0;3 0 0;4 0 0]
            I = A_bin_line_plane(n1,V0,Ptrue2(1,:),Ptrue2(2,:));
            tr = trio([1:gi-1 gi+1:end]);
            b = floor(I(:,tr))+1;
            A = zeros(2*size(b,1),3);
            A(1:2:end,[gi tr])= [[Minn(gi):Maxn(gi)-1]' b];
            A(2:2:end,[gi tr])= [[Minn(gi)+1:Maxn(gi)]' b];
            B =[B;A];
        end
        C = unique(B,'rows')
     %%%%%% Add bins that are next to the ones fillers intersect (Bins
     %%%%%% reacheable with tunelling)
        for j = 1:size(C,1)
            C1(27*j-26:27*j,1) = [ones(9,1)*(C(j,1)-1);ones(9,1)*C(j,1);ones(9,1)*(C(j,1)+1)];
            C11 = [ones(3,1)*(C(j,2)-1);ones(3,1)*C(j,2);ones(3,1)*(C(j,2)+1)];
            C12 = [[C(j,3)-1:C(j,3)+1]';[C(j,3)-1:C(j,3)+1]';[C(j,3)-1:C(j,3)+1]'];
            C1(27*j-26:27*j,2:3) = [C11 C12;C11 C12;C11 C12];
        end
        temp = prod((C1>0 & C1<Nbin+1),2);
        C2 = unique(C1(temp>0,:),'rows');
        
        for ij = 1:size(C2,1)
            list{1} = [list{1} Grid{C2(ij,1),C2(ij,2),C2(ij,3)}]; %list of all CNTs in overlapped bins    
        end
        list2 = unique(cell2mat(list{1})); 
%%%%%%% distance
    if ~isempty(list2)
        [dt, dti, dta] = A_DistSeg2Seg(PP(2*i-1,:),PP(2*i,:),PP(2*list2-1,:),PP(2*list2,:),RVE.dvdw,RVE.D); 
        dCNT(list2,i) = (dt<=RVE.dcut+1e-6) .* dt;
        ja(list2,i) = (dt<=RVE.dcut+1e-6) .* dta;   %distance between junction on CNT a(j) and starting point of CNT a(j)
        ja(i,list2) = (dt<=RVE.dcut+1e-6) .* dti;   %distance between junction on CNT i and starting point of CNT i
    end
%%%%%%% put filler in only bin it intersects 
        for ij = 1:size(C,1)
        Grid{C(ij,1),C(ij,2),C(ij,3)} = [Grid{C(ij,1),C(ij,2),C(ij,3)} {i}];  
        end

end
toc
end

