function [ dCNT, ja ] =A_CNT_contact_Grids_multiple( PP, Pcenter, Plength, RVE)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
% ///////////////////////// This is working but NOT FAST ///////////////////////////
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
% //////////////////////// Here Grids are used /////////////////////////
tic
N = length(Plength);
dCNT = sparse(N,N);
ja = dCNT;
Nbin = RVE.nBin;
ctobin = Nbin/RVE.size(1);
Grid = cell(Nbin,Nbin,Nbin); %Create cell of array with Nbin in each dimension

%% compute distances betwen CNTs
for i=1:N
    list = zeros(1,0);
    Min = min(PP(2*i-1:2*i,:));%[xmin ymin zmin]
    Max = max(PP(2*i-1:2*i,:));%[xmax ymax zmax]
        %find the bins this filler may overlap
        Minn = floor(Min*ctobin)+1;
        Minn = Minn + (Minn < 1);
        Maxn = floor(Max*ctobin)+1;
        Maxn = Maxn - (Maxn > Nbin);
        
        for ij = 1:length(Minn(1):Maxn(1))
            for ji = 1:length(Minn(3):Maxn(3))
            list = [list cell2mat(Grid(Minn(1)+ij-1,:,Minn(3)+ji-1))]; %list of all CNTs in overlapped bins
            end
        end
        list = unique(list);
        nlist = length(list);    
    Lspace = RVE.dcut + 0.5*Plength(i) + (0.5*Plength(list))';     %neighbour range as a colomn vector
    temp = (ones(nlist,1)*Pcenter(i,:)>(Pcenter(list,:)-Lspace*ones(1,3))) & (ones(nlist,1)*Pcenter(i,:)<(Pcenter(list,:)+Lspace*ones(1,3))); 
    temp = prod(double(temp),2);                                % find adjcent CNTs
    temp = temp > 0;
    if ~isempty(find(temp, 1))
        a = find(temp);
        [dt, dti, dta] = A_DistSeg2Seg(PP(2*i-1,:),PP(2*i,:),PP(2*list(a)-1,:),PP(2*list(a),:),RVE.dvdw,RVE.D); 
        dCNT(list(a),i) = (dt<=RVE.dcut+1e-6) .* dt;
        ja(list(a),i) = (dt<=RVE.dcut+1e-6) .* dta;   %distance between junction on CNT a(j) and starting point of CNT a(j)
        ja(i,list(a)) = (dt<=RVE.dcut+1e-6) .* dti;   %distance between junction on CNT i and starting point of CNT i
    end
        for ij = 1:length(Minn(1):Maxn(1))
            for iji = 1:length(Minn(2):Maxn(2))
                for ji = 1:length(Minn(3):Maxn(3))
                list1 =  [cell2mat(Grid(Minn(1)+ij-1,Minn(2)+iji-1,Minn(3)+ji-1)),i]; %list of all CNTs in that overlapped bin
                Grid(Minn(1)+ij-1,Minn(2)+iji-1,Minn(3)+ji-1) = {list1};
                end
            end
        end       
end
toc
end

