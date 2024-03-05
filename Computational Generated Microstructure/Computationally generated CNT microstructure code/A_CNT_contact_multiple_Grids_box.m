function [ dCNT, ja ] =A_CNT_contact_multiple_Grids_box( PP, Plength, RVE)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
% ///////////////////////// This is working THIS IS FASTEST///////////////////////////
% PP has all the CNT ends points as each row with 3 columns
% Pcenter is a 3 columns matrix with all centers
% Plength is a row vector with the CNT lengths
% Angle is a 3 columns matrix (first is 
% dCNT is a N by N matrix with the tunneling disances
% ja is a N by N matrix with each row representing the different junction
% on the CNT of that row (their distance from the starting point of the
% CNT)
% N is the number of CNTs
% //////////////////////// Here A_DistSeg2Seg is used. It computes simultaneously the distance for all
% adjacents CNT to the current CNT /////////////////////////
N = length(Plength);
dCNT = sparse(N,N);
ja = dCNT;

L1 = PP(2:2:2*N,:)-PP(1:2:2*N-1,:); %CNT vector
unit = L1./(Plength'*ones(1,3)); % unit vector of each CNT
Angle(:,1) = unit(:,3); %u1
Angle(:,2) = sqrt(1-(unit(:,3)).^2);%v1
Angle(:,3) = atan2(unit(:,2),unit(:,1)) + 2*pi*(unit(:,2)<0);%w1
%Angle(isnan(Angle(:,3)),3) = 0;
trio = [1 2 3];
Nbin = RVE.nBin
ctobin = Nbin./RVE.size; %[ctobinx, ctobiny, ctobinz]
n = [1,0,0;0,1,0;0,0,1];
Grid = cell(Nbin,Nbin,Nbin); %Create cell of array with Nbin in each dimension

%% compute distances betwen CNTs
for i=1:N
    list = cell(1,1);
    Matrix = [Angle(i,1)*cos(Angle(i,3)) -sin(Angle(i,3)) Angle(i,2)*cos(Angle(i,3));Angle(i,1)*sin(Angle(i,3)) cos(Angle(i,3)) Angle(i,2)*sin(Angle(i,3));-Angle(i,2) 0 Angle(i,1)];
    base = 0.5*[RVE.dcut+RVE.D RVE.dcut+RVE.D -RVE.D-RVE.dcut -RVE.D-RVE.dcut;RVE.D+RVE.dcut -RVE.D-RVE.dcut RVE.D+RVE.dcut -RVE.D-RVE.dcut;-RVE.dcut -RVE.dcut -RVE.dcut -RVE.dcut]; %[X;Y;Z]four points around center but with tunneling
    X2=(Plength(i)+ RVE.dcut)*[Angle(i,2)*cos(Angle(i,3)),Angle(i,2)*sin(Angle(i,3)),Angle(i,1)]'; %add tunneling to length of CNT for CNT vector
    Z = zeros(3,8);
    Z(:,1:2:7) = Matrix*base + (PP(2*i-1,:))'*ones(1,4);%four starting points of the 4 CNTs all around (Matrix*base  find components in reference frame)
    Z(:,2:2:8) = Z(:,1:2:7) + X2*ones(1,4);%four end points of the 4 CNTs all around
    Z = Z'; %[X Y Z]
    BB = [];
    for j =1:4
    Ptrue2 = (ones(2,1)*ctobin).*Z(2*j-1:2*j,:);%to scale everything for each bin with length 1 unit, to avoid rounding errors
    %PP(2*i-1:2*i,:)
    %PP(2*i-1,:)-Plength(i)*[Angle(i,2)*cos(Angle(i,2)),Angle(i,2)*sin(Angle(i,3)),Angle(i,1)]
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
        BB = [BB;B];
    end
        C = unique(BB,'rows');
        C = C + (C < 1);
        C = C - (C > Nbin);
        
         for ij = 1:size(C,1)
            list{1} = [list{1} Grid{C(ij,1),C(ij,2),C(ij,3)}]; %list of all CNTs in overlapped bins         
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
end

