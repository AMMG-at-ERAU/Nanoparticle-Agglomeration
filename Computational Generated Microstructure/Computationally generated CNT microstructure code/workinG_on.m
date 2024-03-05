function [ Ptrue, Center, L, Ptrue1, Center1, L1 ] = workinG_on( RVE)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
% ///////////////////////// This is working ///////////////////////////
% PP has all the CNT ends points as each row with 3 columns
% Pcenter is a 3 columns matrix with all centers
% Pangle is a 2 columns matrix with all angles
% Plength is a row vector with the CNT lengths
% dCNT is a N by N matrix with the tunneling disances
% ja is a N by N matrix with each row representing the different junction
% on the CNT of that row (their distance from the starting point of the
% CNT)
% N is the number of CNTs
% Nbin is number of bin in each axis
% ctobin is the number of bin per unit of length, assuming cubic bins(Nbin/Lx)
% ///////////////////// Grids are used //////////////////////////////
tic
format long
choice = 1e6; %Maximum number of particles expected
i = 1;                                                          % useful CNT's index
ii = 1;
Overlap = 0;
Vpr=0;
%Ptrue=[];
Ptrue = zeros(2*choice,3);
%Center=[];
Center = zeros(choice,3);
%L = [];
L = zeros(1,choice);
%Ptrue1 =[];
Ptrue1 = zeros(2*choice,3);
%Center1 = [];
Center1 = zeros(choice,3);
%L1 = [];
L1 = zeros(1,choice);
trio = [1 2 3];
Nbin = RVE.nBin
ctobin = Nbin/RVE.size(1)
n = [1,0,0;0,1,0;0,0,1];
Grid = cell(Nbin,Nbin,Nbin); %Create cell of array with Nbin in each dimension
RVE.Vf*prod(RVE.size)
%% compute distances betwen CNTs
while Vpr<RVE.Vf*prod(RVE.size)
% Generate CNT
r=rand(1,5);
randP=r(1,1:3).*[1,1,1];

center = RVE.size.*randP;
li = 0.15; %micrometer

u1 = 1.0-2.0*r(1,4);
v1 = sqrt(1.0-u1^2);
w1 = 2*pi*r(1,5);

P1=center-0.5*li*[v1*cos(w1),v1*sin(w1),u1];
P2=center+0.5*li*[v1*cos(w1),v1*sin(w1),u1];
    %Get the CNTs data before PBC
Ptrue(2*ii-1:2*ii,:) = [P1;P2];
Center(ii,:) = center;
L(1,ii) = li;

list = cell(1,1);
CNT.pt = Ptrue(2*ii-1:2*ii,:);
CNT.center = Center(ii,:);
CNT.length = L(1,ii);
% PBC    
[ P, center1, li1] = A_PBC(CNT,RVE); 
P = P';center1 = center1';

    k = 1;
    %Get the CNTs data after PBC
    Center1(i,:) = center1;
    L1(1,i) = li1;
    Ptrue1(2*i-1:2*i,:) = P;
    Ptrue2 = ctobin*P;%to scale everything for each bin with length 1 unit, to avoid rounding errors
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
        C = unique(B,'rows');
        
        for ij = 1:size(C,1)
            list{1} = [list{1} Grid{C(ij,1),C(ij,2),C(ij,3)}]; %list of all CNTs in overlapped bins    
        end
        list2 = unique(cell2mat(list{1}));
        nlist = length(list2);    

    if ~isempty(list2)      
        for j=1:nlist
            [dt, ~] = A_DistBetween2Segment_HC(Ptrue1(2*i-1,:),Ptrue1(2*i,:),Ptrue1(2*list2(j)-1,:),Ptrue1(2*list2(j),:),RVE.dvdw,RVE.D);
            if dt<=0
                ii = ii - 1;
                i = i - 1;
                Overlap = Overlap + 1;
                k = 0;
                break;
            end
            
        end
    end

    Vp=0.25*k*(pi*RVE.D^2)*L1(i);
    if k == 1 %either distance between CNTs greatr than van der waal or iii>=100
        for ij = 1:size(C,1)
        Grid{C(ij,1),C(ij,2),C(ij,3)} = [Grid{C(ij,1),C(ij,2),C(ij,3)} {i}];  
        end
%         if  mod(i,1000) ==0
%             i
%         end
               
    end
    Vpr=Vpr+Vp;
    i = i + 1;
    ii = ii + 1;

end
Ptrue = Ptrue(1:2*(ii-1),:);
Center = Center(1:ii-1,:);
L = L(1,1:ii-1);
Ptrue1 = Ptrue1(1:2*(i-1),:);
Center1 = Center1(1:i-1,:);
L1 = L1(1,1:i-1);
0.25*(pi*RVE.D^2)*sum(L1)
Vpr
Overlap
disp(['// Time for GenerateHardCore = ',num2str(toc),' second ']);
end

