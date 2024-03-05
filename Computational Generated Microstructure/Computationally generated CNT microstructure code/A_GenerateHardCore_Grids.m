function [ Ptrue, Center, L, Ptrue1, Center1, L1 ] = A_GenerateHardCore_Grids( RVE)
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
i = 1;                                                          % useful CNT's index
ii = 1;
total=[];
Overlap = 0;
contact = 0;
Vpr=0;
Ptrue=[];
Center=[];
L = [];
Ptrue1 =[];
Center1 = [];
L1 = [];
Nbin = RVE.nBin
trio = [1 2 3];
ctobin = Nbin/RVE.size(1) %if Lx = Ly = Lz
Grid = cell(Nbin,Nbin,Nbin); %Create cell of array with Nbin in each dimension
RVE.Vf*prod(RVE.size)
%% compute distances betwen CNTs
while Vpr<RVE.Vf*prod(RVE.size)
% Generate CNT
r=rand(1,5);
randP=r(1,1:3).*[1,1,1];

center = RVE.size.*randP;
li = 5.0; %micrometer

u1 = 1.0-2.0*r(1,4);
v1 = sqrt(1.0-u1^2);
w1 = 2*pi*r(1,5);

P1=center-0.5*li*[v1*cos(w1),v1*sin(w1),u1];
P2=center+0.5*li*[v1*cos(w1),v1*sin(w1),u1];
    %Get the CNTs data before PBC
Ptrue(2*ii-1:2*ii,:) = [P1;P2];
Center(ii,:) = center;
L(1,ii) = li;

iii = 1;
k=0;
while k==0 && iii<101
list = zeros(1,0);
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
    %Spatial patitioning
    Min = min(Ptrue1(2*i-1:2*i,:));%[xmin ymin zmin]
    Max = max(Ptrue1(2*i-1:2*i,:));%[xmax ymax zmax]
        %find the bins this filler may overlap
        Minn = floor(Min*ctobin)+1;
        Minn = Minn + (Minn < 1);
        Maxn = floor(Max*ctobin)+1;
        Maxn = Maxn - (Maxn > Nbin);
        
        for ij = 1:length(Minn(1):Maxn(1))
            for ji = 1:length(Minn(3):Maxn(3))
            list = [list cell2mat(Grid(Minn(1)+ij-1,Minn(2):Maxn(2),Minn(3)+ji-1))]; %list of all CNTs in overlapped bins
            end
        end
        list = unique(list);
        nlist = length(list);    
    %Bounding volume
    Lspace = RVE.D + 0.5*L1(i) + (0.5*L1(list))';     %neighbour range as a colomn vector
    temp = (ones(nlist,1)*Center1(i,:)>(Center1(list,:)-Lspace*ones(1,3))) & (ones(nlist,1)*Center1(i,:)<(Center1(list,:)+Lspace*ones(1,3))); 
    temp = prod(double(temp),2);                                % find adjcent CNTs
    temp = temp > 0;
    if ~isempty(find(temp, 1))      
        for j=1:length(find(temp))
        a = find(temp);
%         if iii == 1
%             a
%             Ptrue1(2*list(a(j))-1,:)
%         end
        [dt, R] = A_DistBetween2Segment_HC(Ptrue1(2*i-1,:),Ptrue1(2*i,:),Ptrue1(2*list(a(j))-1,:),Ptrue1(2*list(a(j)),:),RVE.dvdw,RVE.D);
        if dt<RVE.dvdw
%             ade=a(j);
%             Ptrdeb=Ptrue1(2*a(j)-1,:);
            if iii<100
%            Cendeb=Center(ii,:);
%            R
           Center(ii,:)=Center(ii,:)+R;
%            cend=Center(ii,:);
           Ptrue(2*ii-1:2*ii,:) = Ptrue(2*ii-1:2*ii,:)+ones(2,1)*R;
           Overlap = Overlap + 1; 
           total(Overlap) = i;
           if prod(Center(ii,:)>=0 & Center(ii,:)<=RVE.size)==0 %center point out
               iii =101;
               total(total==i)=[]; %don't consider overlap due to that because the CNT is gonna be deleted
               i = i-1;
           ii = ii - 1;
           k = 0;
               break;
           end
           i = i-1;
           ii = ii - 1;
           k = 0;
%             else
%            Center1(i,:)=[];
%            L1(i) = [];
%            Ptrue1(2*i-1:2*i,:) = [];
%            
%            Center(ii,:)=[];
%            L(ii) = [];
%            Ptrue(2*ii-1:2*ii,:) = [];
%            Overlap = Overlap + 1  
            else
                contact = contact+1;
            end
%            i = i-1;
%            ii = ii - 1;
%            k = 0;
           iii = iii + 1;
           break;
        end
       end
    end

    Vp=0.25*k*(pi*RVE.D^2)*L1(i);
    if k == 1 %either distance between CNTs greatr than van der waal or iii>=100
        for ij = 1:length(Minn(1):Maxn(1))
            for iji = 1:length(Minn(2):Maxn(2))
                for ji = 1:length(Minn(3):Maxn(3))
%                     cm1=Minn(1)+ij-1
%                     cm2=Minn(2)+iji-1
%                     cm3=Minn(3)+ji-1
                list1 =  [cell2mat(Grid(Minn(1)+ij-1,Minn(2)+iji-1,Minn(3)+ji-1)),i]; %list of all CNTs in that overlapped bin
                Grid(Minn(1)+ij-1,Minn(2)+iji-1,Minn(3)+ji-1) = {list1};
                end
            end
        end        
    end
    Vpr=Vpr+Vp;
    i = i + 1;
    ii = ii + 1;
end

end
Vpr
Overlap
contact
tt = length(unique(total))
ratio = contact/tt
disp(['// Time for GenerateHardCore = ',num2str(toc),' second ']);
end

