function [ Ptrue, Center, L, Ptrue1, Center1, L1 ] = A_HardCore_BdB_PBC( RVE)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
tic
choice = 1e6; %Maximum number of particles expected
i = 1;                                                          % useful CNT's index
ii = 1;
%unit = [1 0 0;0 1 0;0 0 1];
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
BDB = zeros(choice,6);
% trio = [1 2 3];
% Nbin = RVE.nBin
% ctobin = Nbin./RVE.size %[ctobinx, ctobiny, ctobinz]
% n = [1,0,0;0,1,0;0,0,1];
% Grid = cell(Nbin,Nbin,Nbin); %Create cell of array with Nbin in each dimension
% RVE.Vf*prod(RVE.size)
% %Weibul parameters in micrometer
% aW = 5.6403;
% bW = 2.4;

while Vpr<RVE.Vf*prod(RVE.size)
% Generate CNT
r=rand(1,5);
%r=rand(1,6);
randP=r(1,1:3).*[1,1,1];

center = RVE.size.*randP;
li = 0.4; %micrometer
%li=aW*(-log(1-r(1,6)))^(1/bW);

u1 = 1.0-2.0*r(1,4);
v1 = sqrt(1.0-u1^2);
w1 = 2*pi*r(1,5);
% Matrix = [u1*cos(w1) -sin(w1) v1*cos(w1);u1*sin(w1) cos(w1) v1*sin(w1);-v1 0 u1];%rotation matrix (Local to global)
% Runit = Matrix*unit; %unit vector of CNT after rotation [U V W]
P1=center-0.5*li*[v1*cos(w1),v1*sin(w1),u1];
P2=center+0.5*li*[v1*cos(w1),v1*sin(w1),u1];
    %Get the CNTs data before PBC
Ptrue(2*ii-1:2*ii,:) = [P1;P2];
Center(ii,:) = center;
L(1,ii) = li;


%list = cell(1,1);
CNT.pt = Ptrue(2*ii-1:2*ii,:);
CNT.center = Center(ii,:);
CNT.length = L(1,ii);
% PBC    
[ P, center1, li1] = A_PBC1(CNT,RVE); 
P = P';center1 = center1';
sizC = size(center1,1);
for jj =1:sizC
    k = 1;
    %Get the CNTs data after PBC
    Center1(i,:) = center1(jj,:);
    L1(1,i) = li1(jj);
    Ptrue1(2*i-1:2*i,:) = P(2*jj-1:2*jj,:);
    vt = [min(P(2*jj-1:2*jj,:)) max(P(2*jj-1:2*jj,:))] + 0.5*RVE.D*[-1 -1 -1 1 1 1];
    BD(jj,:) = [0.5*(vt(1)+vt(4)) 0.5*(vt(2)+vt(5)) 0.5*(vt(3)+vt(6)) -vt(1)+vt(4) -vt(2)+vt(5) -vt(3)+vt(6)];%[vtcenter(x,yz) vtlength(x,y,z)]
    if i>1
        
        temp = (ones(i-1,1)*BD(jj,1:3)>=(BDB(1:i-1,1:3)-0.5*(ones(i-1,1)*BD(jj,4:6)+BDB(1:i-1,4:6)))) & (ones(i-1,1)*BD(jj,1:3)<=(BDB(1:i-1,1:3)+0.5*(ones(i-1,1)*BD(jj,4:6)+BDB(1:i-1,4:6))));
        temp = prod(double(temp),2);                                % find adjcent GNPs
        temp = temp > 0;
        if ~isempty(find(temp, 1))
            a = find(temp);
            dt = A_DistBetween2Segment_multiple(Ptrue1(2*i-1,:),Ptrue1(2*i,:),Ptrue1(2*a-1,:),Ptrue1(2*a,:),RVE.dvdw,RVE.D);
            if prod(dt>0)==0%prod(dt>=RVE.dvdw)==0
                i = i-1;
                ii = ii - 1;
                Overlap = Overlap + 1;
                k = 0;
                sizC = 1;
                break;
            end
        end
    end
    
end

for jj1 =1:sizC
    
    Vp = 0;
    if k == 1 %either distance between CNTs greatr than van der waal or iii>=100
        
    Center1(i,:) = center1(jj1,:);
    L1(1,i) = li1(jj1);
    Ptrue1(2*i-1:2*i,:) = P(2*jj1-1:2*jj1,:);  
    BDB(i,:) = BD(jj1,:);
    Vp=0.25*(pi*RVE.D^2)*L1(i);  
    end
    Vpr=Vpr+Vp;
    i = i + 1;
end
    ii = ii + 1;

end
%since we preallocate choice slots of data points in the following
%variable, the non used slots will be zero. We need to take those away
Ptrue = Ptrue(1:2*(ii-1),:);
Center = Center(1:ii-1,:);
L = L(1,1:ii-1);
Ptrue1 = Ptrue1(1:2*(i-1),:);
Center1 = Center1(1:i-1,:);
L1 = L1(1,1:i-1);

Vpr
Overlap
disp(['// Time for GenerateHardCore = ',num2str(toc),' second ']);
end

