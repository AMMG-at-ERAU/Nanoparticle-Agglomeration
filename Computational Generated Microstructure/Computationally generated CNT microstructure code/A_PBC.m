function [ P, center, L ] = A_PBC(CNT,RVE)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%This  just cut off the CNT for every CNT, not the PBC and return
%the obtained CNT, their respective length and angle
%P(3,n);L(1,n/2);Angle(2,n/2),th is first row, phi is second row
%P is a 3 lines matrix
%center is a 3 lines vector
center=CNT.center';
L = CNT.length;
P1=CNT.pt(1,:);%starting point of CNT
P2=CNT.pt(2,:);%ending point of CNT
temp=(P1>=0 & P1<=RVE.size & P2>=0 & P2<=RVE.size);
temp1= P2>=0 & P2<=RVE.size;
temp2= P1>=0 & P1<=RVE.size;
P=[P1;P2];
if prod(temp)>0 %no end points out of RVE
    P=[P1' P2'];
    return
else
    if prod(temp1)==0 %end points is out 
        P=[P1;P2];
            [Pi, ~]=A_Line_Plane(RVE,P2,P1); %intersection and translation vector
            P=[P(1:end-1,:);Pi]; %forget the section outside

    end
    if prod(temp2)==0%  start point is out
        P1=P(1,:);
        P2=P(2,:);
            [Pi, ~]=A_Line_Plane(RVE,P1,P2); %intersection and translation vector
            P=[Pi;P(2:end,:)];
    end
    P= P';
end
n = size(P,2);
L = P(:,1:2:n-1)-P(:,2:2:n);
L = sqrt(sum(L.*L,1));
center = 0.5*(P(:,1:2:n-1)+P(:,2:2:n));
end