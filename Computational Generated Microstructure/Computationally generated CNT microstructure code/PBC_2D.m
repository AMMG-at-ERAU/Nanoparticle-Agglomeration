function [ P, center, L ] = PBC_2D(pt,RVE)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%This computes the Periodic Boundary Condition for every CNT, and return
%the obtained CNT, their respective length and angle
%P(3,n);L(1,n/2);Angle(2,n/2),th is first row, phi is second row
%%center = CNT.center';
%Angle=CNT.angle';
%%L = CNT.length;
center = pt.center';
L = pt.length;
P1 = pt.P(1,:);%starting point of CNT
P2 = pt.P(2,:);%ending point of CNT
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
        while prod(temp1)==0 %while end point is still outside RVE
            [Pi, Ti]=Line_Plane_2D(RVE,P2,P1); %intersection and translation vector
            P2=P2+Ti;
            P1=Pi+Ti;
            P=[P(1:end-1,:);Pi;P1;P2];
            temp1= P2>=0 & P2<=RVE.size;
        end
    end
    if prod(temp2)==0%  start point is out
        P1=P(1,:);
        P2=P(2,:);
        while prod(temp2)==0 %while start point is still outside RVE
            [Pi, Ti]=Line_Plane_2D(RVE,P1,P2); %intersection and translation vector
            P1=P1+Ti;
            P2=Pi+Ti;
            P=[P1;P2;Pi;P(2:end,:)];
            temp2= P1>=0 & P1<=RVE.size;
        end
    end
    P= P';
end
n = size(P,2);
L = P(:,1:2:n-1)-P(:,2:2:n);
L = sqrt(sum(L.*L,1));
center = 0.5*(P(:,1:2:n-1)+P(:,2:2:n));
% %Angle=Angle*ones(1,n/2);
end

