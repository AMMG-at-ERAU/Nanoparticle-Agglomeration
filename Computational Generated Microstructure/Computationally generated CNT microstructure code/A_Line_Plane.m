function [P, T]=A_Line_Plane(RVE,P1,P2)
%New_plane_line_intersect computes the intersection of a plane with BCs and a segment
% Inputs: 
%       ni: normal vector of the Plane i
%       V0: A point that belongs to simultaneously the Planes 1,3,6
%       V1: A point that belongs to simultaneously the Planes 2,4,5
%       P1: Initial point the segment P0P1
%       P2: End point of the segment P0P1
%       Lx: Length of the cubicoid side
%
%Outputs:
%      P0    is the point of interection 
%      T     is the translation vector for the periodic Boundary Condition
%  
%
% Example:
% Determine the intersection of following the plane x+y+z+3=0 with the segment P0P1:
% The plane is represented by the normal vector n=[1 1 1]
% and an arbitrary point that lies on the plane, ex: V0=[1 1 -5]
% The segment is represented by the following two points
% P0=[-5 1 -1]
%P1=[1 2 3]   
% [I,check]=plane_line_intersect([1 1 1],[1 1 -5],[-5 1 -1],[1 2 3]);

%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%If you have any comments or face any problems, please feel free to leave
%your comments and i will try to reply to you as fast as possible.
V0=RVE.V(1,:);
V1=RVE.V(2,:);
n1=RVE.n(1,:);
n2=RVE.n(2,:);
n3=RVE.n(3,:);
n4=RVE.n(4,:);
n5=RVE.n(5,:);
n6=RVE.n(6,:);
Lx=RVE.size(1);
Ly=RVE.size(2);
Lz=RVE.size(3);
%%%% Plane 1
u = P2-P1;
w = P1 - V0;
D = n1*u';
N = -n1*w';
sI = N / D;

Pn=[];%Take the intersection points (maximum 2)
Tn=[];%take the respective translation vector
if (abs(D) < 10^-7 || sI < 0 || sI > 1)   % The segment is parallel to plane or The intersection point  lies outside the segment, so there is no intersection
 
else
    
    P0 = P1+ sI.*u;  %The intersection point
    T=[0 0 Lz];      %The translation vector for periodic BC
    if (P0(1)>=0 && P0(1)<=Lx && P0(2)>=0 && P0(2)<=Ly) % Make sure that the intersection point is between the BC of the plane
    Pn=[Pn;P0];
    Tn=[Tn;T];
    end
end

%%%% Plane 2
u = P2-P1;
w = P1 - V1;
D = n2*u';
N = -n2*w';
sI = N / D;
if (abs(D) < 10^-7 || sI < 0 || sI > 1)   % The segment is parallel to plane or The intersection point  lies outside the segment, so there is no intersection

else
    P0 = P1+ sI.*u;  %The intersection point
    T=[0 0 -Lz];      %The translation vector for periodic BC
    if (P0(1)>=0 && P0(1)<=Lx && P0(2)>=0 && P0(2)<=Ly)
    Pn=[Pn;P0];
    Tn=[Tn;T];
    end
end


%%%% Plane 3
u = P2-P1;
w = P1 - V0;
D = n3*u';
N = -n3*w';
sI = N / D;
if (abs(D) < 10^-7 || sI < 0 || sI > 1)   % The segment is parallel to plane or The intersection point  lies outside the segment, so there is no intersection

else
    
    P0 = P1+ sI.*u;  %The intersection point
    T=[Lx 0 0];      %The translation vector for periodic BC
    if (P0(3)>=0 && P0(3)<=Lz && P0(2)>=0 && P0(2)<=Ly)
    Pn=[Pn;P0];
    Tn=[Tn;T];
    end
end


%%%% Plane 4
u = P2-P1;
w = P1 - V1;
D = n4*u';
N = -n4*w';
sI = N / D;
if (abs(D) < 10^-7 || sI < 0 || sI > 1)   % The segment is parallel to plane or The intersection point  lies outside the segment, so there is no intersection

else
    P0 = P1+ sI.*u;  %The intersection point
    T=[0 -Ly 0];      %The translation vector for periodic BC
    if (P0(1)>=0 && P0(1)<=Lx && P0(3)>=0 && P0(3)<=Lz)
    Pn=[Pn;P0];
    Tn=[Tn;T];
    end
end


%%%% Plane 5
u = P2-P1;
w = P1 - V1;
D = n5*u';
N = -n5*w';
sI = N / D;

if (abs(D) < 10^-7 || sI < 0 || sI > 1)   % The segment is parallel to plane or The intersection point  lies outside the segment, so there is no intersection

else
    P0 = P1+ sI.*u;  %The intersection point
    T=[-Lx 0 0];      %The translation vector for periodic BC
    if (P0(3)>=0 && P0(3)<=Lz && P0(2)>=0 && P0(2)<=Ly)
    Pn=[Pn;P0];
    Tn=[Tn;T];
    end
end


%%%% Plane 6
u = P2-P1;
w = P1 - V0;
D = n6*u';
N = -n6*w';
sI = N / D;

if (abs(D) < 10^-7 || sI < 0 || sI > 1)   % The segment is parallel to plane or The intersection point  lies outside the segment, so there is no intersection

else
    
    P0 = P1+ sI.*u;  %The intersection point
    T=[0 Ly 0];      %The translation vector for periodic BC
    if (P0(1)>=0 && P0(1)<=Lx && P0(3)>=0 && P0(3)<=Lz)
    Pn=[Pn;P0];
    Tn=[Tn;T];
    end
end

Pn1 = bsxfun(@minus, Pn, P1);
dn=(Pn1(:,1)).^2+(Pn1(:,2)).^2+(Pn1(:,3)).^2;
[~,d]=min(dn);
P=Pn(d,:);
T=Tn(d,:);
end












