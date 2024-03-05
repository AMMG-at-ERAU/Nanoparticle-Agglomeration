% Computes the minimum distance between two line segments. Code
% is adapted for Matlab from Dan Sunday's Geometry Algorithms originally
% written in C++
% http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment

% Usage: Input the start and end x,y,z coordinates for two line segments. 
% p1, p2 are [x,y,z] coordinates of first line segment and p3,p4 are for
% second line segment. 

% Output: scalar minimum distance between the two segments.

%  Example:
%	P1 = [0 0 0];     P2 = [1 0 0];
%   P3 = [0 1 0];     P4 = [1 1 0];
%	dist = DistBetween2Segment(P1, P2, P3, P4)


function [distance, dti, dta] = A_DistSeg2Seg(p1, p2, p3, p4,dvdw,Dia)
% distance of a CNT to n neghboring CNTs                   tunnelling distance with CNT diameter substracted
% dti           distance between intersection on p1p2 and p1
% dta           distance between intersection on p3p4 and p3
% p1 1x3 matrix starting point of CNT i
% p2 1x3 matrix end point of CNT i
% p3 nx3 matrix starting point of neighboring CNTs
% p4 nx3 matrix end point of neighboring CNTs
format long
    n = size(p3,1);
    pt = zeros(n,12);
    u = ones(n,1)*(p1 - p2); %nx3 matrix
    v = p3 - p4;             %nx3 matrix
    w = ones(n,1)*p2 - p4;   %nx3 matrix
%     w1 = ones(n,1)*p2 - p3;
%     w2 = ones(n,1)*p1 - p3;
%     w3 = ones(n,1)*p1 - p4;
    
    pt(:,1) = sum(u.*u,2); %a is the 1st column
    pt(:,2) = sum(u.*v,2); %b is the 2nd column
    pt(:,3) = sum(v.*v,2); %c is the 3rd column
    pt(:,4) = sum(u.*w,2); %d is the 4th column
    pt(:,5) = sum(v.*w,2); %e is the 5th column
    pt(:,6) = pt(:,1).*pt(:,3) - pt(:,2).*pt(:,2); %D is the 6th column
    pt(:,7) = pt(:,6); %sD is the 7th column
    pt(:,8) = pt(:,6); %tD is the 8th column
    
    SMALL_NUM = 0.0000000001;
    % compute the line parameters of the two closest points 
    t = pt(:,6) < SMALL_NUM; %first condition D < SMALL NUM
    pt(t,9) = 0.0; %sN is the 9th column
    pt(t,7) = 1.0;
    pt(t,10) = pt(t,5);  %tN is the 10th column
    pt(t,8) = pt(t,3);
    
    pt(~t,9) = pt(~t,2).*pt(~t,5) - pt(~t,3).*pt(~t,4);
    pt(~t,10) = pt(~t,1).*pt(~t,5) - pt(~t,2).*pt(~t,4);
    
    t1 = pt(:,9) < 0.0;
    t2 = pt(:,9) > pt(:,7);
    pt(~t & t1,9) = 0.0;
    pt(~t & t1,10) = pt(~t & t1,5);
    pt(~t & t1,8) = pt(~t & t1,3);
    
    pt(~t & (~t1 & t2),9) = pt(~t & (~t1 & t2),7);
    pt(~t & (~t1 & t2),10) = pt(~t & (~t1 & t2),5) + pt(~t & (~t1 & t2),2);
    pt(~t & (~t1 & t2),8) = pt(~t & (~t1 & t2),3);
    
    t = pt(:,10) < 0.0;
    t1 = pt(:,10) > pt(:,8);
    
    pt(t,10) = 0.0;
    t2 = -pt(:,4) < 0.0;
    t3 = -pt(:,4) > pt(:,1);
    pt(t & t2,9) = 0.0;
    pt(t & (~t2 & t3),9) = pt(t & (~t2 & t3),7);
    pt(t & (~t2 & ~t3),9) = -pt(t & (~t2 & ~t3),4);
    pt(t & (~t2 & ~t3),7) = pt(t & (~t2 & ~t3),1);
    
    pt(~t & t1,10) = pt(~t & t1,8);
    t2 = (-pt(:,4) + pt(:,2)) < 0.0;
    t3 = (-pt(:,4) + pt(:,2)) > pt(:,1);   
    pt((~t & t1) & t2,9) = 0.0;
    pt((~t & t1) & (~t2 & t3),9) = pt((~t & t1) & (~t2 & t3),7);
    pt((~t & t1) & (~t2 & ~t3),9) = -pt((~t & t1) & (~t2 & ~t3),4) + pt((~t & t1) & (~t2 & ~t3),2);
    pt((~t & t1) & (~t2 & ~t3),7) = pt((~t & t1) & (~t2 & ~t3),1);
     
    % finally do the division to get sc and tc
    t = abs(pt(:,9)) < SMALL_NUM;
    pt(t,11) = 0.0; %sc is the 11th column
    pt(~t,11) = pt(~t,9)./pt(~t,7);
    
    t = abs(pt(:,10)) < SMALL_NUM;
    pt(t,12) = 0.0; %tc is the 11th column
    pt(~t,12) = pt(~t,10)./pt(~t,8);    
    
    % get the difference of the two closest points
    dP = w + ((pt(:,11)*ones(1,3)) .* u) - ((pt(:,12)*ones(1,3)) .* v);
    distance = sqrt(sum(dP.*dP,2)) - Dia; %substract CNT diameter
    %distance = min([sqrt(dot(dP,dP,2)) - Dia, sqrt(dot(w,w,2)) - Dia , sqrt(dot(w1,w1,2))- Dia , sqrt(dot(w2,w2,2))- Dia , sqrt(dot(w3,w3,2))- Dia],[],2); %substract CNT diameter
    v1 = ones(n,1)*p2 + ((pt(:,11)*ones(1,3)) .* u); % Closest point on object p1p2 
    v1 = ones(n,1) * p1 - v1;            % vector from p1 to the intersection on p1p2
    v2 = p4 + ((pt(:,12)*ones(1,3)) .* v);           % Closest point on object p3p4
    v2 = p3 - v2;                        % vector from p3 to intersec point on p3p4
    dti = sqrt(sum(v1.*v1,2))+1e-5;%0.1;        %I have added 0.1 so that the value won't be zer0 when the intersection is p1 or p3. I will subtract it
    dta = sqrt(sum(v2.*v2,2))+1e-5;%0.1;        %every time before using the value

    %distance(distance <=1e-9) = -10; %Use this for A_ResNetwork_Soft_leafs.m (So there can be between CNT resistor)
    %distance(distance <dvdw) = -10;
    %distance(distance <dvdw) = dvdw; %Use this for A_ResNetwork_Soft_leafs_notouch.m (So there is no contact between CNT resistor)
    distance(distance <=0) = 1e-6;

    
end