function [ P1, P2, distance ] = SegmentEllipse( p1, p2, p5, p6, p7, Center, a, b, alpha, dcut, d_tun,D1 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% ellipse is  defined by the equation
%(x*cos(alpha) + y * sin(alpha))^2 / a^2 + (x*sin(alpha) - y *
%        cos(alpha))^2 / b^2 = 1
% the center p7=[h;k]
% the radius along x-axis a and the radius along y-axis b
% the counterclockwise angle alpha (in radian)
%line is defined by the point column vectors p1 and p2

% Copyright: Audrey Gbaguidi
% Parameters of line
x1=p1(1);
y1=p1(2);
x2=p2(1);
y2=p2(2);
% Parameters of ellipse and ellipse with max tunnelling
ai=[a+0.5*D1,a+dcut];
bi=[b+0.5*D1,b+dcut];
h = p7(1);
k = p7(2);
A = ( ( (cos(alpha))^2 )./ ( ai.^2 )  + ( (sin(alpha))^2 )./ ( bi.^2 ) );
B = - 2 * cos (alpha) * sin(alpha) * (-1./(ai.^2) + 1./(bi.^2));
C = ( ((sin(alpha))^2)./ (ai.^2)  + ((cos(alpha))^2)./ (bi.^2) );

% Parameters of quadratic equation D*t^2+E*t+F = 0
D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;

det = E.^2-4*F.*D;  %determinant

   ti=[(-E(2)+sqrt(det(2)))/(2*D(2)) (-E(2)-sqrt(det(2)))/(2*D(2))];

   
if ( det(2)>=0 && ( any(ti>=0 & ti<=1) || (ti(1)<0 && ti(2)>1) || (ti(1)>1 && ti(2)<0) ) ) %CNT segment is intersecting or inside max tunnelling ellipse

if det(1)>=0 %CNT line is intersecting the ellipse
   
   t=[(-E(1)+sqrt(det(1)))/(2*D(1)) (-E(1)-sqrt(det(1)))/(2*D(1))];

%     P11 = p1+(p2-p1)*t(1); % intersection point 1
%     P21 = p1+(p2-p1)*t(2); % intersection point 2
   
   if any(t>=0 & t<=1) | (t(1)<0 & t(2)>1) | (t(1)>1 & t(2)<0) %CNT segment contact with ellipse or CNT segment inside ellipse
      [P1,P2]=DistBetween2Segment_modified_pt(p1,p2,p5,p6,d_tun); %intersection between segment and ellipse segment
      distance=d_tun;
   elseif any(t<0) %distance between p1 and ellipse if CNT is outside ellipse (tunnelling)
     [P1, P2, distance] = dist_point_ellipse( p1, p7,Center, a, b, alpha, dcut);
   elseif any(t>1) %distance between p2 and ellipse is CNT is outside ellipse (tunnelling)
     [P1, P2, distance] = dist_point_ellipse( p2, p7,Center, a, b, alpha, dcut);
   end
else %if line never intercepts the ellipse
%look for the two lines parralel to the CNT segment and intersecting
%the ellipse on the points pi1 and pi2. FInd the pi1 and pi2
m=(y1-y2)/(x1-x2); %slope of the CNT segment
%P2=p7; %since there is no intersection the tunelling point on the ellipse is the center
if (abs(m)==inf) %verical line

    r=2*C(1)/B(1);
    yi1=k+1/sqrt(-C(1)+A(1)*r^2);
    yi2=k-1/sqrt(-C(1)+A(1)*r^2);
    xi1=h-r*(1/sqrt(-C(1)+A(1)*r^2));
    xi2=h-r*(-1/sqrt(-C(1)+A(1)*r^2));
%     P11=[xi1;yi1];
%     P21=[xi2;yi2];
    [P1, P2, distance]= distPointToLineSegment_ellipse( p1, p2, Center, [xi1 xi2; yi1 yi2], p7, a, b, alpha, dcut);
else
    r=(2*m*C(1)+B(1))/(2*A(1)+m*B(1));
    yi1=k+1/sqrt(A(1)*r^2-r*B(1)+C(1));
    yi2=k-1/sqrt(A(1)*r^2-r*B(1)+C(1));
    xi1=h-r*(yi1-k);
    xi2=h-r*(yi2-k);
%     P11=[xi1;yi1];
%     P21=[xi2;yi2];
    [P1, P2, distance]= distPointToLineSegment_ellipse( p1, p2, Center, [xi1 xi2; yi1 yi2], p7, a, b, alpha, dcut);
end

end 
if distance>dcut %make sure distance is always less than dcut due to dtun with strain
    P1=[NaN;NaN];
end
else
    P1=[NaN;NaN];
    P2=[NaN;NaN];
    distance=inf;
end
%    %%%check
% x1=P1(1);
% y1=P1(2);
% x2=P2(1);
% y2=P2(2);
% A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2
% A*(x2-h)^2+B*(x2-h)*(y2-k)+C*(y2-k)^2
% %

end