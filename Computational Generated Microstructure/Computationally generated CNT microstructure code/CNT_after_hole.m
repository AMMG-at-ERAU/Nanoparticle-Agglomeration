function [Pout,center,L] = CNT_after_hole(RVE,P,C)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
P1 = P(1,:);
P2 = P(2,:);

P1n = P1 - C;
P2n = P2 - C;
d1 = sqrt(P1n*P1n'); %distance end point 1 to center
d2 = sqrt(P2n*P2n'); %distance end point 2 to center
if ((d1<=0.5*RVE.hole_D) && (d2<=0.5*RVE.hole_D))%both ends are inside the circle
    Pout = [];
    center = [];
    L = [];
else %find intersection segment ellipse
    x1 = P1(1); y1 = P1(3);
    x2 = P2(1); y2 = P2(3);
    xc = C(1); yc = C(3);
    dx = x2 - x1;
    dy = y2 - y1;
    
    A = dx^2 + dy^2;
    B = 2*(x1*dx + y1*dy - dx*xc - dy*yc);
    C = (x1 - xc)^2 + (y1 - yc)^2 - (0.5*RVE.hole_D)^2;
    
    det = B^2 - 4*A*C;
    if ((A<=1e-7) || (det <= 0)) % no intersection (det <0) or one intersection point (det = 0) (on or out of the segment)
        Pout = P;
    else   % two intersection points (on or off line segment)
        t1 = (-B - sqrt(det))/(2*A);
        t2 = (-B + sqrt(det))/(2*A);
        Pt1 = P1 + (P2 - P1)*t1;
        Pt2 = P1 + (P2 - P1)*t2;
        
        if (t1<0 || t1>1) && (t2<0 || t2>1) % both points out of segment
            Pout = P;
        elseif (t1>=0 && t1<=1) && (t2>=0 && t2<=1) %both points on segment
            if (t1 == 0)
                Pout = [Pt2;P2];
            elseif (t2 == 1)
                Pout = [P1;Pt1];
            else
                Pout = [P1;Pt1;Pt2;P2];
            end
        elseif (t1>=0 && t1<=1) %only first point on segment
            Pout = [P1;Pt1];
        else
            Pout = [Pt2;P2];
        end
    end
   Pout = Pout';
   n = size(Pout,2);
   L = Pout(:,1:2:n-1)-Pout(:,2:2:n);
   L = sqrt(sum(L.*L,1));
   center = 0.5*(Pout(:,1:2:n-1)+Pout(:,2:2:n));
   Pout = Pout';
   center = center';
end
end

