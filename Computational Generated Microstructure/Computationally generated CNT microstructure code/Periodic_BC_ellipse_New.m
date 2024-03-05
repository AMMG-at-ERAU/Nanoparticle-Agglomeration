function [ Pg,Pc,cent,Ar ] = Periodic_BC_ellipse_New( Lx, Ly, p7, a, b, alpha )
%PBC with segment of ellipse
%if it had gone out on the corner, we should count the corner ones only
%once. For that we flag that case with ex=0
%Input
%Lx, Ly : respectively the x and y length of the RVE
%a, b, alpha : respectively the major, minor axis and angle of the original ellipse
%p7 : the center of the original ellipse define as a column vector
%Output
%Pg  : the centers of all the ellipses created due to PBC, the original ellipse
%center is included. Each colum for the coordinates of each center
%([xg;yg])
ex=1;
% xg=p7(1);
% yg=p7(2);
Pg=[p7];
Pc=[];
Ar=[];

% Parameters of ellipse and ellipse
h = p7(1);
k = p7(2);
A = ( ( (cos(alpha))^2 )./ ( a.^2 )  + ( (sin(alpha))^2 )./ ( b.^2 ) );
B = - 2 * cos (alpha) * sin(alpha) * (-1./(a.^2) + 1./(b.^2));
C = ( ((sin(alpha))^2)./ (a.^2)  + ((cos(alpha))^2)./ (b.^2) );


t = ellipse_boundary( Lx, Ly, p7, a, b, alpha );%intersection point 
if all(isnan(t))
    Pc=[p7(1)-a*cos(alpha) p7(1)+a*cos(alpha);p7(2)-a*sin(alpha) p7(2)+a*sin(alpha)];
    cent=p7;
    Ar=[pi*a*b];
elseif  all(t(:,1)>=0) & all(t(:,1)<=1)%on line segment
    if all(t(:,3)<0) %case b1
        Pg=[Pg [h;k+Ly] [h+Lx;k+Ly]];
        
        %Area clculation
        xi=Lx*t(1,1)-h;
        xii=Lx*t(2,1)-h;
        yi=0-k;
        xu=0-h;
        yu=Ly*t(1,3)-k;
        yuu=Ly*t(2,3)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[pi*a*b-area area-area1 area1];
        %PBC
        t2=0.5*(t(1,1)+t(2,1));
        x2=t2*Lx;
        y2=0;
        t1=0.5*(t(1,3)+t(2,3));
        x1=0;
        y1=t1*Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x2 x22;y2 y22] [x1 x2;y1+Ly y2+Ly] [x11+Lx x1+Lx;y11+Ly y1+Ly]];
    elseif all(t(:,3)>0) %case b5
        Pg=[Pg [h;k+Ly] [h+Lx;k]];
        
        %Area clculation
        xi=Lx*t(1,1)-h;
        xii=Lx*t(2,1)-h;
        yi=0-k;
        xu=0-h;
        yu=Ly*t(1,3)-k;
        yuu=Ly*t(2,3)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end        
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[area1-area area pi*a*b-area1];
        %PBC
        
        t1=0.5*(t(1,3)+t(2,3));
        x1=0;
        y1=t1*Ly;
        t2=0.5*(t(1,1)+t(2,1));
        x2=t2*Lx;
        y2=0;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x1 x2;y1 y2] [x2 x22;y2+Ly y22+Ly] [x11+Lx x1+Lx;y11 y1]];
    elseif all(t(:,4)<0) %case b2
        Pg=[Pg [h;k+Ly] [h-Lx;k+Ly]];
        %Area clculation
        xi=Lx*t(1,1)-h;
        xii=Lx*t(2,1)-h;
        yi=0-k;
        xu=Lx-h;
        yu=Ly*t(1,4)-k;
        yuu=Ly*t(2,4)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[pi*a*b-area area-area1 area1];
        %PBC
        t1=0.5*(t(1,1)+t(2,1));
        x1=t1*Lx;
        y1=0;
        t2=0.5*(t(1,4)+t(2,4));
        x2=Lx;
        y2=t2*Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x11 x1;y11 0] [x1 Lx;Ly y2+Ly] [0 x22-Lx;y2+Ly y22+Ly]]; 
    elseif all(t(:,4)>0) %case b6
        Pg=[Pg [h-Lx;k] [h;k+Ly]];
        %Area clculation
        xi=Lx*t(1,1)-h;
        xii=Lx*t(2,1)-h;
        yi=0-k;
        xu=Lx-h;
        yu=Ly*t(1,4)-k;
        yuu=Ly*t(2,4)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);
        if al-al1>0
            al=2*pi-al+al1;
            al1=0;
        end         
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[area1-area pi*a*b-area1 area];
        %PBC
        t1=0.5*(t(1,1)+t(2,1));
        x1=t1*Lx;
        y1=0;
        t2=0.5*(t(1,4)+t(2,4));
        x2=Lx;
        y2=t2*Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x1 Lx;0 y2] [0 x22-Lx;y2 y22] [x11 x1;y11+Ly Ly]];
    else %case a
        Pg=[Pg [h;k+Ly]];
        %Area clculation
        xi=Lx*t(1,1)-h;
        xii=Lx*t(2,1)-h;
        yi=0-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[pi*a*b-area area];
        %PBC
        t1=0.5*(t(1,1)+t(2,1));
        x1=t1*Lx;
        y1=0;
        x2=x1;
        y2=b/10;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x1 x22;0 y22] [x11 x1;y11+Ly Ly]];        
        
    end
elseif all(t(:,2)>=0) & all(t(:,2)<=1)%on line segment
    if all(t(:,3)>1) %case c2
        Pg=[Pg [h;k-Ly] [h+Lx;k-Ly]];
        
        %Area clculation
        xi=Lx*t(1,2)-h;
        xii=Lx*t(2,2)-h;
        yi=Ly-k;
        xu=0-h;
        yu=Ly*t(1,3)-k;
        yuu=Ly*t(2,3)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[area area1-area pi*a*b-area1];
        %PBC        
        t2=0.5*(t(1,2)+t(2,2));
        x2=t2*Lx;
        y2=Ly;
        t1=0.5*(t(1,3)+t(2,3));
        x1=0;
        y1=t1*Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x2 x22;y2 y22] [0 x2;y1-Ly 0] [x11+Lx Lx;y11-Ly y1-Ly]];
    elseif all(t(:,3)<1) %case c4
        Pg=[Pg [h;k-Ly] [h+Lx;k]];
        %Area clculation
        xi=Lx*t(1,2)-h;
        xii=Lx*t(2,2)-h;
        yi=Ly-k;
        xu=0-h;
        yu=Ly*t(1,3)-k;
        yuu=Ly*t(2,3)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[area-area1 pi*a*b-area area1];
        %PBC 
        t1=0.5*(t(1,3)+t(2,3));
        x1=0;
        y1=t1*Ly;
        t2=0.5*(t(1,2)+t(2,2));
        x2=t2*Lx;
        y2=Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [0 x2;y1 Ly] [x2 x22;0 y22-Ly] [x11+Lx Lx;y11 y1]];
    elseif all(t(:,4)>1) %case c3
        Pg=[Pg [h;k-Ly] [h-Lx;k-Ly]];
        %Area clculation
        xi=Lx*t(1,2)-h;
        xii=Lx*t(2,2)-h;
        yi=Ly-k;
        xu=Lx-h;
        yu=Ly*t(1,4)-k;
        yuu=Ly*t(2,4)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);
        if al-al1>0
            al=2*pi-al+al1;
            al1=0;
        end        
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[area area1-area pi*a*b-area1];
        %PBC         
        t1=0.5*(t(1,2)+t(2,2));
        x1=t1*Lx;
        y1=Ly;
        t2=0.5*(t(1,4)+t(2,4));
        x2=Lx;
        y2=t2*Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x11 x1;y11 Ly] [x1 Lx;0 y2-Ly] [0 x22-Lx;y2-Ly y22-Ly]]; 
    elseif all(t(:,4)<1) %case c5
        Pg=[Pg [h-Lx;k] [h;k-Ly]];
        %Area clculation
        xi=Lx*t(1,2)-h;
        xii=Lx*t(2,2)-h;
        yi=Ly-k;
        xu=Lx-h;
        yu=Ly*t(1,4)-k;
        yuu=Ly*t(2,4)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end         
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);        
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[area-area1 area1 pi*a*b-area];
        %PBC         
        t1=0.5*(t(1,2)+t(2,2));
        x1=t1*Lx;
        y1=Ly;
        t2=0.5*(t(1,4)+t(2,4));
        x2=Lx;
        y2=t2*Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x1 Lx;Ly y2] [0 x22-Lx;y2 y22] [x11 x1;y11-Ly 0]];
    else %case c1
        Pg=[Pg [h;k-Ly]];
        %Area clculation
        xi=Lx*t(1,2)-h;
        xii=Lx*t(2,2)-h;
        yi=Ly-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end        
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[area pi*a*b-area];
        %PBC        
        t1=0.5*(t(1,2)+t(2,2));
        x1=t1*Lx;
        y1=Ly;
        x2=x1;
        y2=b/10+Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x11 x1;y11 Ly] [x1 x22;0 y22-Ly]];        
        
    end        
elseif all(t(:,3)>=0) & all(t(:,3)<=1)%on line segment
    if all(t(:,2)<0) %case d2
        Pg=[Pg [h+Lx;k] [h+Lx;k-Ly]];
        %Area clculation
        xi=Lx*t(1,2)-h;
        xii=Lx*t(2,2)-h;
        yi=Ly-k;
        xu=0-h;
        yu=Ly*t(1,3)-k;
        yuu=Ly*t(2,3)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[area1 area-area1 pi*a*b-area];
        %PBC         
        t1=0.5*(t(1,2)+t(2,2));
        x1=t1*Lx;
        y1=Ly;
        t2=0.5*(t(1,3)+t(2,3));
        x2=0;
        y2=t2*Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [0 x22;y2 y22] [x1+Lx Lx;Ly y2] [x11+Lx x1+Lx;y11-Ly 0]];
    elseif all(t(:,1)<0) %case d3
        Pg=[Pg [h+Lx;k] [h+Lx;k+Ly]];
        %Area clculation
        xi=Lx*t(1,1)-h;
        xii=Lx*t(2,1)-h;
        yi=0-k;
        xu=0-h;
        yu=Ly*t(1,3)-k;
        yuu=Ly*t(2,3)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);       
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[pi*a*b-area1 area1-area area];
        %PBC        
        t1=0.5*(t(1,1)+t(2,1));
        x1=t1*Lx;
        y1=0;
        t2=0.5*(t(1,3)+t(2,3));
        x2=0;
        y2=t2*Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [0 x22;y2 y22] [x1+Lx Lx;0 y2] [x11+Lx x1+Lx;y11+Ly Ly]];
    else %case d1
        Pg=[Pg [h+Lx;k]];
        %Area clculation
        xi=0-h;
        yi=Ly*t(1,3)-k;
        yii=Ly*t(2,3)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xi*cos(alpha)+yii*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yii*cos(alpha)-xi*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        if al-al1>0
            al=2*pi-al+al1;
            al1=0;
        end        
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[pi*a*b-area area];
        %PBC         
        t1=0.5*(t(1,3)+t(2,3));
        x1=0;
        y1=t1*Ly;
        x2=b/10;
        y2=y1;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [0 x22;y1 y22] [x11+Lx Lx;y11 y1]];        
        
    end    
elseif all(t(:,4)>=0) & all(t(:,4)<=1)%on line segment
    if all(t(:,1)>1) %case e3
        Pg=[Pg [h-Lx;k] [h-Lx;k+Ly]];
        %Area clculation
        xi=Lx*t(1,1)-h;
        xii=Lx*t(2,1)-h;
        yi=0-k;
        xu=Lx-h;
        yu=Ly*t(1,4)-k;
        yuu=Ly*t(2,4)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end              
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);    
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end          
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[pi*a*b-area1 area1-area area];
        %PBC         
        t1=0.5*(t(1,4)+t(2,4));
        x1=Lx;
        y1=t1*Ly;
        t2=0.5*(t(1,1)+t(2,1));
        x2=t2*Lx;
        y2=0;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x11 Lx;y11 y1] [0 x2-Lx;y1 0] [x2-Lx x22-Lx;Ly y22+Ly]];
    elseif all(t(:,2)>1) %case e2
        Pg=[Pg [h-Lx;k] [h-Lx;k-Ly]];
        %Area clculation
        xi=Lx*t(1,2)-h;
        xii=Lx*t(2,2)-h;
        yi=Ly-k;
        xu=Lx-h;
        yu=Ly*t(1,4)-k;
        yuu=Ly*t(2,4)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end              
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        al=mod(atan2(yit(3),xit(3)),2*pi);
        al1=mod(atan2(yit(4),xit(4)),2*pi);       
        area1=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[area1 area-area1 pi*a*b-area];
        %PBC        
        t1=0.5*(t(1,4)+t(2,4));
        x1=Lx;
        y1=t1*Ly;
        t2=0.5*(t(1,2)+t(2,2));
        x2=t2*Lx;
        y2=Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x11 Lx;y11 y1] [0 x2-Lx;y1 Ly] [x2-Lx x22-Lx;0 y22-Ly]];
    else %case e1
        Pg=[Pg [h-Lx;k]];
        %Area clculation
        xi=Lx-h;
        yi=Ly*t(1,4)-k;
        yii=Ly*t(2,4)-k;
        xit=[xi*cos(alpha)+yi*sin(alpha) xi*cos(alpha)+yii*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yii*cos(alpha)-xi*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);
        if al-al1>0
            al=2*pi-al+al1;
            al1=0;
        end        
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        Ar=[area pi*a*b-area];
        %PBC        
        t1=0.5*(t(1,4)+t(2,4));
        x1=Lx;
        y1=t1*Ly;
        x2=Lx+b/10;
        y2=y1;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x11 Lx;y11 y1] [0 x22-Lx;y1 y22]];        
        
    end   
else
    if any(t(:,1)<0) & any(t(:,3)<0) %case b3
        Pg=[Pg [h+Lx;k+Ly] [h+Lx;k] [h;k+Ly] ];
        %Area clculation
        xi=Lx*t(1,1)-h;
        xii=Lx*t(2,1)-h;
        yi=0-k;
        xu=0-h;
        yu=Ly*t(1,3)-k;
        yuu=Ly*t(2,3)-k;
        areat=0.5*(xi+h)*(yu+k);
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);  
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end 
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        %areac=pi*a*b-area;
        al=mod(atan2(yit(1),xit(1)),2*pi); %to not change the original value        
        alll=mod(atan2(yit(3),xit(3)),2*pi);
        all1=mod(atan2(yit(4),xit(4)),2*pi);   
        if alll-all1>0
            alll=2*pi-alll+all1;
            all1=0;
        end        
        area1=0.5*a*b*abs((alll-all1)-sin(alll-all1));
        %area1c=pi*a*b-area1;
        alll=mod(atan2(yit(3),xit(3)),2*pi);        
        if al-alll<0
            al=2*pi+al-alll;
            alll=0;
        end         
        area2=0.5*a*b*abs((al-alll)-sin(al-alll));
        %area2c=pi*a*b-area2;
        Ar=[pi*a*b-area2+areat area+areat+area1-area2 area2-areat-area area2-areat-area1];
        %PBC                
        x1=0;
        y1=0;
        x2=Lx;
        y2=Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [0 x22;0 y22] [x11+Lx Lx;y11+Ly Ly]];
        
        x1=0;
        y1=0;
        x2=-Lx;
        y2=Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x22+Lx Lx;y22 0] [0 x11;Ly y11+Ly]];
    elseif any(t(:,1)>1) & any( t(:,4)<0) %case b4 
        Pg=[Pg [h-Lx;k+Ly] [h-Lx;k] [h;k+Ly] ];
        %Area clculation
        xi=Lx*t(1,1)-h;
        xii=Lx*t(2,1)-h;
        yi=0-k;
        xu=Lx-h;
        yu=Ly*t(1,4)-k;
        yuu=Ly*t(2,4)-k;
        areat=0.5*(Lx-xii-h)*(yu+k);
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);  
        if al-al1<0
            al=2*pi+al-al1;
            al1=0;
        end 
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        %areac=pi*a*b-area;
        al1=mod(atan2(yit(2),xit(2)),2*pi); %to not change the original value        
        alll=mod(atan2(yit(3),xit(3)),2*pi);
        all1=mod(atan2(yit(4),xit(4)),2*pi);   
        if alll-all1<0
            alll=2*pi+alll-all1;
            all1=0;
        end        
        area1=0.5*a*b*abs((alll-all1)-sin(alll-all1));
        %area1c=pi*a*b-area1;
        alll=mod(atan2(yit(3),xit(3)),2*pi);        
        if al1-alll>0
            al1=2*pi-al1+alll;
            alll=0;
        end         
        area2=0.5*a*b*abs((al1-alll)-sin(al1-alll));
        %area2c=pi*a*b-area2;
        Ar=[pi*a*b-area2+areat area+areat+area1-area2 area2-areat-area area2-areat-area1];
        %PBC       
        x1=Lx;
        y1=0;
        x2=0;
        y2=Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;

        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x22 Lx;y22 0] [0 x11-Lx;Ly y11+Ly]];
        
        x1=Lx;
        y1=0;
        x2=2*Lx;
        y2=Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [0 x22-Lx;0 y22] [x11 Lx;y11+Ly Ly]];
    elseif any(t(:,2)<0) & any(t(:,3)>1) %case c7 
        Pg=[Pg [h+Lx;k-Ly] [h+Lx;k] [h;k-Ly] ];
        %Area clculation
        xi=Lx*t(1,2)-h;
        xii=Lx*t(2,2)-h;
        yi=Ly-k;
        xu=0-h;
        yu=Ly*t(1,3)-k;
        yuu=Ly*t(2,3)-k;
        areat=0.5*(xi+h)*(Ly-yuu-k);
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);  
        if al-al1>0
            al=2*pi-al+al1;
            al1=0;
        end 
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        %areac=pi*a*b-area;
        al=mod(atan2(yit(1),xit(1)),2*pi); %to not change the original value        
        alll=mod(atan2(yit(3),xit(3)),2*pi);
        all1=mod(atan2(yit(4),xit(4)),2*pi);   
        if alll-all1>0
            alll=2*pi-alll+all1;
            all1=0;
        end        
        area1=0.5*a*b*abs((alll-all1)-sin(alll-all1));
        %area1c=pi*a*b-area1;
        all1=mod(atan2(yit(4),xit(4)),2*pi);         
        if al-all1>0
            al=2*pi-al+all1;
            all1=0;
        end         
        area2=0.5*a*b*abs((al-all1)-sin(al-all1));
        %area2c=pi*a*b-area2;
        Ar=[pi*a*b-area2+areat area+areat+area1-area2 area2-areat-area area2-areat-area1];
        %PBC              
        x1=0;
        y1=Ly;
        x2=Lx;
        y2=0;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [0 x22;Ly y22] [x11+Lx Lx;y11-Ly 0]];
        
        x1=0;
        y1=Ly;
        x2=Lx;
        y2=2*Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x11+Lx Lx;y11 Ly] [0 x22;0 y22-Ly]]; 
    else %case c6
        Pg=[Pg [h-Lx;k-Ly] [h-Lx;k] [h;k-Ly] ];
        %Area clculation
        xi=Lx*t(1,2)-h;
        xii=Lx*t(2,2)-h;
        yi=Ly-k;
        xu=Lx-h;
        yu=Ly*t(1,4)-k;
        yuu=Ly*t(2,4)-k;
        areat=0.5*(Lx-xii-h)*(Ly-yuu-k);
        xit=[xi*cos(alpha)+yi*sin(alpha) xii*cos(alpha)+yi*sin(alpha) xu*cos(alpha)+yu*sin(alpha) xu*cos(alpha)+yuu*sin(alpha)];%in ellipse axis coordinaates
        yit=(a/b)*[yi*cos(alpha)-xi*sin(alpha) yi*cos(alpha)-xii*sin(alpha) yu*cos(alpha)-xu*sin(alpha) yuu*cos(alpha)-xu*sin(alpha)];%in ellipse axis coordinaate
        al=mod(atan2(yit(1),xit(1)),2*pi);
        al1=mod(atan2(yit(2),xit(2)),2*pi);  
        if al-al1>0
            al=2*pi-al+al1;
            al1=0;
        end 
        area=0.5*a*b*abs((al-al1)-sin(al-al1));
        %areac=pi*a*b-area;
        al1=mod(atan2(yit(2),xit(2)),2*pi); %to not change the original value        
        alll=mod(atan2(yit(3),xit(3)),2*pi);
        all1=mod(atan2(yit(4),xit(4)),2*pi);   
        if alll-all1<0
            alll=2*pi+alll-all1;
            all1=0;
        end        
        area1=0.5*a*b*abs((alll-all1)-sin(alll-all1));
        %area1c=pi*a*b-area1;
        all1=mod(atan2(yit(4),xit(4)),2*pi);         
        if al1-all1<0
            al1=2*pi+al1-all1;
            all1=0;
        end         
        area2=0.5*a*b*abs((al1-all1)-sin(al1-all1));
        %area2c=pi*a*b-area2;
        Ar=[pi*a*b-area2+areat area+areat+area1-area2 area2-areat-area area2-areat-area1];
        %PBC                    
        x1=Lx;
        y1=Ly;
        x2=2*Lx;
        y2=2*Ly;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [x11 Lx;y11 Ly] [0 x22-Lx;0 y22-Ly]];
        x1=Lx;
        y1=Ly;
        x2=2*Lx;
        y2=0;
        % Parameters of quadratic equation D*t^2+E*t+F = 0
        D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
        E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
        F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;
        
        det = E^2-4*F*D;  %determinant
        t1=(-E-sqrt(det))/(2*D);
        x11=x1+(x2-x1)*t1;
        y11=y1+(y2-y1)*t1;
        t2=(-E+sqrt(det))/(2*D);
        x22=x1+(x2-x1)*t2;
        y22=y1+(y2-y1)*t2;
        Pc=[Pc [0 x22-Lx;Ly y22] [x11 Lx;y11-Ly 0]];

    end            
end
cent=0.5*(Pc(:,1:2:end)+Pc(:,2:2:end));

end