function [ t ] = ellipse_boundary( Lx, Ly, p7, a, b, alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% ellipse is  defined by the equation
%(x*cos(alpha) + y * sin(alpha))^2 / a^2 + (x*sin(alpha) - y *
%        cos(alpha))^2 / b^2 = 1
% the center p7=[h;k]
% the radius along x-axis a and the radius along y-axis b
% the counterclockwise angle alpha (in radian)
%line is defined by the point column vectors p1 (point closer to the
%perpendicular axis) and p2 (the one farther)
% t(:,1) intersection points with y = 0
% t(:,2) intersection points with y = Ly
% t(:,3) intersection points with x = 0
% t(:,4) intersection points with x = Lx

% Copyright: Audrey Gbaguidi
t=zeros(2,4);
% Parameters of line y=0
x1=0;
y1=0;
x2=Lx;
y2=0;
% Parameters of ellipse and ellipse
h = p7(1);
k = p7(2);
A = ( ( (cos(alpha))^2 )./ ( a.^2 )  + ( (sin(alpha))^2 )./ ( b.^2 ) );
B = - 2 * cos (alpha) * sin(alpha) * (-1./(a.^2) + 1./(b.^2));
C = ( ((sin(alpha))^2)./ (a.^2)  + ((cos(alpha))^2)./ (b.^2) );

% Parameters of quadratic equation D*t^2+E*t+F = 0
D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;

det = E^2-4*F*D;  %determinant
if det>10^-6 % two intersecting points
t(:,1)=[(-E+sqrt(det))/(2*D);(-E-sqrt(det))/(2*D)];
else
t(:,1)=[NaN;NaN];
end

% Parameters of line y=Ly
x1=0;
y1=Ly;
x2=Lx;
y2=Ly;
% Parameters of quadratic equation D*t^2+E*t+F = 0
D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;

det = E^2-4*F*D;  %determinant
if det>10^-6 % two intersecting points
t(:,2)=[(-E+sqrt(det))/(2*D);(-E-sqrt(det))/(2*D)];
else
t(:,2)=[NaN;NaN];
end

% Parameters of line x=0
x1=0;
y1=0;
x2=0;
y2=Ly;
% Parameters of quadratic equation D*t^2+E*t+F = 0
D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;

det = E^2-4*F*D;  %determinant
if det>10^-6 % two intersecting points
t(:,3)=[(-E+sqrt(det))/(2*D);(-E-sqrt(det))/(2*D)];
else
t(:,3)=[NaN;NaN];
end


% Parameters of line x=Lx
x1=Lx;
y1=0;
x2=Lx;
y2=Ly;
% Parameters of quadratic equation D*t^2+E*t+F = 0
D = A*(x2-x1)^2+B*(x2-x1)*(y2-y1)+C*(y2-y1)^2;
E = 2*A*(x1-h)*(x2-x1)+B*((x1-h)*(y2-y1)+(x2-x1)*(y1-k))+2*C*(y1-k)*(y2-y1);
F = A*(x1-h)^2+B*(x1-h)*(y1-k)+C*(y1-k)^2-1;

det = E^2-4*F*D;  %determinant
if det>10^-6 % two intersecting points
t(:,4)=[(-E+sqrt(det))/(2*D);(-E-sqrt(det))/(2*D)];
else
t(:,4)=[NaN;NaN];
end
end



