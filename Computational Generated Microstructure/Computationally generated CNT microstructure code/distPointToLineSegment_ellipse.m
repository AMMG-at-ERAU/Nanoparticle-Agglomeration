function [P1,P2,d]= distPointToLineSegment_ellipse( p1, p2,Center, pt, p7, a, b, alpha, dcut)
% r = distPointToLineSegment( xy0, xy1, xyP )
%%%%%xy0=[x;y]
P11=zeros(2,2);
d1=zeros(1,2);
for i=1:2
vx = p1(1)-pt(1,i);
vy = p1(2)-pt(2,i);
ux = p2(1)-p1(1);
uy = p2(2)-p1(2);
a1=p2-p1;
lenSqr= (ux*ux+uy*uy);
detP= -vx*ux + -vy*uy;
if( detP < 0 ) %%% point located before p1
    % calculate shortest distance between p1 and ellipse
    [P11(:,i), P2, d1(i)] = dist_point_ellipse( p1, p7,Center, a, b, alpha, dcut);
%   d1(i) = norm(p1-p7(:,i),2);
%   P11(:,i)=p1;
elseif( detP > lenSqr ) %%% point located after p2
    % calculate shortest distance between p2 and ellipse
    [P11(:,i), P2, d1(i)] = dist_point_ellipse( p2, p7,Center, a, b, alpha, dcut);
%   d1(i) = norm(p2-p7(:,i),2);
%   P11(:,i)=p2;
else %%% point between p1 and p2
  d1(i) = abs(ux*vy-uy*vx)/sqrt(lenSqr);
  P11(:,i)=(detP/lenSqr)*a1+p1;
  P2=Center;
end
if d1(i)<=dcut
   d=d1(i);
   P1=P11(:,i);
   P2=Center;
   break
end

end

% if any(d1<dcut)
%     P2=p7;
% if d1(1)<d1(2)
%    d=d1(1);
%    P1=P11(:,1);
% else
%    d=d1(2);
%    P1=P11(:,2);
% end
% else
%    P1=[NaN;NaN];
%    P2=P1;
%    d=min(d1);
% end
end



