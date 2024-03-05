function [ dCNT, ja ] = Hybrid_contact_samarth( Resistor_GNP, Resistor_Center_GNP, Center_GNP, rad_GNP, RVE)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
% ///////////////////////// This is working 2nd FAST///////////////////////////
% PP has all the CNT ends points as each row with 3 columns
% Pcenter is a 3 columns matrix with all centers
% Plength is a row vector with the CNT lengths
% dCNT is a N by N matrix with the tunneling disances
% ja is a N by N matrix with each row representing the different junction
% on the CNT of that row (their distance from the starting point of the
% CNT)
% N is the number of CNTs
% //////////////////////// Here A_DistSeg2Seg is used. It computes simultaneously the distance for all
% adjacents CNT to the current CNT /////////////////////////
Ncnt = 0;
Ngnp = size(Center_GNP,1);
dCNT = sparse(Ncnt + Ngnp,Ncnt + Ngnp);
rr = RVE.rr(1:Ncnt + Ngnp,1:Ncnt + Ngnp); %tunelling distance
ja = dCNT;
Resistor_GNP = Resistor_GNP';

%% compute distances betwen GNPs
if Ngnp~=0
nber = RVE.nbr;
xg=Center_GNP(1,1);
yg=Center_GNP(1,3);
ae=rad_GNP(1,1)+RVE.dcut;
be=rad_GNP(2,1)+RVE.dcut;
alpha=rad_GNP(3,1);
%ellipse parameters with tunelling softcore
ABC(1,1) = ( ( (cos(alpha))^2 )/ ( ae^2 )  + ( (sin(alpha))^2 )/ ( be^2 ) );
ABC(2,1) = - 2 * cos (alpha) * sin(alpha) * (-1/(ae^2) + 1/(be^2));
ABC(3,1) = ( ((sin(alpha))^2)/ (ae^2)  + ((cos(alpha))^2)/ (be^2) );

for i=2:Ngnp
    p1 = [Resistor_GNP(2*i-1,1);Resistor_GNP(2*i-1,3)];%end point of ellipse i
    %p2 = [Resistor_GNP(2*i,1);Resistor_GNP(2*i,3)];%end point of ellipse i
    xg=Center_GNP(i,1);
    yg=Center_GNP(i,3);
    ae=rad_GNP(1,i);
    be=rad_GNP(2,i);
    alpha=rad_GNP(3,i);
    %ellipse parameters without tunneling softcore
    tt = linspace(0,2*pi,nber);
    X=ae*sin(tt)*cos(alpha)-be*cos(tt)*sin(alpha)+xg;
    Y=ae*sin(tt)*sin(alpha)+be*cos(tt)*cos(alpha)+yg;
    for j=1:i-1
        P2 = [Center_GNP(j,1);Center_GNP(j,3)];%center of j ellipse
        P2c = [Resistor_Center_GNP(j,1);Resistor_Center_GNP(j,3)];%fictitious center of j ellipse
        d_tun=rr(j + Ncnt,i + Ncnt);
        p3 = [Resistor_GNP(2*j-1,1);Resistor_GNP(2*j-1,3)];%end point of ellipse j
        p4 = [Resistor_GNP(2*j,1);Resistor_GNP(2*j,3)];%end point of ellipse j
        ae=rad_GNP(1,j);
        be=rad_GNP(2,j);
        alpha=rad_GNP(3,j);
        check=ABC(1,j).*(X-P2(1)).^2+ABC(2,j).*(X-P2(1)).*(Y-P2(2))+ABC(3,j).*(Y-P2(2)).^2-1;
        nb=find(check<=0);
        pnew = [X(nb);Y(nb)];%nodes on the i ellipse inside the j ellipse with tunelling
        if (~isempty(nb) & all(pnew(1,:)>=0) & all(pnew(1,:)<=RVE.size(1)) & all(pnew(2,:)>=0) & all(pnew(1,:)<=RVE.size(3))) %make sure that the points are inside the RVE
            [dt] = dist_multiple_points_ellipse( pnew, P2, ae, be, alpha,d_tun);
            P1c = [Resistor_Center_GNP(i,1);Resistor_Center_GNP(i,3)];%fictitious center of i ellipse
            dti = sqrt((P1c(1) - p1(1))^2 + (P1c(2) - p1(2))^2) + 1e-5;%I have added 0.1 so that the value won't be zer0 when the intersection is p1 or p3. I will subtract it
            dta = sqrt((P2c(1) - p3(1))^2 + (P2(2) - p3(2))^2) + 1e-5;%every time before using the value
            dCNT(Ncnt + j,Ncnt + i) = (dt<=RVE.dcut+1e-6) .* dt;
            ja(Ncnt + j,Ncnt + i) = (dt<=RVE.dcut+1e-6) .* dta;   %distance between junction on CNT a(j) and starting point of CNT a(j)
            ja(Ncnt + i,Ncnt + j) = (dt<=RVE.dcut+1e-6) .* dti;   %distance between junction on CNT i and starting point of CNT i
        end
    end
    %ellipse parameters with tunneling softcore
    ae=rad_GNP(1,i)+RVE.dcut;
    be=rad_GNP(2,i)+RVE.dcut;
    alpha=rad_GNP(3,i);
    ABC(1,i) = ( ( (cos(alpha))^2 )/ ( ae^2 )  + ( (sin(alpha))^2 )/ ( be^2 ) );
    ABC(2,i) = - 2 * cos (alpha) * sin(alpha) * (-1/(ae^2) + 1/(be^2));
    ABC(3,i) = ( ((sin(alpha))^2)/ (ae^2)  + ((cos(alpha))^2)/ (be^2) );
end
end
end

