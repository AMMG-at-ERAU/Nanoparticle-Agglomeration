function [ gh1,gh2,PPg ] = GNP_generate_strain_patterned( RVE )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
ae = RVE.ae;
be = RVE.be;
Lx = RVE.size(1);
Ly = RVE.size(3); %it really is Lz

XYc=zeros(2,RVE.nbr);
gh1=[];
gh2=[];
ABC=[];
Pgg=[];
PPg=[];
%Pcg=[];
rad_g=[];
Ae=[];
Center=[];
bb=[];
Vpg=0;
t=1;
tt = linspace(0,2*pi,RVE.nbr);
k=linspace(0,2*pi);
%%%% Random sequential Addition
%P = [Lx/6 2*Lx/6 3*Lx/6 4*Lx/6 5*Lx/6 Lx/5 2*Lx/5 3*Lx/5 4*Lx/5 Lx/6 2*Lx/6 3*Lx/6 4*Lx/6 5*Lx/6;Ly/4 Ly/4 Ly/4 Ly/4 Ly/4 2*Ly/4 2*Ly/4 2*Ly/4 2*Ly/4 3*Ly/4 3*Ly/4 3*Ly/4 3*Ly/4 3*Ly/4];
%%% Patern 44 ellipses a = 3, b = 1
%P1a = (Lx/22)+(Lx/11)*[0 1 2 3 4 5 6 7 8 9 10];
%P2a = (Ly/8)+(Ly/4)*[0 1 2 3];
%P = [P1a,P1a,P1a,P1a;ones(1,11)*P2a(1),ones(1,11)*P2a(2),ones(1,11)*P2a(3),ones(1,11)*P2a(4)]; 

% %% Patern 6 by 22 ellipses a = 3, b = 1, Lx = 32, Ly = 160
P1a = (Lx/22)+(Lx/11)*(0:1:10);
P2a = (Ly/24)+(Ly/12)*(0:1:11);
P = [P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a;
    ones(1,11)*P2a(1),ones(1,11)*P2a(2),ones(1,11)*P2a(3),ones(1,11)*P2a(4),ones(1,11)*P2a(5),ones(1,11)*P2a(6),ones(1,11)*P2a(7),ones(1,11)*P2a(8),ones(1,11)*P2a(9),ones(1,11)*P2a(10),ones(1,11)*P2a(11),ones(1,11)*P2a(12)];


%% 49 circles r = sart(3)
%P1a = (Lx/14)+(Lx/7)*[0 1 2 3 4 5 6];
%P2a = (Ly/14)+(Lx/7)*[0 1 2 3 4 5 6];
%P = [P1a,P1a,P1a,P1a,P1a,P1a,P1a;ones(1,7)*P2a(1),ones(1,7)*P2a(2),ones(1,7)*P2a(3),ones(1,7)*P2a(4),ones(1,7)*P2a(5),ones(1,7)*P2a(6),ones(1,7)*P2a(7)]; 

% %% Patern 6 by 22 circles a = sqrt(3), b = sqrt(3), Lx = 32, Ly = 160
%P1a = (Lx/12)+(Lx/6)*(0:1:5);
%P2a = (Ly/44)+(Ly/22)*(0:1:21);
%P = [P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a,P1a;
%    ones(1,6)*P2a(1),ones(1,6)*P2a(2),ones(1,6)*P2a(3),ones(1,6)*P2a(4),ones(1,6)*P2a(5),ones(1,6)*P2a(6),ones(1,6)*P2a(7),ones(1,6)*P2a(8),ones(1,6)*P2a(9),ones(1,6)*P2a(10),ones(1,6)*P2a(11),ones(1,6)*P2a(12),ones(1,6)*P2a(13),ones(1,6)*P2a(14),ones(1,6)*P2a(15),ones(1,6)*P2a(16),ones(1,6)*P2a(17),ones(1,6)*P2a(18),ones(1,6)*P2a(19),ones(1,6)*P2a(20),ones(1,6)*P2a(21),ones(1,6)*P2a(22)]; 


    %plot(P(1,:),P(2,:),'*')
iii = 1;
while Vpg<RVE.perc*RVE.Vf*Lx*Ly
    
    if t>1
        check=-1; %to enter the while loop
        while any(check<0)
            reak=0;
            rg=rand(3,1);
            %%%center coordinates such as ellipse never crosses the boundaries
        xg=P(1,iii);%Lx*rg(1,:);
        yg=P(2,iii);%Ly*rg(2,:);
            %     mi=0;
            %     ma=pi/6;
            %     alpha=mi + (ma - mi)*rg(3,:); %%ellipse angle
            % %     mi=pi/2;
            % %     ma=pi/2+6*pi/6;
            % %     alpha=mi + (ma-mi)*rg(3,:);
            alpha=pi/2;%*rg(3,:); %%ellipse angle
            Pg=[xg yg]; %ellipse center
            radius_g=[ae;be;alpha];
            
            if (xg>=ae & xg<=Lx-ae) & (yg>=ae & yg<=Ly-ae)
                P7=Pg';
                cent=P7;
                Pc=[xg-ae*cos(alpha) xg+ae*cos(alpha);yg-ae*sin(alpha) yg+ae*sin(alpha)];
                Sur=pi*ae*be;
                siz=1;
            else
                [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( Lx, Ly, Pg', ae, be, alpha );
                siz=size(P7,2);
            end
            
            for j=1:siz
                XYc(1,:)=ae*sin(tt)*cos(alpha)-be*cos(tt)*sin(alpha)+P7(1,j);
                XYc(2,:)=ae*sin(tt)*sin(alpha)+be*cos(tt)*cos(alpha)+P7(2,j);
                
                XYc1=[P7(:,j) XYc];%include the center of the ellipse so the check starts by the center
                for i=1:RVE.nbr %since 0 and 2*pi are the same take away 2*pi
                    %plug each of the nber points in all the ellipses previously
                    %computed equations to see if they are inside them
                    check=ABC(1,:).*(XYc1(1,i)-Pgg(1,:)).^2+ABC(2,:).*(XYc1(1,i)-Pgg(1,:)).*(XYc1(2,i)-Pgg(2,:))+ABC(3,:).*(XYc1(2,i)-Pgg(2,:)).^2-1;
                    if any(check < 0)
                        reak=1;
                        break
                    end
                end
                if reak==1
                    break
                end
                
            end
        end
        %ellipse parameters
        A = ( ( (cos(alpha))^2 )/ ( ae^2 )  + ( (sin(alpha))^2 )/ ( be^2 ) );
        B = - 2 * cos (alpha) * sin(alpha) * (-1/(ae^2) + 1/(be^2));
        C = ( ((sin(alpha))^2)/ (ae^2)  + ((cos(alpha))^2)/ (be^2) );
    else
        
        rg=rand(3,1);
        %%%center coordinates such as ellipse never crosses the boundaries
        xg=P(1,iii);%Lx*rg(1,:);
        yg=P(2,iii);%Ly*rg(2,:);
        %     mi=0;
        %     ma=pi/6;
        %     alpha=mi + (ma - mi)*rg(3,:); %%ellipse angle
        %     mi=pi/2;
        %     ma=pi/2+6*pi/6;
        %     alpha=mi + (ma-mi)*rg(3,:);
        alpha=pi/2;%*rg(3,:); %%ellipse angle
        Pg=[xg yg];
        radius_g=[ae;be;alpha];
        
        if (xg>=ae & xg<=Lx-ae) & (yg>=ae & yg<=Ly-ae)
            P7=Pg';
            cent=P7;
            Pc=[xg-ae*cos(alpha) xg+ae*cos(alpha);yg-ae*sin(alpha) yg+ae*sin(alpha)];
            Sur=pi*ae*be;
            siz=1;
        else
            [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( Lx, Ly, Pg', ae, be, alpha );
            siz=size(P7,2);
        end
        
        %ellipse parameters
        A = ( ( (cos(alpha))^2 )/ ( ae^2 )  + ( (sin(alpha))^2 )/ ( be^2 ) );
        B = - 2 * cos (alpha) * sin(alpha) * (-1/(ae^2) + 1/(be^2));
        C = ( ((sin(alpha))^2)/ (ae^2)  + ((cos(alpha))^2)/ (be^2) );
    end
    gh1=[gh1 [xg;yg]]; %center of first ellipses before BC are applied just for strain
    gh2=[gh2 alpha]; %angle of that first ellipse
    ABC=[ABC [A;B;C]*ones(1,siz)];
    Pgg=[Pgg P7];%ellipse real center coordinates
    PPg=[PPg Pc];%segments representing ellipse for Resistor Network
    Center=[Center cent];%center of segments for Resistor Network
    Ae=[Ae Sur];%area of ellipse
    rad_g=[rad_g radius_g*ones(1,siz)];
    Vp = pi*ae*be;
    Vpg=Vpg+Vp;
    t=t+siz;
    % for i=1:siz
    % xx=ae*sin(k)*cos(alpha)-be*cos(k)*sin(alpha)+P7(1,i);
    % yy=ae*sin(k)*sin(alpha)+be*cos(k)*cos(alpha)+P7(2,i);
    % plot(xx,yy,'r')                      %%ellipse
    % hold on
    % end
    iii = iii + 1;
end
end

