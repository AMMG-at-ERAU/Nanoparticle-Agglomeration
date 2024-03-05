function [ Pgg,PPg,Center,Ae,rad_g ] = GNP_generate( RVE )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
ae = RVE.ae;
be = RVE.be;
Lx = RVE.size(1);
Ly = RVE.size(3); %it really is Lz

XYc=zeros(2,RVE.nbr);
ABC=[];
Pgg=[];
PPg=[];
rad_g=[];
Center=[];
Ae=[];
Vpg=0;
t=1;
tt = linspace(0,2*pi,RVE.nbr);
k=linspace(0,2*pi);
%%%% Random sequential Addition

while Vpg<RVE.perc*RVE.Vf*Lx*Ly
    
    if t>1
        check=-1; %to enter the while loop
        while any(check<0)
            reak=0;
            rg=rand(3,1);
            %%%center coordinates such as ellipse never crosses the boundaries
            xg=Lx*rg(1,:);
            yg=Ly*rg(2,:);
            %     mi=0;
            %     ma=pi/3;
            %     alpha=mi + (ma - mi)*rg(3,:); %%ellipse angle
            alpha=pi*rg(3,:); %%ellipse angle
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
        xg=Lx*rg(1,:);
        yg=Ly*rg(2,:);
        %     mi=0;
        %     ma=pi/3;
        %     alpha=mi + (ma - mi)*rg(3,:); %%ellipse angle
        alpha=pi*rg(3,:); %%ellipse angle
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
end
end

