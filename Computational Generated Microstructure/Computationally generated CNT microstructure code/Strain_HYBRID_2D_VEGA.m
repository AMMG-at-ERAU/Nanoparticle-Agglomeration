myclusterSd=parcluster('local');
myclusterSd.NumWorkers=60;
parpool(myclusterSd,60,'SpmdEnabled',false);
Number = 5000;

parfor ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.093; %CNT area fraction
    RVE.Vf_agg = 0.10; %percentage of agglomerate CNTs


    RVE.perc = 0.2/RVE.Vf; %percentage of CNT fraction for GNP particles

    
    RVE.agg_angle = 60; %*/-angle around agglomerate center CNT in degree
    RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
    RVE.pois=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    


    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 4; %semi major axis of GNP
    RVE.be = 0.5; %semi minor axis of GNP
    RVE.nbr = 30*10+1;%Number of points on the ellipse contour
    RVE.D = 50e-3; %CNT diameter (micrometer)
    RVE.sigma = 1e4;  %CNT intrinsic conductivity (S/meter)
    RVE.Re = 2.8e2; %graphene sheet resistance ohms/m^2 (in reality just ohms)
    RVE.dvdw = 0.34e-3;   %van der waals (micrometer)
    RVE.dcut = 1.4e-3;                           %cut off tunnelling distance %micrometer
    RVE.me=9.1094*10^-31; %Kg
    RVE.lambda=1*1.602e-19;%Joules
    RVE.M=460;
    RVE.hp=6.62607004e-34; %meter^2 Kg/s
    RVE.e = 1.602e-19;                           %Coulomb
    
    %% Generate Softcore
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_M( RVE);% weibul length
    [ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_cst_agglo( RVE);%constant CNT length
    Ptrue = Ptrue';%[x y z;x y z]
    Center = Center';%[x y z;x y z]
    
    %%%Generate Hardcore GNP
    [ Pgg, gh, PPg ] = GNP_generate_strain( RVE );%Pgg (ellipse center no BC), gh (angle no BC), PPg (resistor with BC)
    Ngnp = length(gh);
    if Ngnp~=0
    Pgg = [Pgg(1,:);zeros(1,Ngnp);Pgg(2,:)]';
    end
    
    Linit = RVE.size;
    p = 1;
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(2*length(Length) + Ngnp,2*length(Length) + Ngnp);
    
    for strain = -0.006:0.006:0.006
        RVE.rr = rri*(1+strain);
        RVE.size = Linit .* [1-strain*RVE.pois,0,1+strain];
        RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
        Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
        
        %Reorientation of CNT after strain
        n=size(Center,1);
        deform = ones(n,1)*[1-strain*RVE.pois,0,1+strain];
        Center1 = Center .*deform;
        Ud = [(1-strain*RVE.pois)*cos(Angle'), 0*ones(n,1), (1+strain)*sin(Angle')];%[cos 0 sin]
        Ptrue1(1:2:2*n-1,:) = Center1 - 0.5*(Length'*ones(1,3)).*Ud;
        Ptrue1(2:2:2*n,:) = Center1 + 0.5*(Length'*ones(1,3)).*Ud;
        
        PiP=[];
        Pcenter=[];
        Plength=[];
        pt = [];
        for i = 1:n
            pt.P = Ptrue1(2*i-1:2*i,:);
            pt.center = Center1(i,:);
            pt.length = Length(i);
            [P, C, L] = PBC_2D(pt,RVE);
            PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
            Pcenter=[Pcenter C]; %%3 rows matrix with all centers
            Plength=[Plength L]; %%row with all lengths
        end
        
        Ptrue1 = PiP';
        Center1 = Pcenter';
        Length1 = Plength;
        
        %Reorientation of GNP after strain   
        Pggb=[];
        PPgb=[];
        rad_gb=[];
        Centerb=[];
        if Ngnp~=0
        deform = ones(Ngnp,1)*[1-strain*RVE.pois,0,1+strain];
        Pgg1 = Pgg .*deform;
        gh1 = atan(((1+strain)/(1-strain*RVE.pois))*tan(gh));
        for i=1:Ngnp
            xg=Pgg1(i,1);
            yg=Pgg1(i,3);
            Pg=[xg yg];
            alpha=gh1(i);
            radius_g=[RVE.ae;RVE.be;alpha];
            if (xg>=RVE.ae & xg<=RVE.size(1)-RVE.ae) & (yg>=RVE.ae & yg<=RVE.size(3)-RVE.ae)
                P7=Pg';
                cent=P7;
                Pc=[xg-RVE.ae*cos(alpha) xg+RVE.ae*cos(alpha);yg-RVE.ae*sin(alpha) yg+RVE.ae*sin(alpha)];
                Sur=pi*RVE.ae*RVE.be;
                siz=1;
            else
                [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( RVE.size(1), RVE.size(3), Pg', RVE.ae, RVE.be, alpha );
                siz=size(P7,2);
            end
            Pggb=[Pggb P7];%ellipse real center coordinates
            PPgb=[PPgb Pc];%segments representing ellipse for Resistor Network
            Centerb=[Centerb cent];%center of segments for Resistor Network
            rad_gb=[rad_gb radius_g*ones(1,siz)];
        end
        
        ng=size(PPgb,2);
        Pggb = [Pggb(1,:);zeros(1,0.5*ng);Pggb(2,:)]';
        PPgb = [PPgb(1,:);zeros(1,ng);PPgb(2,:)]';
        Centerb = [Centerb(1,:);zeros(1,0.5*ng);Centerb(2,:)]';
        
        end
        
        
        [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPgb, Centerb, Pggb, rad_gb, RVE);%Find all tunneling distances and the junctions on CNTs
        Ncnt = length(Length1);
        PP = [Ptrue1;PPgb];
        [Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
        

        if perco ==0
%             disp('// Warning: No percolation cluster exist!');
            result(1,1:3) = 0;
            fLAG(1,1:3) = NaN;
            junct(1,1:9) = NaN;
            break
        else 
            bn = zeros(length(Gn),1);
            bn(length(Gn))= 1;
            setup = [];
            setup.droptol = 1e-15;                       %tolerance
            maxit = 3e6;
            setup.type = 'ict';
            
            L1=ichol(Gn,setup);
            
            [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
            flag1
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            flag2
            Res2 = xn(length(Gn));
%             disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct(1,3*p-2:3*p) = Nj;
            p = p+1;
        end
        
    end 
    result = (result-result(1,2))./result(1,2);
    fLAG = fLAG*1;
    junct = junct*1;
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct(ii,:) = junct(1:end);
%     disp(['// Change % = ',num2str(100*result1(ii,:)),' %']);
%     disp(['// Flag = ',num2str(FLAG(ii,:)),' %']);
%     disp(['// Junct = ',num2str(Junct(ii,:)),' %']);
end
save('Hybrid_Strain_2D_Vf9.30p_G20p_agg_10p_angle_60d_rad_2L-over-100_5000_P100.mat','result1','Junct','FLAG')



parfor ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.0920; %CNT area fraction
    RVE.Vf_agg = 0.20; %percentage of agglomerate CNTs


    RVE.perc = 0.2/RVE.Vf; %percentage of CNT fraction for GNP particles

    
    RVE.agg_angle = 60; %*/-angle around agglomerate center CNT in degree
    RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
    RVE.pois=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    


    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 4; %semi major axis of GNP
    RVE.be = 0.5; %semi minor axis of GNP
    RVE.nbr = 30*10+1;%Number of points on the ellipse contour
    RVE.D = 50e-3; %CNT diameter (micrometer)
    RVE.sigma = 1e4;  %CNT intrinsic conductivity (S/meter)
    RVE.Re = 2.8e2; %graphene sheet resistance ohms/m^2 (in reality just ohms)
    RVE.dvdw = 0.34e-3;   %van der waals (micrometer)
    RVE.dcut = 1.4e-3;                           %cut off tunnelling distance %micrometer
    RVE.me=9.1094*10^-31; %Kg
    RVE.lambda=1*1.602e-19;%Joules
    RVE.M=460;
    RVE.hp=6.62607004e-34; %meter^2 Kg/s
    RVE.e = 1.602e-19;                           %Coulomb
    
    %% Generate Softcore
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_M( RVE);% weibul length
    [ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_cst_agglo( RVE);%constant CNT length
    Ptrue = Ptrue';%[x y z;x y z]
    Center = Center';%[x y z;x y z]
    
    %%%Generate Hardcore GNP
    [ Pgg, gh, PPg ] = GNP_generate_strain( RVE );%Pgg (ellipse center no BC), gh (angle no BC), PPg (resistor with BC)
    Ngnp = length(gh);
    if Ngnp~=0
    Pgg = [Pgg(1,:);zeros(1,Ngnp);Pgg(2,:)]';
    end
    
    Linit = RVE.size;
    p = 1;
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(2*length(Length) + Ngnp,2*length(Length) + Ngnp);
    
    for strain = -0.006:0.006:0.006
        RVE.rr = rri*(1+strain);
        RVE.size = Linit .* [1-strain*RVE.pois,0,1+strain];
        RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
        Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
        
        %Reorientation of CNT after strain
        n=size(Center,1);
        deform = ones(n,1)*[1-strain*RVE.pois,0,1+strain];
        Center1 = Center .*deform;
        Ud = [(1-strain*RVE.pois)*cos(Angle'), 0*ones(n,1), (1+strain)*sin(Angle')];%[cos 0 sin]
        Ptrue1(1:2:2*n-1,:) = Center1 - 0.5*(Length'*ones(1,3)).*Ud;
        Ptrue1(2:2:2*n,:) = Center1 + 0.5*(Length'*ones(1,3)).*Ud;
        
        PiP=[];
        Pcenter=[];
        Plength=[];
        pt = [];
        for i = 1:n
            pt.P = Ptrue1(2*i-1:2*i,:);
            pt.center = Center1(i,:);
            pt.length = Length(i);
            [P, C, L] = PBC_2D(pt,RVE);
            PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
            Pcenter=[Pcenter C]; %%3 rows matrix with all centers
            Plength=[Plength L]; %%row with all lengths
        end
        
        Ptrue1 = PiP';
        Center1 = Pcenter';
        Length1 = Plength;
        
        %Reorientation of GNP after strain   
        Pggb=[];
        PPgb=[];
        rad_gb=[];
        Centerb=[];
        if Ngnp~=0
        deform = ones(Ngnp,1)*[1-strain*RVE.pois,0,1+strain];
        Pgg1 = Pgg .*deform;
        gh1 = atan(((1+strain)/(1-strain*RVE.pois))*tan(gh));
        for i=1:Ngnp
            xg=Pgg1(i,1);
            yg=Pgg1(i,3);
            Pg=[xg yg];
            alpha=gh1(i);
            radius_g=[RVE.ae;RVE.be;alpha];
            if (xg>=RVE.ae & xg<=RVE.size(1)-RVE.ae) & (yg>=RVE.ae & yg<=RVE.size(3)-RVE.ae)
                P7=Pg';
                cent=P7;
                Pc=[xg-RVE.ae*cos(alpha) xg+RVE.ae*cos(alpha);yg-RVE.ae*sin(alpha) yg+RVE.ae*sin(alpha)];
                Sur=pi*RVE.ae*RVE.be;
                siz=1;
            else
                [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( RVE.size(1), RVE.size(3), Pg', RVE.ae, RVE.be, alpha );
                siz=size(P7,2);
            end
            Pggb=[Pggb P7];%ellipse real center coordinates
            PPgb=[PPgb Pc];%segments representing ellipse for Resistor Network
            Centerb=[Centerb cent];%center of segments for Resistor Network
            rad_gb=[rad_gb radius_g*ones(1,siz)];
        end
        
        ng=size(PPgb,2);
        Pggb = [Pggb(1,:);zeros(1,0.5*ng);Pggb(2,:)]';
        PPgb = [PPgb(1,:);zeros(1,ng);PPgb(2,:)]';
        Centerb = [Centerb(1,:);zeros(1,0.5*ng);Centerb(2,:)]';
        
        end
        
        
        [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPgb, Centerb, Pggb, rad_gb, RVE);%Find all tunneling distances and the junctions on CNTs
        Ncnt = length(Length1);
        PP = [Ptrue1;PPgb];
        [Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
        

        if perco ==0
%             disp('// Warning: No percolation cluster exist!');
            result(1,1:3) = 0;
            fLAG(1,1:3) = NaN;
            junct(1,1:9) = NaN;
            break
        else 
            bn = zeros(length(Gn),1);
            bn(length(Gn))= 1;
            setup = [];
            setup.droptol = 1e-15;                       %tolerance
            maxit = 3e6;
            setup.type = 'ict';
            
            L1=ichol(Gn,setup);
            
            [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
            flag1
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            flag2
            Res2 = xn(length(Gn));
%             disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct(1,3*p-2:3*p) = Nj;
            p = p+1;
        end
        
    end 
    result = (result-result(1,2))./result(1,2);
    fLAG = fLAG*1;
    junct = junct*1;
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct(ii,:) = junct(1:end);
%     disp(['// Change % = ',num2str(100*result1(ii,:)),' %']);
%     disp(['// Flag = ',num2str(FLAG(ii,:)),' %']);
%     disp(['// Junct = ',num2str(Junct(ii,:)),' %']);
end
save('Hybrid_Strain_2D_Vf9.20p_G20p_agg_20p_angle_60d_rad_2L-over-100_5000_P100.mat','result1','Junct','FLAG')



parfor ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.0955; %CNT area fraction
    RVE.Vf_agg = 0.30; %percentage of agglomerate CNTs


    RVE.perc = 0.2/RVE.Vf; %percentage of CNT fraction for GNP particles

    
    RVE.agg_angle = 60; %*/-angle around agglomerate center CNT in degree
    RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
    RVE.pois=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    


    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 4; %semi major axis of GNP
    RVE.be = 0.5; %semi minor axis of GNP
    RVE.nbr = 30*10+1;%Number of points on the ellipse contour
    RVE.D = 50e-3; %CNT diameter (micrometer)
    RVE.sigma = 1e4;  %CNT intrinsic conductivity (S/meter)
    RVE.Re = 2.8e2; %graphene sheet resistance ohms/m^2 (in reality just ohms)
    RVE.dvdw = 0.34e-3;   %van der waals (micrometer)
    RVE.dcut = 1.4e-3;                           %cut off tunnelling distance %micrometer
    RVE.me=9.1094*10^-31; %Kg
    RVE.lambda=1*1.602e-19;%Joules
    RVE.M=460;
    RVE.hp=6.62607004e-34; %meter^2 Kg/s
    RVE.e = 1.602e-19;                           %Coulomb
    
    %% Generate Softcore
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_M( RVE);% weibul length
    [ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_cst_agglo( RVE);%constant CNT length
    Ptrue = Ptrue';%[x y z;x y z]
    Center = Center';%[x y z;x y z]
    
    %%%Generate Hardcore GNP
    [ Pgg, gh, PPg ] = GNP_generate_strain( RVE );%Pgg (ellipse center no BC), gh (angle no BC), PPg (resistor with BC)
    Ngnp = length(gh);
    if Ngnp~=0
    Pgg = [Pgg(1,:);zeros(1,Ngnp);Pgg(2,:)]';
    end
    
    Linit = RVE.size;
    p = 1;
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(2*length(Length) + Ngnp,2*length(Length) + Ngnp);
    
    for strain = -0.006:0.006:0.006
        RVE.rr = rri*(1+strain);
        RVE.size = Linit .* [1-strain*RVE.pois,0,1+strain];
        RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
        Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
        
        %Reorientation of CNT after strain
        n=size(Center,1);
        deform = ones(n,1)*[1-strain*RVE.pois,0,1+strain];
        Center1 = Center .*deform;
        Ud = [(1-strain*RVE.pois)*cos(Angle'), 0*ones(n,1), (1+strain)*sin(Angle')];%[cos 0 sin]
        Ptrue1(1:2:2*n-1,:) = Center1 - 0.5*(Length'*ones(1,3)).*Ud;
        Ptrue1(2:2:2*n,:) = Center1 + 0.5*(Length'*ones(1,3)).*Ud;
        
        PiP=[];
        Pcenter=[];
        Plength=[];
        pt = [];
        for i = 1:n
            pt.P = Ptrue1(2*i-1:2*i,:);
            pt.center = Center1(i,:);
            pt.length = Length(i);
            [P, C, L] = PBC_2D(pt,RVE);
            PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
            Pcenter=[Pcenter C]; %%3 rows matrix with all centers
            Plength=[Plength L]; %%row with all lengths
        end
        
        Ptrue1 = PiP';
        Center1 = Pcenter';
        Length1 = Plength;
        
        %Reorientation of GNP after strain   
        Pggb=[];
        PPgb=[];
        rad_gb=[];
        Centerb=[];
        if Ngnp~=0
        deform = ones(Ngnp,1)*[1-strain*RVE.pois,0,1+strain];
        Pgg1 = Pgg .*deform;
        gh1 = atan(((1+strain)/(1-strain*RVE.pois))*tan(gh));
        for i=1:Ngnp
            xg=Pgg1(i,1);
            yg=Pgg1(i,3);
            Pg=[xg yg];
            alpha=gh1(i);
            radius_g=[RVE.ae;RVE.be;alpha];
            if (xg>=RVE.ae & xg<=RVE.size(1)-RVE.ae) & (yg>=RVE.ae & yg<=RVE.size(3)-RVE.ae)
                P7=Pg';
                cent=P7;
                Pc=[xg-RVE.ae*cos(alpha) xg+RVE.ae*cos(alpha);yg-RVE.ae*sin(alpha) yg+RVE.ae*sin(alpha)];
                Sur=pi*RVE.ae*RVE.be;
                siz=1;
            else
                [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( RVE.size(1), RVE.size(3), Pg', RVE.ae, RVE.be, alpha );
                siz=size(P7,2);
            end
            Pggb=[Pggb P7];%ellipse real center coordinates
            PPgb=[PPgb Pc];%segments representing ellipse for Resistor Network
            Centerb=[Centerb cent];%center of segments for Resistor Network
            rad_gb=[rad_gb radius_g*ones(1,siz)];
        end
        
        ng=size(PPgb,2);
        Pggb = [Pggb(1,:);zeros(1,0.5*ng);Pggb(2,:)]';
        PPgb = [PPgb(1,:);zeros(1,ng);PPgb(2,:)]';
        Centerb = [Centerb(1,:);zeros(1,0.5*ng);Centerb(2,:)]';
        
        end
        
        
        [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPgb, Centerb, Pggb, rad_gb, RVE);%Find all tunneling distances and the junctions on CNTs
        Ncnt = length(Length1);
        PP = [Ptrue1;PPgb];
        [Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
        

        if perco ==0
%             disp('// Warning: No percolation cluster exist!');
            result(1,1:3) = 0;
            fLAG(1,1:3) = NaN;
            junct(1,1:9) = NaN;
            break
        else 
            bn = zeros(length(Gn),1);
            bn(length(Gn))= 1;
            setup = [];
            setup.droptol = 1e-15;                       %tolerance
            maxit = 3e6;
            setup.type = 'ict';
            
            L1=ichol(Gn,setup);
            
            [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
            flag1
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            flag2
            Res2 = xn(length(Gn));
%             disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct(1,3*p-2:3*p) = Nj;
            p = p+1;
        end
        
    end 
    result = (result-result(1,2))./result(1,2);
    fLAG = fLAG*1;
    junct = junct*1;
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct(ii,:) = junct(1:end);
%     disp(['// Change % = ',num2str(100*result1(ii,:)),' %']);
%     disp(['// Flag = ',num2str(FLAG(ii,:)),' %']);
%     disp(['// Junct = ',num2str(Junct(ii,:)),' %']);
end
save('Hybrid_Strain_2D_Vf9.55p_G20p_agg_30p_angle_60d_rad_2L-over-100_5000_P100.mat','result1','Junct','FLAG')


parfor ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.10; %CNT area fraction
    RVE.Vf_agg = 0.40; %percentage of agglomerate CNTs


    RVE.perc = 0.2/RVE.Vf; %percentage of CNT fraction for GNP particles

    
    RVE.agg_angle = 60; %*/-angle around agglomerate center CNT in degree
    RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
    RVE.pois=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    


    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 4; %semi major axis of GNP
    RVE.be = 0.5; %semi minor axis of GNP
    RVE.nbr = 30*10+1;%Number of points on the ellipse contour
    RVE.D = 50e-3; %CNT diameter (micrometer)
    RVE.sigma = 1e4;  %CNT intrinsic conductivity (S/meter)
    RVE.Re = 2.8e2; %graphene sheet resistance ohms/m^2 (in reality just ohms)
    RVE.dvdw = 0.34e-3;   %van der waals (micrometer)
    RVE.dcut = 1.4e-3;                           %cut off tunnelling distance %micrometer
    RVE.me=9.1094*10^-31; %Kg
    RVE.lambda=1*1.602e-19;%Joules
    RVE.M=460;
    RVE.hp=6.62607004e-34; %meter^2 Kg/s
    RVE.e = 1.602e-19;                           %Coulomb
    
    %% Generate Softcore
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_M( RVE);% weibul length
    [ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_cst_agglo( RVE);%constant CNT length
    Ptrue = Ptrue';%[x y z;x y z]
    Center = Center';%[x y z;x y z]
    
    %%%Generate Hardcore GNP
    [ Pgg, gh, PPg ] = GNP_generate_strain( RVE );%Pgg (ellipse center no BC), gh (angle no BC), PPg (resistor with BC)
    Ngnp = length(gh);
    if Ngnp~=0
    Pgg = [Pgg(1,:);zeros(1,Ngnp);Pgg(2,:)]';
    end
    
    Linit = RVE.size;
    p = 1;
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(2*length(Length) + Ngnp,2*length(Length) + Ngnp);
    
    for strain = -0.006:0.006:0.006
        RVE.rr = rri*(1+strain);
        RVE.size = Linit .* [1-strain*RVE.pois,0,1+strain];
        RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
        Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
        
        %Reorientation of CNT after strain
        n=size(Center,1);
        deform = ones(n,1)*[1-strain*RVE.pois,0,1+strain];
        Center1 = Center .*deform;
        Ud = [(1-strain*RVE.pois)*cos(Angle'), 0*ones(n,1), (1+strain)*sin(Angle')];%[cos 0 sin]
        Ptrue1(1:2:2*n-1,:) = Center1 - 0.5*(Length'*ones(1,3)).*Ud;
        Ptrue1(2:2:2*n,:) = Center1 + 0.5*(Length'*ones(1,3)).*Ud;
        
        PiP=[];
        Pcenter=[];
        Plength=[];
        pt = [];
        for i = 1:n
            pt.P = Ptrue1(2*i-1:2*i,:);
            pt.center = Center1(i,:);
            pt.length = Length(i);
            [P, C, L] = PBC_2D(pt,RVE);
            PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
            Pcenter=[Pcenter C]; %%3 rows matrix with all centers
            Plength=[Plength L]; %%row with all lengths
        end
        
        Ptrue1 = PiP';
        Center1 = Pcenter';
        Length1 = Plength;
        
        %Reorientation of GNP after strain   
        Pggb=[];
        PPgb=[];
        rad_gb=[];
        Centerb=[];
        if Ngnp~=0
        deform = ones(Ngnp,1)*[1-strain*RVE.pois,0,1+strain];
        Pgg1 = Pgg .*deform;
        gh1 = atan(((1+strain)/(1-strain*RVE.pois))*tan(gh));
        for i=1:Ngnp
            xg=Pgg1(i,1);
            yg=Pgg1(i,3);
            Pg=[xg yg];
            alpha=gh1(i);
            radius_g=[RVE.ae;RVE.be;alpha];
            if (xg>=RVE.ae & xg<=RVE.size(1)-RVE.ae) & (yg>=RVE.ae & yg<=RVE.size(3)-RVE.ae)
                P7=Pg';
                cent=P7;
                Pc=[xg-RVE.ae*cos(alpha) xg+RVE.ae*cos(alpha);yg-RVE.ae*sin(alpha) yg+RVE.ae*sin(alpha)];
                Sur=pi*RVE.ae*RVE.be;
                siz=1;
            else
                [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( RVE.size(1), RVE.size(3), Pg', RVE.ae, RVE.be, alpha );
                siz=size(P7,2);
            end
            Pggb=[Pggb P7];%ellipse real center coordinates
            PPgb=[PPgb Pc];%segments representing ellipse for Resistor Network
            Centerb=[Centerb cent];%center of segments for Resistor Network
            rad_gb=[rad_gb radius_g*ones(1,siz)];
        end
        
        ng=size(PPgb,2);
        Pggb = [Pggb(1,:);zeros(1,0.5*ng);Pggb(2,:)]';
        PPgb = [PPgb(1,:);zeros(1,ng);PPgb(2,:)]';
        Centerb = [Centerb(1,:);zeros(1,0.5*ng);Centerb(2,:)]';
        
        end
        
        
        [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPgb, Centerb, Pggb, rad_gb, RVE);%Find all tunneling distances and the junctions on CNTs
        Ncnt = length(Length1);
        PP = [Ptrue1;PPgb];
        [Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
        

        if perco ==0
%             disp('// Warning: No percolation cluster exist!');
            result(1,1:3) = 0;
            fLAG(1,1:3) = NaN;
            junct(1,1:9) = NaN;
            break
        else 
            bn = zeros(length(Gn),1);
            bn(length(Gn))= 1;
            setup = [];
            setup.droptol = 1e-15;                       %tolerance
            maxit = 3e6;
            setup.type = 'ict';
            
            L1=ichol(Gn,setup);
            
            [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
            flag1
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            flag2
            Res2 = xn(length(Gn));
%             disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct(1,3*p-2:3*p) = Nj;
            p = p+1;
        end
        
    end 
    result = (result-result(1,2))./result(1,2);
    fLAG = fLAG*1;
    junct = junct*1;
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct(ii,:) = junct(1:end);
%     disp(['// Change % = ',num2str(100*result1(ii,:)),' %']);
%     disp(['// Flag = ',num2str(FLAG(ii,:)),' %']);
%     disp(['// Junct = ',num2str(Junct(ii,:)),' %']);
end
save('Hybrid_Strain_2D_Vf10p_G20p_agg_40p_angle_60d_rad_2L-over-100_5000_P100.mat','result1','Junct','FLAG')



parfor ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.1030; %CNT area fraction
    RVE.Vf_agg = 0.50; %percentage of agglomerate CNTs


    RVE.perc = 0.2/RVE.Vf; %percentage of CNT fraction for GNP particles

    
    RVE.agg_angle = 60; %*/-angle around agglomerate center CNT in degree
    RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
    RVE.pois=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    


    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 4; %semi major axis of GNP
    RVE.be = 0.5; %semi minor axis of GNP
    RVE.nbr = 30*10+1;%Number of points on the ellipse contour
    RVE.D = 50e-3; %CNT diameter (micrometer)
    RVE.sigma = 1e4;  %CNT intrinsic conductivity (S/meter)
    RVE.Re = 2.8e2; %graphene sheet resistance ohms/m^2 (in reality just ohms)
    RVE.dvdw = 0.34e-3;   %van der waals (micrometer)
    RVE.dcut = 1.4e-3;                           %cut off tunnelling distance %micrometer
    RVE.me=9.1094*10^-31; %Kg
    RVE.lambda=1*1.602e-19;%Joules
    RVE.M=460;
    RVE.hp=6.62607004e-34; %meter^2 Kg/s
    RVE.e = 1.602e-19;                           %Coulomb
    
    %% Generate Softcore
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_M( RVE);% weibul length
    [ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_cst_agglo( RVE);%constant CNT length
    Ptrue = Ptrue';%[x y z;x y z]
    Center = Center';%[x y z;x y z]
    
    %%%Generate Hardcore GNP
    [ Pgg, gh, PPg ] = GNP_generate_strain( RVE );%Pgg (ellipse center no BC), gh (angle no BC), PPg (resistor with BC)
    Ngnp = length(gh);
    if Ngnp~=0
    Pgg = [Pgg(1,:);zeros(1,Ngnp);Pgg(2,:)]';
    end
    
    Linit = RVE.size;
    p = 1;
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(2*length(Length) + Ngnp,2*length(Length) + Ngnp);
    
    for strain = -0.006:0.006:0.006
        RVE.rr = rri*(1+strain);
        RVE.size = Linit .* [1-strain*RVE.pois,0,1+strain];
        RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
        Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
        
        %Reorientation of CNT after strain
        n=size(Center,1);
        deform = ones(n,1)*[1-strain*RVE.pois,0,1+strain];
        Center1 = Center .*deform;
        Ud = [(1-strain*RVE.pois)*cos(Angle'), 0*ones(n,1), (1+strain)*sin(Angle')];%[cos 0 sin]
        Ptrue1(1:2:2*n-1,:) = Center1 - 0.5*(Length'*ones(1,3)).*Ud;
        Ptrue1(2:2:2*n,:) = Center1 + 0.5*(Length'*ones(1,3)).*Ud;
        
        PiP=[];
        Pcenter=[];
        Plength=[];
        pt = [];
        for i = 1:n
            pt.P = Ptrue1(2*i-1:2*i,:);
            pt.center = Center1(i,:);
            pt.length = Length(i);
            [P, C, L] = PBC_2D(pt,RVE);
            PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
            Pcenter=[Pcenter C]; %%3 rows matrix with all centers
            Plength=[Plength L]; %%row with all lengths
        end
        
        Ptrue1 = PiP';
        Center1 = Pcenter';
        Length1 = Plength;
        
        %Reorientation of GNP after strain   
        Pggb=[];
        PPgb=[];
        rad_gb=[];
        Centerb=[];
        if Ngnp~=0
        deform = ones(Ngnp,1)*[1-strain*RVE.pois,0,1+strain];
        Pgg1 = Pgg .*deform;
        gh1 = atan(((1+strain)/(1-strain*RVE.pois))*tan(gh));
        for i=1:Ngnp
            xg=Pgg1(i,1);
            yg=Pgg1(i,3);
            Pg=[xg yg];
            alpha=gh1(i);
            radius_g=[RVE.ae;RVE.be;alpha];
            if (xg>=RVE.ae & xg<=RVE.size(1)-RVE.ae) & (yg>=RVE.ae & yg<=RVE.size(3)-RVE.ae)
                P7=Pg';
                cent=P7;
                Pc=[xg-RVE.ae*cos(alpha) xg+RVE.ae*cos(alpha);yg-RVE.ae*sin(alpha) yg+RVE.ae*sin(alpha)];
                Sur=pi*RVE.ae*RVE.be;
                siz=1;
            else
                [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( RVE.size(1), RVE.size(3), Pg', RVE.ae, RVE.be, alpha );
                siz=size(P7,2);
            end
            Pggb=[Pggb P7];%ellipse real center coordinates
            PPgb=[PPgb Pc];%segments representing ellipse for Resistor Network
            Centerb=[Centerb cent];%center of segments for Resistor Network
            rad_gb=[rad_gb radius_g*ones(1,siz)];
        end
        
        ng=size(PPgb,2);
        Pggb = [Pggb(1,:);zeros(1,0.5*ng);Pggb(2,:)]';
        PPgb = [PPgb(1,:);zeros(1,ng);PPgb(2,:)]';
        Centerb = [Centerb(1,:);zeros(1,0.5*ng);Centerb(2,:)]';
        
        end
        
        
        [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPgb, Centerb, Pggb, rad_gb, RVE);%Find all tunneling distances and the junctions on CNTs
        Ncnt = length(Length1);
        PP = [Ptrue1;PPgb];
        [Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
        

        if perco ==0
%             disp('// Warning: No percolation cluster exist!');
            result(1,1:3) = 0;
            fLAG(1,1:3) = NaN;
            junct(1,1:9) = NaN;
            break
        else 
            bn = zeros(length(Gn),1);
            bn(length(Gn))= 1;
            setup = [];
            setup.droptol = 1e-15;                       %tolerance
            maxit = 3e6;
            setup.type = 'ict';
            
            L1=ichol(Gn,setup);
            
            [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
            flag1
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            flag2
            Res2 = xn(length(Gn));
%             disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct(1,3*p-2:3*p) = Nj;
            p = p+1;
        end
        
    end 
    result = (result-result(1,2))./result(1,2);
    fLAG = fLAG*1;
    junct = junct*1;
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct(ii,:) = junct(1:end);
%     disp(['// Change % = ',num2str(100*result1(ii,:)),' %']);
%     disp(['// Flag = ',num2str(FLAG(ii,:)),' %']);
%     disp(['// Junct = ',num2str(Junct(ii,:)),' %']);
end
save('Hybrid_Strain_2D_Vf10.30p_G20p_agg_50p_angle_60d_rad_2L-over-100_5000_P100.mat','result1','Junct','FLAG')


parfor ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.1110; %CNT area fraction
    RVE.Vf_agg = 0.60; %percentage of agglomerate CNTs


    RVE.perc = 0.2/RVE.Vf; %percentage of CNT fraction for GNP particles

    
    RVE.agg_angle = 60; %*/-angle around agglomerate center CNT in degree
    RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
    RVE.pois=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    


    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 4; %semi major axis of GNP
    RVE.be = 0.5; %semi minor axis of GNP
    RVE.nbr = 30*10+1;%Number of points on the ellipse contour
    RVE.D = 50e-3; %CNT diameter (micrometer)
    RVE.sigma = 1e4;  %CNT intrinsic conductivity (S/meter)
    RVE.Re = 2.8e2; %graphene sheet resistance ohms/m^2 (in reality just ohms)
    RVE.dvdw = 0.34e-3;   %van der waals (micrometer)
    RVE.dcut = 1.4e-3;                           %cut off tunnelling distance %micrometer
    RVE.me=9.1094*10^-31; %Kg
    RVE.lambda=1*1.602e-19;%Joules
    RVE.M=460;
    RVE.hp=6.62607004e-34; %meter^2 Kg/s
    RVE.e = 1.602e-19;                           %Coulomb
    
    %% Generate Softcore
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_M( RVE);% weibul length
    [ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_cst_agglo( RVE);%constant CNT length
    Ptrue = Ptrue';%[x y z;x y z]
    Center = Center';%[x y z;x y z]
    
    %%%Generate Hardcore GNP
    [ Pgg, gh, PPg ] = GNP_generate_strain( RVE );%Pgg (ellipse center no BC), gh (angle no BC), PPg (resistor with BC)
    Ngnp = length(gh);
    if Ngnp~=0
    Pgg = [Pgg(1,:);zeros(1,Ngnp);Pgg(2,:)]';
    end
    
    Linit = RVE.size;
    p = 1;
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(2*length(Length) + Ngnp,2*length(Length) + Ngnp);
    
    for strain = -0.006:0.006:0.006
        RVE.rr = rri*(1+strain);
        RVE.size = Linit .* [1-strain*RVE.pois,0,1+strain];
        RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
        Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
        
        %Reorientation of CNT after strain
        n=size(Center,1);
        deform = ones(n,1)*[1-strain*RVE.pois,0,1+strain];
        Center1 = Center .*deform;
        Ud = [(1-strain*RVE.pois)*cos(Angle'), 0*ones(n,1), (1+strain)*sin(Angle')];%[cos 0 sin]
        Ptrue1(1:2:2*n-1,:) = Center1 - 0.5*(Length'*ones(1,3)).*Ud;
        Ptrue1(2:2:2*n,:) = Center1 + 0.5*(Length'*ones(1,3)).*Ud;
        
        PiP=[];
        Pcenter=[];
        Plength=[];
        pt = [];
        for i = 1:n
            pt.P = Ptrue1(2*i-1:2*i,:);
            pt.center = Center1(i,:);
            pt.length = Length(i);
            [P, C, L] = PBC_2D(pt,RVE);
            PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
            Pcenter=[Pcenter C]; %%3 rows matrix with all centers
            Plength=[Plength L]; %%row with all lengths
        end
        
        Ptrue1 = PiP';
        Center1 = Pcenter';
        Length1 = Plength;
        
        %Reorientation of GNP after strain   
        Pggb=[];
        PPgb=[];
        rad_gb=[];
        Centerb=[];
        if Ngnp~=0
        deform = ones(Ngnp,1)*[1-strain*RVE.pois,0,1+strain];
        Pgg1 = Pgg .*deform;
        gh1 = atan(((1+strain)/(1-strain*RVE.pois))*tan(gh));
        for i=1:Ngnp
            xg=Pgg1(i,1);
            yg=Pgg1(i,3);
            Pg=[xg yg];
            alpha=gh1(i);
            radius_g=[RVE.ae;RVE.be;alpha];
            if (xg>=RVE.ae & xg<=RVE.size(1)-RVE.ae) & (yg>=RVE.ae & yg<=RVE.size(3)-RVE.ae)
                P7=Pg';
                cent=P7;
                Pc=[xg-RVE.ae*cos(alpha) xg+RVE.ae*cos(alpha);yg-RVE.ae*sin(alpha) yg+RVE.ae*sin(alpha)];
                Sur=pi*RVE.ae*RVE.be;
                siz=1;
            else
                [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( RVE.size(1), RVE.size(3), Pg', RVE.ae, RVE.be, alpha );
                siz=size(P7,2);
            end
            Pggb=[Pggb P7];%ellipse real center coordinates
            PPgb=[PPgb Pc];%segments representing ellipse for Resistor Network
            Centerb=[Centerb cent];%center of segments for Resistor Network
            rad_gb=[rad_gb radius_g*ones(1,siz)];
        end
        
        ng=size(PPgb,2);
        Pggb = [Pggb(1,:);zeros(1,0.5*ng);Pggb(2,:)]';
        PPgb = [PPgb(1,:);zeros(1,ng);PPgb(2,:)]';
        Centerb = [Centerb(1,:);zeros(1,0.5*ng);Centerb(2,:)]';
        
        end
        
        
        [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPgb, Centerb, Pggb, rad_gb, RVE);%Find all tunneling distances and the junctions on CNTs
        Ncnt = length(Length1);
        PP = [Ptrue1;PPgb];
        [Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
        

        if perco ==0
%             disp('// Warning: No percolation cluster exist!');
            result(1,1:3) = 0;
            fLAG(1,1:3) = NaN;
            junct(1,1:9) = NaN;
            break
        else 
            bn = zeros(length(Gn),1);
            bn(length(Gn))= 1;
            setup = [];
            setup.droptol = 1e-15;                       %tolerance
            maxit = 3e6;
            setup.type = 'ict';
            
            L1=ichol(Gn,setup);
            
            [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
            flag1
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            flag2
            Res2 = xn(length(Gn));
%             disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct(1,3*p-2:3*p) = Nj;
            p = p+1;
        end
        
    end 
    result = (result-result(1,2))./result(1,2);
    fLAG = fLAG*1;
    junct = junct*1;
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct(ii,:) = junct(1:end);
%     disp(['// Change % = ',num2str(100*result1(ii,:)),' %']);
%     disp(['// Flag = ',num2str(FLAG(ii,:)),' %']);
%     disp(['// Junct = ',num2str(Junct(ii,:)),' %']);
end
save('Hybrid_Strain_2D_Vf11.10p_G20p_agg_60p_angle_60d_rad_2L-over-100_5000_P100.mat','result1','Junct','FLAG')


parfor ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.1180; %CNT area fraction
    RVE.Vf_agg = 0.70; %percentage of agglomerate CNTs


    RVE.perc = 0.2/RVE.Vf; %percentage of CNT fraction for GNP particles

    
    RVE.agg_angle = 60; %*/-angle around agglomerate center CNT in degree
    RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
    RVE.pois=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    


    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 4; %semi major axis of GNP
    RVE.be = 0.5; %semi minor axis of GNP
    RVE.nbr = 30*10+1;%Number of points on the ellipse contour
    RVE.D = 50e-3; %CNT diameter (micrometer)
    RVE.sigma = 1e4;  %CNT intrinsic conductivity (S/meter)
    RVE.Re = 2.8e2; %graphene sheet resistance ohms/m^2 (in reality just ohms)
    RVE.dvdw = 0.34e-3;   %van der waals (micrometer)
    RVE.dcut = 1.4e-3;                           %cut off tunnelling distance %micrometer
    RVE.me=9.1094*10^-31; %Kg
    RVE.lambda=1*1.602e-19;%Joules
    RVE.M=460;
    RVE.hp=6.62607004e-34; %meter^2 Kg/s
    RVE.e = 1.602e-19;                           %Coulomb
    
    %% Generate Softcore
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_M( RVE);% weibul length
    [ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_cst_agglo( RVE);%constant CNT length
    Ptrue = Ptrue';%[x y z;x y z]
    Center = Center';%[x y z;x y z]
    
    %%%Generate Hardcore GNP
    [ Pgg, gh, PPg ] = GNP_generate_strain( RVE );%Pgg (ellipse center no BC), gh (angle no BC), PPg (resistor with BC)
    Ngnp = length(gh);
    if Ngnp~=0
    Pgg = [Pgg(1,:);zeros(1,Ngnp);Pgg(2,:)]';
    end
    
    Linit = RVE.size;
    p = 1;
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(2*length(Length) + Ngnp,2*length(Length) + Ngnp);
    
    for strain = -0.006:0.006:0.006
        RVE.rr = rri*(1+strain);
        RVE.size = Linit .* [1-strain*RVE.pois,0,1+strain];
        RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
        Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
        
        %Reorientation of CNT after strain
        n=size(Center,1);
        deform = ones(n,1)*[1-strain*RVE.pois,0,1+strain];
        Center1 = Center .*deform;
        Ud = [(1-strain*RVE.pois)*cos(Angle'), 0*ones(n,1), (1+strain)*sin(Angle')];%[cos 0 sin]
        Ptrue1(1:2:2*n-1,:) = Center1 - 0.5*(Length'*ones(1,3)).*Ud;
        Ptrue1(2:2:2*n,:) = Center1 + 0.5*(Length'*ones(1,3)).*Ud;
        
        PiP=[];
        Pcenter=[];
        Plength=[];
        pt = [];
        for i = 1:n
            pt.P = Ptrue1(2*i-1:2*i,:);
            pt.center = Center1(i,:);
            pt.length = Length(i);
            [P, C, L] = PBC_2D(pt,RVE);
            PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
            Pcenter=[Pcenter C]; %%3 rows matrix with all centers
            Plength=[Plength L]; %%row with all lengths
        end
        
        Ptrue1 = PiP';
        Center1 = Pcenter';
        Length1 = Plength;
        
        %Reorientation of GNP after strain   
        Pggb=[];
        PPgb=[];
        rad_gb=[];
        Centerb=[];
        if Ngnp~=0
        deform = ones(Ngnp,1)*[1-strain*RVE.pois,0,1+strain];
        Pgg1 = Pgg .*deform;
        gh1 = atan(((1+strain)/(1-strain*RVE.pois))*tan(gh));
        for i=1:Ngnp
            xg=Pgg1(i,1);
            yg=Pgg1(i,3);
            Pg=[xg yg];
            alpha=gh1(i);
            radius_g=[RVE.ae;RVE.be;alpha];
            if (xg>=RVE.ae & xg<=RVE.size(1)-RVE.ae) & (yg>=RVE.ae & yg<=RVE.size(3)-RVE.ae)
                P7=Pg';
                cent=P7;
                Pc=[xg-RVE.ae*cos(alpha) xg+RVE.ae*cos(alpha);yg-RVE.ae*sin(alpha) yg+RVE.ae*sin(alpha)];
                Sur=pi*RVE.ae*RVE.be;
                siz=1;
            else
                [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( RVE.size(1), RVE.size(3), Pg', RVE.ae, RVE.be, alpha );
                siz=size(P7,2);
            end
            Pggb=[Pggb P7];%ellipse real center coordinates
            PPgb=[PPgb Pc];%segments representing ellipse for Resistor Network
            Centerb=[Centerb cent];%center of segments for Resistor Network
            rad_gb=[rad_gb radius_g*ones(1,siz)];
        end
        
        ng=size(PPgb,2);
        Pggb = [Pggb(1,:);zeros(1,0.5*ng);Pggb(2,:)]';
        PPgb = [PPgb(1,:);zeros(1,ng);PPgb(2,:)]';
        Centerb = [Centerb(1,:);zeros(1,0.5*ng);Centerb(2,:)]';
        
        end
        
        
        [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPgb, Centerb, Pggb, rad_gb, RVE);%Find all tunneling distances and the junctions on CNTs
        Ncnt = length(Length1);
        PP = [Ptrue1;PPgb];
        [Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
        

        if perco ==0
%             disp('// Warning: No percolation cluster exist!');
            result(1,1:3) = 0;
            fLAG(1,1:3) = NaN;
            junct(1,1:9) = NaN;
            break
        else 
            bn = zeros(length(Gn),1);
            bn(length(Gn))= 1;
            setup = [];
            setup.droptol = 1e-15;                       %tolerance
            maxit = 3e6;
            setup.type = 'ict';
            
            L1=ichol(Gn,setup);
            
            [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
            flag1
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            flag2
            Res2 = xn(length(Gn));
%             disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct(1,3*p-2:3*p) = Nj;
            p = p+1;
        end
        
    end 
    result = (result-result(1,2))./result(1,2);
    fLAG = fLAG*1;
    junct = junct*1;
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct(ii,:) = junct(1:end);
%     disp(['// Change % = ',num2str(100*result1(ii,:)),' %']);
%     disp(['// Flag = ',num2str(FLAG(ii,:)),' %']);
%     disp(['// Junct = ',num2str(Junct(ii,:)),' %']);
end
save('Hybrid_Strain_2D_Vf11.80p_G20p_agg_70p_angle_60d_rad_2L-over-100_5000_P100.mat','result1','Junct','FLAG')


parfor ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.1370; %CNT area fraction
    RVE.Vf_agg = 0.80; %percentage of agglomerate CNTs


    RVE.perc = 0.2/RVE.Vf; %percentage of CNT fraction for GNP particles

    
    RVE.agg_angle = 60; %*/-angle around agglomerate center CNT in degree
    RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
    RVE.pois=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    


    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 4; %semi major axis of GNP
    RVE.be = 0.5; %semi minor axis of GNP
    RVE.nbr = 30*10+1;%Number of points on the ellipse contour
    RVE.D = 50e-3; %CNT diameter (micrometer)
    RVE.sigma = 1e4;  %CNT intrinsic conductivity (S/meter)
    RVE.Re = 2.8e2; %graphene sheet resistance ohms/m^2 (in reality just ohms)
    RVE.dvdw = 0.34e-3;   %van der waals (micrometer)
    RVE.dcut = 1.4e-3;                           %cut off tunnelling distance %micrometer
    RVE.me=9.1094*10^-31; %Kg
    RVE.lambda=1*1.602e-19;%Joules
    RVE.M=460;
    RVE.hp=6.62607004e-34; %meter^2 Kg/s
    RVE.e = 1.602e-19;                           %Coulomb
    
    %% Generate Softcore
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_M( RVE);% weibul length
    [ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_cst_agglo( RVE);%constant CNT length
    Ptrue = Ptrue';%[x y z;x y z]
    Center = Center';%[x y z;x y z]
    
    %%%Generate Hardcore GNP
    [ Pgg, gh, PPg ] = GNP_generate_strain( RVE );%Pgg (ellipse center no BC), gh (angle no BC), PPg (resistor with BC)
    Ngnp = length(gh);
    if Ngnp~=0
    Pgg = [Pgg(1,:);zeros(1,Ngnp);Pgg(2,:)]';
    end
    
    Linit = RVE.size;
    p = 1;
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(2*length(Length) + Ngnp,2*length(Length) + Ngnp);
    
    for strain = -0.006:0.006:0.006
        RVE.rr = rri*(1+strain);
        RVE.size = Linit .* [1-strain*RVE.pois,0,1+strain];
        RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
        Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
        
        %Reorientation of CNT after strain
        n=size(Center,1);
        deform = ones(n,1)*[1-strain*RVE.pois,0,1+strain];
        Center1 = Center .*deform;
        Ud = [(1-strain*RVE.pois)*cos(Angle'), 0*ones(n,1), (1+strain)*sin(Angle')];%[cos 0 sin]
        Ptrue1(1:2:2*n-1,:) = Center1 - 0.5*(Length'*ones(1,3)).*Ud;
        Ptrue1(2:2:2*n,:) = Center1 + 0.5*(Length'*ones(1,3)).*Ud;
        
        PiP=[];
        Pcenter=[];
        Plength=[];
        pt = [];
        for i = 1:n
            pt.P = Ptrue1(2*i-1:2*i,:);
            pt.center = Center1(i,:);
            pt.length = Length(i);
            [P, C, L] = PBC_2D(pt,RVE);
            PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
            Pcenter=[Pcenter C]; %%3 rows matrix with all centers
            Plength=[Plength L]; %%row with all lengths
        end
        
        Ptrue1 = PiP';
        Center1 = Pcenter';
        Length1 = Plength;
        
        %Reorientation of GNP after strain   
        Pggb=[];
        PPgb=[];
        rad_gb=[];
        Centerb=[];
        if Ngnp~=0
        deform = ones(Ngnp,1)*[1-strain*RVE.pois,0,1+strain];
        Pgg1 = Pgg .*deform;
        gh1 = atan(((1+strain)/(1-strain*RVE.pois))*tan(gh));
        for i=1:Ngnp
            xg=Pgg1(i,1);
            yg=Pgg1(i,3);
            Pg=[xg yg];
            alpha=gh1(i);
            radius_g=[RVE.ae;RVE.be;alpha];
            if (xg>=RVE.ae & xg<=RVE.size(1)-RVE.ae) & (yg>=RVE.ae & yg<=RVE.size(3)-RVE.ae)
                P7=Pg';
                cent=P7;
                Pc=[xg-RVE.ae*cos(alpha) xg+RVE.ae*cos(alpha);yg-RVE.ae*sin(alpha) yg+RVE.ae*sin(alpha)];
                Sur=pi*RVE.ae*RVE.be;
                siz=1;
            else
                [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( RVE.size(1), RVE.size(3), Pg', RVE.ae, RVE.be, alpha );
                siz=size(P7,2);
            end
            Pggb=[Pggb P7];%ellipse real center coordinates
            PPgb=[PPgb Pc];%segments representing ellipse for Resistor Network
            Centerb=[Centerb cent];%center of segments for Resistor Network
            rad_gb=[rad_gb radius_g*ones(1,siz)];
        end
        
        ng=size(PPgb,2);
        Pggb = [Pggb(1,:);zeros(1,0.5*ng);Pggb(2,:)]';
        PPgb = [PPgb(1,:);zeros(1,ng);PPgb(2,:)]';
        Centerb = [Centerb(1,:);zeros(1,0.5*ng);Centerb(2,:)]';
        
        end
        
        
        [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPgb, Centerb, Pggb, rad_gb, RVE);%Find all tunneling distances and the junctions on CNTs
        Ncnt = length(Length1);
        PP = [Ptrue1;PPgb];
        [Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
        

        if perco ==0
%             disp('// Warning: No percolation cluster exist!');
            result(1,1:3) = 0;
            fLAG(1,1:3) = NaN;
            junct(1,1:9) = NaN;
            break
        else 
            bn = zeros(length(Gn),1);
            bn(length(Gn))= 1;
            setup = [];
            setup.droptol = 1e-15;                       %tolerance
            maxit = 3e6;
            setup.type = 'ict';
            
            L1=ichol(Gn,setup);
            
            [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
            flag1
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            flag2
            Res2 = xn(length(Gn));
%             disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct(1,3*p-2:3*p) = Nj;
            p = p+1;
        end
        
    end 
    result = (result-result(1,2))./result(1,2);
    fLAG = fLAG*1;
    junct = junct*1;
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct(ii,:) = junct(1:end);
%     disp(['// Change % = ',num2str(100*result1(ii,:)),' %']);
%     disp(['// Flag = ',num2str(FLAG(ii,:)),' %']);
%     disp(['// Junct = ',num2str(Junct(ii,:)),' %']);
end
save('Hybrid_Strain_2D_Vf13.70p_G20p_agg_80p_angle_60d_rad_2L-over-100_5000_P100.mat','result1','Junct','FLAG')


parfor ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.2105; %CNT area fraction
    RVE.Vf_agg = 0.90; %percentage of agglomerate CNTs


    RVE.perc = 0.2/RVE.Vf; %percentage of CNT fraction for GNP particles

    
    RVE.agg_angle = 60; %*/-angle around agglomerate center CNT in degree
    RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
    RVE.pois=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    


    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 4; %semi major axis of GNP
    RVE.be = 0.5; %semi minor axis of GNP
    RVE.nbr = 30*10+1;%Number of points on the ellipse contour
    RVE.D = 50e-3; %CNT diameter (micrometer)
    RVE.sigma = 1e4;  %CNT intrinsic conductivity (S/meter)
    RVE.Re = 2.8e2; %graphene sheet resistance ohms/m^2 (in reality just ohms)
    RVE.dvdw = 0.34e-3;   %van der waals (micrometer)
    RVE.dcut = 1.4e-3;                           %cut off tunnelling distance %micrometer
    RVE.me=9.1094*10^-31; %Kg
    RVE.lambda=1*1.602e-19;%Joules
    RVE.M=460;
    RVE.hp=6.62607004e-34; %meter^2 Kg/s
    RVE.e = 1.602e-19;                           %Coulomb
    
    %% Generate Softcore
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_M( RVE);% weibul length
    [ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
    %[ Ptrue, Center, Length, Angle, ~, ~, ~ ] = CNT_generate_cst_agglo( RVE);%constant CNT length
    Ptrue = Ptrue';%[x y z;x y z]
    Center = Center';%[x y z;x y z]
    
    %%%Generate Hardcore GNP
    [ Pgg, gh, PPg ] = GNP_generate_strain( RVE );%Pgg (ellipse center no BC), gh (angle no BC), PPg (resistor with BC)
    Ngnp = length(gh);
    if Ngnp~=0
    Pgg = [Pgg(1,:);zeros(1,Ngnp);Pgg(2,:)]';
    end
    
    Linit = RVE.size;
    p = 1;
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(2*length(Length) + Ngnp,2*length(Length) + Ngnp);
    
    for strain = -0.006:0.006:0.006
        RVE.rr = rri*(1+strain);
        RVE.size = Linit .* [1-strain*RVE.pois,0,1+strain];
        RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
        Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
        
        %Reorientation of CNT after strain
        n=size(Center,1);
        deform = ones(n,1)*[1-strain*RVE.pois,0,1+strain];
        Center1 = Center .*deform;
        Ud = [(1-strain*RVE.pois)*cos(Angle'), 0*ones(n,1), (1+strain)*sin(Angle')];%[cos 0 sin]
        Ptrue1(1:2:2*n-1,:) = Center1 - 0.5*(Length'*ones(1,3)).*Ud;
        Ptrue1(2:2:2*n,:) = Center1 + 0.5*(Length'*ones(1,3)).*Ud;
        
        PiP=[];
        Pcenter=[];
        Plength=[];
        pt = [];
        for i = 1:n
            pt.P = Ptrue1(2*i-1:2*i,:);
            pt.center = Center1(i,:);
            pt.length = Length(i);
            [P, C, L] = PBC_2D(pt,RVE);
            PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
            Pcenter=[Pcenter C]; %%3 rows matrix with all centers
            Plength=[Plength L]; %%row with all lengths
        end
        
        Ptrue1 = PiP';
        Center1 = Pcenter';
        Length1 = Plength;
        
        %Reorientation of GNP after strain   
        Pggb=[];
        PPgb=[];
        rad_gb=[];
        Centerb=[];
        if Ngnp~=0
        deform = ones(Ngnp,1)*[1-strain*RVE.pois,0,1+strain];
        Pgg1 = Pgg .*deform;
        gh1 = atan(((1+strain)/(1-strain*RVE.pois))*tan(gh));
        for i=1:Ngnp
            xg=Pgg1(i,1);
            yg=Pgg1(i,3);
            Pg=[xg yg];
            alpha=gh1(i);
            radius_g=[RVE.ae;RVE.be;alpha];
            if (xg>=RVE.ae & xg<=RVE.size(1)-RVE.ae) & (yg>=RVE.ae & yg<=RVE.size(3)-RVE.ae)
                P7=Pg';
                cent=P7;
                Pc=[xg-RVE.ae*cos(alpha) xg+RVE.ae*cos(alpha);yg-RVE.ae*sin(alpha) yg+RVE.ae*sin(alpha)];
                Sur=pi*RVE.ae*RVE.be;
                siz=1;
            else
                [ P7,Pc,cent,Sur ] = Periodic_BC_ellipse_New( RVE.size(1), RVE.size(3), Pg', RVE.ae, RVE.be, alpha );
                siz=size(P7,2);
            end
            Pggb=[Pggb P7];%ellipse real center coordinates
            PPgb=[PPgb Pc];%segments representing ellipse for Resistor Network
            Centerb=[Centerb cent];%center of segments for Resistor Network
            rad_gb=[rad_gb radius_g*ones(1,siz)];
        end
        
        ng=size(PPgb,2);
        Pggb = [Pggb(1,:);zeros(1,0.5*ng);Pggb(2,:)]';
        PPgb = [PPgb(1,:);zeros(1,ng);PPgb(2,:)]';
        Centerb = [Centerb(1,:);zeros(1,0.5*ng);Centerb(2,:)]';
        
        end
        
        
        [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPgb, Centerb, Pggb, rad_gb, RVE);%Find all tunneling distances and the junctions on CNTs
        Ncnt = length(Length1);
        PP = [Ptrue1;PPgb];
        [Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
        

        if perco ==0
%             disp('// Warning: No percolation cluster exist!');
            result(1,1:3) = 0;
            fLAG(1,1:3) = NaN;
            junct(1,1:9) = NaN;
            break
        else 
            bn = zeros(length(Gn),1);
            bn(length(Gn))= 1;
            setup = [];
            setup.droptol = 1e-15;                       %tolerance
            maxit = 3e6;
            setup.type = 'ict';
            
            L1=ichol(Gn,setup);
            
            [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
            flag1
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            flag2
            Res2 = xn(length(Gn));
%             disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct(1,3*p-2:3*p) = Nj;
            p = p+1;
        end
        
    end 
    result = (result-result(1,2))./result(1,2);
    fLAG = fLAG*1;
    junct = junct*1;
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct(ii,:) = junct(1:end);
%     disp(['// Change % = ',num2str(100*result1(ii,:)),' %']);
%     disp(['// Flag = ',num2str(FLAG(ii,:)),' %']);
%     disp(['// Junct = ',num2str(Junct(ii,:)),' %']);
end
save('Hybrid_Strain_2D_Vf21.05p_G20p_agg_90p_angle_60d_rad_2L-over-100_5000_P100.mat','result1','Junct','FLAG')

delete(gcp('nocreate'))
% delete(gcp('nocreate'))
% Junt = (Junct-Junct(:,1)*ones(1,size(Junct,2)))./(Junct(:,1)*ones(1,size(Junct,2)));
%disp(['// Change % = ',num2str(100*mean(result1,1)),' %']);
% disp(['// Change % = ',num2str(100*mean(Rff1)),' %']);
% disp(['// Change % = ',num2str(100*mean(Junt)),' %']);
