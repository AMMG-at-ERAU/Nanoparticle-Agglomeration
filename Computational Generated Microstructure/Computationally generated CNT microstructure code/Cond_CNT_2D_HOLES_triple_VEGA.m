myclusterHoleCondStriple=parcluster('local');
myclusterHoleCondStriple.NumWorkers=36;
parpool(myclusterHoleCondStriple,35,'SpmdEnabled',true);

Number = 5000;


parfor ii = 1:Number
    
    format long
    stream = RandStream('mrg32k3a','seed',sum(clock)+200*ii);
    RandStream.setGlobalStream(stream); %restart from begining
    
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.12; %CNT area fraction
    RVE.radius_hole = RVE.size(1)/20;
    Area_frac_hole = pi*RVE.radius_hole^2/(RVE.size(1)*RVE.size(3))
    RVE.Vf_holes = Area_frac_hole*(2:2:20); %area fraction of holes
    RVE.Vf_agg = 0; %percentage of agglomerate CNTs
    
    
    RVE.perc = 0.13/RVE.Vf; %percentage of CNT fraction for GNP particles
    
    
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
    
    p = 1;
    
    %% Generate Softcore
    [ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_M( RVE); %Weibull length
    %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
    %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_cst_agglo( RVE); %Constant length
    Ptrue1 = Ptrue1';%[x y z;x y z]
    Center1 = Center1';%[x y z;x y z]
    
    %%%Generate Hardcore GNP
    [ Pgg,PPg,CenterPPg,Ae,rad_g ] = GNP_generate( RVE ); %Pgg (ellipse center), PPg (resistor), CenterPPg (resistor center), Ae (area of ellipse), rad_g (ae,be,alpha)
    ng=size(PPg,2)
    Ngnp = size(Pgg,2)
    if ng~=0
        Pgg = [Pgg(1,:);zeros(1,0.5*ng);Pgg(2,:)]';
        PPg = [PPg(1,:);zeros(1,ng);PPg(2,:)]';
        CenterPPg = [CenterPPg(1,:);zeros(1,0.5*ng);CenterPPg(2,:)]';
    end
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(length(Length1) + 0.5*ng,length(Length1) + 0.5*ng);
    RVE.rr = rri;
    Ncnt = length(Length1);
    
    for damage = 0:1:length(RVE.Vf_holes)
        
        if (damage ==0)
            [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPg, CenterPPg, Pgg, rad_g, RVE);%Find all tunneling distances and the junctions on CNTs
        else
            Vf = RVE.Vf_holes(damage);
            Radius = RVE.radius_hole;
            hole = Generate_holes(RVE, Vf, Radius);
            [ dCNT, ja ] = Hybrid_contact_holes( Ptrue1, Center1, Length1, PPg, CenterPPg, Pgg, rad_g, hole, Radius, RVE)
        end
        
        PP = [Ptrue1;PPg];
        %[Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
        [Nj, Nj1, Nj2, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj_Nj1( Ncnt, dCNT, ja, PP, RVE);
        
        if perco ==0
            result(1,1:1+length(RVE.Vf_holes)) = 0;
            fLAG(1,1:1+length(RVE.Vf_holes)) = NaN;
            junct_CNT(1,1:1+length(RVE.Vf_holes)) = NaN;
            junct1_CNT(1,1:1+length(RVE.Vf_holes)) = NaN;
            junct2_CNT(1,1:1+length(RVE.Vf_holes)) = NaN;
            junct_GNP(1,1:1+length(RVE.Vf_holes)) = NaN;
            junct1_GNP(1,1:1+length(RVE.Vf_holes)) = NaN;
            junct2_GNP(1,1:1+length(RVE.Vf_holes)) = NaN;
            junct(1,1:1+length(RVE.Vf_holes)) = NaN;
            junct1(1,1:1+length(RVE.Vf_holes)) = NaN;
            junct2(1,1:1+length(RVE.Vf_holes)) = NaN;
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
            %flag1
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            %flag2
            Res2 = xn(length(Gn));
            %         disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            %
            %         Gtotal2 = 1/(Res2);
            %disp(['// Conductance = ',num2str(Gtotal2),' S'])
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct_CNT(1,p) = Nj(1);
            junct1_CNT(1,p) = Nj1(1);
            junct2_CNT(1,p) = Nj2(1);
            junct_GNP(1,p) = Nj(2);
            junct1_GNP(1,p) = Nj1(2);
            junct2_GNP(1,p) = Nj2(2);
            junct(1,p) = Nj(3);
            junct1(1,p) = Nj1(3);
            junct2(1,p) = Nj2(3);
            p = p+1;
            
        end
    end
    
    Data(ii) = RVE;
    result = result*1;
    fLAG = fLAG*1;
    junct_CNT = junct_CNT*1;
    junct1_CNT = junct1_CNT*1;
    junct2_CNT = junct2_CNT*1;
    junct_GNP = junct_GNP*1;
    junct1_GNP = junct1_GNP*1;
    junct2_GNP = junct2_GNP*1;
    junct = junct*1;
    junct1 = junct1*1;
    junct2 = junct2*1;
    
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct_CNT(ii,:) = junct_CNT(1:end);
    Junct1_CNT(ii,:) = junct1_CNT(1:end);
    Junct2_CNT(ii,:) = junct2_CNT(1:end);
    Junct_GNP(ii,:) = junct_GNP(1:end);
    Junct1_GNP(ii,:) = junct1_GNP(1:end);
    Junct2_GNP(ii,:) = junct2_GNP(1:end);
    Junct(ii,:) = junct(1:end);
    Junct1(ii,:) = junct1(1:end);
    Junct2(ii,:) = junct2(1:end);
    
    
end
Data = Data(1);
save('Hybrid_Hole_Cond_2D_Vf12p_G13p_agg_0p_Hole_rad_Lover20_number_(2_2_20)_5000.mat','result1','Junct_CNT','Junct1_CNT','Junct2_CNT','Junct_GNP','Junct1_GNP','Junct2_GNP','Junct','Junct1','Junct2','FLAG','Data')

delete(gcp('nocreate'))



