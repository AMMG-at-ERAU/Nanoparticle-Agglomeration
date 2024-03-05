Number = 1;%Monte carlo simulation samples nnumber

count = zeros(1,6);
Total = zeros(Number,5);

for ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    
    
    RVE.Vf = 0.0930; %CNT area fraction
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
    
    Gtotal2 = 0;
    flag1 = 0;
    flag2 = 0;
    
    %% Generate Softcore
    %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_M( RVE); %Weibull length
    [ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
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
    
    [ dCNT, ja ] = Hybrid_contact( Ptrue1, Center1, Length1, PPg, CenterPPg, Pgg, rad_g, RVE);%Find all tunneling distances and the junctions on CNTs
    Ncnt = length(Length1);
    PP = [Ptrue1;PPg];
    [Nj, Gn, perco ] = ResNetwork_Hybrid_Soft_leafs_notouch_Nj( Ncnt, dCNT, ja, PP, RVE); %find the network resistance matrix
    
    
    if perco ==0
        disp('// Warning: No percolation cluster exist!');
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
        disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
        
        Gtotal2 = 1/(Res2);
        %disp(['// Conductance = ',num2str(Gtotal2),' S'])
    end
    Total (ii,:) = [Gtotal2 flag2 Nj];% [G flag2 Njcnt Njgnp Njall]
end

    kkk=0;
    for i=1:Number %change here
        if Total(i,1)~=0
            kkk=kkk+1;
        end
        Result1(i)=Total(i,1);
        Result(i)=mean(Result1(1:i));
    end
    count(1,:) = [kkk/i mean(Total,1)];

%save('Hybrid_Cond_2D_Vf9.30p_G20p_agg_10p_angle_60d_rad_2L-over-100_5000_P100.mat','count')

