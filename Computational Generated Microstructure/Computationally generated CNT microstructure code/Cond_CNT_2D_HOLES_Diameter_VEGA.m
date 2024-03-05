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
    
    
    
    RVE.Vf = 0.20;
    RVE.radius_hole = (1:1:5)*RVE.size(1)/20;
    Area_frac_hole = pi*RVE.radius_hole.^2/(RVE.size(1)*RVE.size(3))
    RVE.Vf_holes = Area_frac_hole; %area fraction of holes
    RVE.Vf_agg = 0; %percentage of agglomerate CNTs
    
    
    RVE.agg_angle = 10; %*/-angle around agglomerate center CNT in degree
    RVE.perc = 0.5; %percentage of CNT fraction for GNP particles
    RVE.agg_rad = 12*RVE.size(1)/100;%radius of each agg sphere
    RVE.miu=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    
    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 4; %semi major axis of GNP
    RVE.be = 0.5; %semi minor axis of GNP
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
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(length(Length1),length(Length1));
    RVE.rr = rri;
    
    for damage = 0:1:length(RVE.Vf_holes)
        
        if (damage ==0)
            [dCNT, ja] = CNT_contact_multiple(Ptrue1, Center1, Length1, RVE);%Find all tunneling distances and the junctions on CNTs
        else
            Vf = RVE.Vf_holes(damage);
            Radius = RVE.radius_hole(damage);
            hole = Generate_holes(RVE, Vf, Radius);
            [dCNT, ja] = CNT_contact_multiple_holes(Ptrue1, Center1, Length1, hole, Radius, RVE);%Find all tunneling distances and the junctions on CNTs
        end
        
%[~, Nj, Gn, perco ] = ResNetwork_Soft_leafs_notouch_Nj(dCNT, ja, Ptrue1, RVE); %find the network resistance matrix
        [~, Nj, Nj1, Nj2, Gn, perco ] = ResNetwork_Soft_leafs_notouch_Nj_Nj1(dCNT, ja, Ptrue1, RVE); %find the network resistance matrix
        
        
        if perco ==0
            %disp('// Warning: No percolation cluster exist!');
            result(1,1:1+length(RVE.Vf_holes)) = 0;
            fLAG(1,1:1+length(RVE.Vf_holes)) = NaN;
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
            Res2 = xn(length(Gn))
            %disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            
            %Gtotal2 = 1/(Res2);
            %disp(['// Conductance = ',num2str(Gtotal2),' S'])
            result(1,p)=Res2;
            fLAG(1,p) = flag2;
            junct(1,p) = Nj;
            junct1(1,p) = Nj1;
            junct2(1,p) = Nj2;
            p = p+1;
        end
        
    end
    Data(ii) = RVE;
    result = result*1;
    fLAG = fLAG*1;
    junct = junct*1;
    junct1 = junct1*1;
    junct2 = junct2*1;
    
    result1(ii,:) = result(1:end);
    FLAG(ii,:) = fLAG(1:end);
    Junct(ii,:) = junct(1:end);
    Junct1(ii,:) = junct1(1:end);
    Junct2(ii,:) = junct2(1:end);
end
Data = Data(1);
save('CNT_Hole_Cond_2D_Vf_20p_agg_0p_Hole_rad_(1_1_5)Lover20_number_1_5000.mat','result1','Junct','Junct1','Junct2','FLAG','Data')

delete(gcp('nocreate'))