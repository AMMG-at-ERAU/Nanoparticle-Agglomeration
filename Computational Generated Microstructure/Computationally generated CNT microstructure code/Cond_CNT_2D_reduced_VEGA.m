
clear;clc

Number = 1;

Total = zeros(Number,4);
parfor ii = 1:Number
    rng(ii)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    RVE.Vf = 0.12; %CNT area fraction
    RVE.perc = 0.5; %percentage of CNT fraction for GNP particles
    RVE.Vf_agg = 0; %percentage of agglomerate CNTs
    RVE.agg_rad = RVE.size(1)/48;%radius of each agg sphere
    RVE.agg_angle = 10; %*/-angle around agglomerate center CNT in degree
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
    
    %% Generate Softcore
    [ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo( RVE); %Weibull length
    %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_cst_agglo( RVE); %Constant length
    Ptrue1 = Ptrue1';%[x y z;x y z]
    Center1 = Center1';%[x y z;x y z]
    
    mi = 1*RVE.dvdw;
    ma = 1*RVE.dcut;
    rri = mi + (ma-mi).*rand(length(Length1),length(Length1));
    RVE.rr = rri;
    
    [dCNT, ja] = CNT_contact_multiple(Ptrue1, Center1, Length1, RVE);%Find all tunneling distances and the junctions on CNTs
    poursave_CNT( dCNT, ja, Ptrue1, Length1, RVE, sprintf('CNT_Cond_2D_Vf20p_L25000n_c_1.4n_Wbli5000n_D50n_%d.mat', ii) );
   
end




