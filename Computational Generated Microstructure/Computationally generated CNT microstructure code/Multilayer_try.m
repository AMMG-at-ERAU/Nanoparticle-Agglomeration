clc
clear all

Number = 1;     %Monte carlo simulation samples nnumber

layer = 1;

RVEsize = cell(layer, 1);
endpoints = cell(layer, 1);
centerpoints = cell(layer, 1);
lengthCNT = cell(layer, 1);

count = zeros(1,2);

    for ii = 1:Number
        for  i = 1:layer
        %rng(ii)
        format long
        RVE = [];
        RVE.size = [25,25,25];                       %[Lx,Ly,Lz] mimcrometer
        RVE.size1 = [25,25,0.3]; 
        RVEsize{i} = RVE.size;
        RVE.layer = i;
        
        
        RVE.Vf = 0.05;
        RVE.Vf_agg = 0.0; %percentage of agglomerate CNTs
        
        
        RVE.agg_angle = 0; %*/-angle around agglomerate center CNT in degree
        RVE.perc = 0; %percentage of CNT fraction for GNP particles
        RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
        RVE.miu=0.3;   %poisson's ratio of polymer
        RVE.li = 2;  %CNT length (micrometer)
        
        RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
        RVE.Dir = 'z';                               %direction of strain and current, [x y z]
        RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
        RVE.Wa = 2; %weibul a (micrometer)
        RVE.Wb = 2; %weibul b (micrometer)
        RVE.ae = 4; %semi major axis of GNP
        RVE.be = 0.5; %semi minor axis of GNP
        RVE.D = 25e-3; %CNT diameter (micrometer)
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
%         if i==1
%[  Ptrue, Center, Length, Angle, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_M_3_percent_0deg( RVE);
        %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_M( RVE); %Weibull length
        %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
         %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_cst_agglo_multilayer_use_this_0deg( RVE); %Constant length
        [  Ptrue, Center, Length, Angle, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_angle( RVE);
        
        Ptrue1= Ptrue1';%[x y z;x y z]
        Center1 = Center1';%[x y z;x y z]
        Length1 = Length1';
        siz(i) = length(Center1);
        size(i) = length(Ptrue1);

        end
        
    end
    
        
Ptrue1 = cell2mat(endpoints);
Center1 = cell2mat(centerpoints);
Length1 = cell2mat(lengthCNT);

%% Generate Softcore
        %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_M( RVE); %Weibull length
        %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
        [ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_cst_agglo( RVE); %Constant length
        Ptrue1 = Ptrue1';%[x y z;x y z]
        Center1 = Center1';%[x y z;x y z]
        
        Ncnt = length(Length1);
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
        Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
        figure(2)
        plot(X',Y','r')
        hold on

        for i=1:Ncnt
        Xx(i,:)=Ptrue1(2*i-1:2*i,1);
        Yy(i,:)=Ptrue1(2*i-1:2*i,3);
        end
        plot(Xx',Yy','k')
        axis([0 25 0 25])