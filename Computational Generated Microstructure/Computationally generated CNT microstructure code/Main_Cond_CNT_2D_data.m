
Number = 1; %Monte carlo simulation samples nnumber

%Vof = [0.02:0.005:0.25 0.3:0.05:0.4];
Vof = [0.245:0.005:0.25 0.3:0.05:0.4];

lVof = length(Vof);
count = zeros(lVof,3);
for jj = 1:lVof
    Total = zeros(Number,4);
    vad = Vof(jj);
    for ii = 1:Number
        rng(ii)
        format long
        RVE = [];
        RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
        
        
        
        RVE.Vf = vad;
        RVE.Vf_agg = 0.90; %percentage of agglomerate CNTs
        
        
        RVE.agg_angle = 0; %*/-angle around agglomerate center CNT in degree
        RVE.perc = 0.5; %percentage of CNT fraction for GNP particles
        RVE.agg_rad = 2*RVE.size(1)/100;%radius of each agg sphere
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
        
        Gtotal2 = 0;
        flag1 = 0;
        flag2 = 0;
        
        %% Generate Softcore
        %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_M( RVE); %Weibull length
        [ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_angle_M( RVE); %Weibull length + Alignment angle
        %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_cst_agglo( RVE); %Constant length
        Ptrue1 = Ptrue1';%[x y z;x y z]
        Center1 = Center1';%[x y z;x y z]
        
        mi = 1*RVE.dvdw;
        ma = 1*RVE.dcut;
        rri = mi + (ma-mi).*rand(length(Length1),length(Length1));
        RVE.rr = rri;
        
        [dCNT, ja] = CNT_contact_multiple(Ptrue1, Center1, Length1, RVE);%Find all tunneling distances and the junctions on CNTs
        [~, Nj, Gn, perco ] = ResNetwork_Soft_leafs_notouch_Nj(dCNT, ja, Ptrue1, RVE); %find the network resistance matrix
        
        
        if perco ==0
            %disp('// Warning: No percolation cluster exist!');
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
            %disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            
            Gtotal2 = 1/(Res2);
            %disp(['// Conductance = ',num2str(Gtotal2),' S'])
        end
        Total (ii,:) = [Gtotal2 Nj flag1 flag2];
    end
    
    kkk=0;
    for i=1:Number %change here
        if Total(i,1)~=0
            kkk=kkk+1;
        end
        Result1(i)=Total(i,1);
        Result(i)=mean(Result1(1:i));
    end
    count(jj,:) = [Vof(jj) kkk/i Result(i)];

    disp(['// Volume fraction right now = ',num2str(vad)])

save('CNT_Cond_data_2D_agg_90p_angle_0d_rad_2L-over-100_5000.mat','count')
end
