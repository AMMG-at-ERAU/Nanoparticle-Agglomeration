% mycluster1=parcluster('local');
% mycluster1.NumWorkers=60;
% parpool(mycluster1,60,'SpmdEnabled',false);
clc
clear all
Number = 1;

Total = zeros(Number,5);
for ii = 1:Number
    rng(ii+1)
    format long
    RVE = [];
    RVE.size = [25,0,25];                       %[Lx,Ly,Lz] mimcrometer
    
    RVE.Vf = 0.10; %CNT area fraction
    RVE.Vf_agg = 0.90; %percentage of agglomerate CNTs
    RVE.agg_angle = 0; %*/-angle around agglomerate center CNT in degree

    RVE.perc = 0.2/RVE.Vf; %percentage of CNT fraction for GNP particles

    

    RVE.agg_rad = 1*RVE.size(1)/100;%radius of each agg sphere
    RVE.pois=0.3;   %poisson's ratio of polymer
    RVE.li = 5;  %CNT length (micrometer)
    
    RVE.n= [0 0 1;1 0 0]; %vectors normal to four RVE planes (nz, nx)
    RVE.Dir = 'z';                               %direction of strain and current, [x y z]
    RVE.V= [0 0 0;RVE.size];      %Points common to all four planes of 2D RVE
    RVE.Wa = 5.6403; %weibul a (micrometer)
    RVE.Wb = 2.4; %weibul b (micrometer)
    RVE.ae = 1; %semi major axis of GNP
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
    
    %% Generate Softcore CNT
    [ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo( RVE);
    %[ ~, ~, ~, ~, Ptrue1, Center1, Length1 ] = CNT_generate_agglo_angle( RVE); %Weibull length
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
    
    X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0];
Y=[0,RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),0;0,0];
figure
plot(X',Y','r')
hold on

% k=linspace(0,2*pi);
% for i=1:Ngnp
% xx=RVE.ae*sin(k)*cos(rad_g(3,i))-RVE.be*cos(k)*sin(rad_g(3,i))+Pgg(i,1);
% yy=RVE.ae*sin(k)*sin(rad_g(3,i))+RVE.be*cos(k)*cos(rad_g(3,i))+Pgg(i,3);
% plot(xx,yy,'r')                      %%ellipse
% hold on
% end
    for i=1:Ncnt
    Xx(i,:)=Ptrue1(2*i-1:2*i,1);
    Yy(i,:)=Ptrue1(2*i-1:2*i,3);
    end
    plot(Xx',Yy','k')
    
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
        flag1
        [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
        flag2
        Res2 = xn(length(Gn));
        disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
        
        Gtotal2 = 1/(Res2);
        disp(['// Conductance = ',num2str(Gtotal2),' S'])
    end
    Total (ii,:) = [Gtotal2 flag2 Nj];% [G flag2 Njcnt Njgnp Njall]
end
%save('Hybrid_Cond_2D_Vf4.864p_agg_60p_rad_25L-over-100_15000.mat','Total')

%delete(gcp('nocreate'))



