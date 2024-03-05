%%here we have the main algorithm
clear all;clc

Number = 1;

Total = zeros(Number,4);
tstart = tic;
for ii = 1:Number
%% Define parameters of the RVE (Lx,Ly,Lz)
% rho_m=1.1;%g/cm^3
%     rho_f=2.1;%g/cm^3
%     Wt=3/100;
%     Vf=1/(1+(-1+1/Wt)*(rho_f/rho_m))
format long
% stream = RandStream('mrg32k3a','seed',sum(clock)+200*ii);
% RandStream.setGlobalStream(stream); %restart from begining
rng(2);
RVE = [];
RVE.Vf = 0.005;%0.0159;                               %Volume fraction
RVE.size = [0.3,0.3,0.3];                       %[Lx,Ly,Lz] mimcrometer
RVE.dvdw = 0.47e-3;                           %Van der Waals distance, micrometer
RVE.me = 9.1094e-31;                         %electron mass Kg
RVE.DE = 1.0*1.602e-19;                        %Joules 
RVE.M = 460;
RVE.hp = 6.62607004e-34;                     %plank constant meter^2 Kg/s
RVE.lambda = 0.5*1.6e-19;                    %J
RVE.e = 1.602e-19;                           %Coulomb
RVE.weibull = [21 2.4];                  %weibull distribution micrometer
RVE.D = 1e-3;                               %CNT diameter, micrometer
RVE.li = 0.1;                               %CNT length, micrometer
RVE.dcut = 2e-3;                           %cut off tunnelling distance %micrometer
RVE.sigma = 1.0e4;                             %intrinsic conductivity S/m
RVE.n=[0 0 1;0 0 1;1 0 0;0 1 0;1 0 0;0 1 0]; %vectors normal to RVE planes
RVE.Dir = 'z';                               %direction of strain and current, [x y z]
RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
RVE.pois = 0.3;                              %poisson's ratio
RVE.nBin = 15;                                %number of grids/bins in each direction (0.2*Lcnt)
xu = 0;

Gtotal2 = 0;
flag1 = 0;
flag2 = 0;
%% form CNT network
[ ~,~,~,~,Ptrue, Center, Length] = A_HardCore_multiple_Grids_small_boxp_PBC(RVE);    %Generate the CNT

% %Plot CNTs
% n=size(Ptrue,1);
% X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0;0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0;0,0;RVE.size(1),RVE.size(1);0,0;RVE.size(1),RVE.size(1)];
% Y=[0,RVE.size(2);RVE.size(2),RVE.size(2);RVE.size(2),0;0,0;0,RVE.size(2);RVE.size(2),RVE.size(2);RVE.size(2),0;0,0;0,0;0,0;RVE.size(2),RVE.size(2);RVE.size(2),RVE.size(2)];
% Z=[0,0;0,0;0,0;0,0;RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);0,RVE.size(3);0,RVE.size(3);0,RVE.size(3);0,RVE.size(3)];
% 
% for i=1:n/2
%     Xx(i,:)=Ptrue(2*i-1:2*i,1)';
%     Yy(i,:)=Ptrue(2*i-1:2*i,2)';
%     Zz(i,:)=Ptrue(2*i-1:2*i,3)';
% end
% figure
% plot3(X',Y',Z','r')
% hold on
% plot3(Xx',Yy',Zz','k')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')

tic
[dCNT, ja] = A_CNT_contact_multiple_Grids_boxp( Ptrue, Length, RVE); %Find all tunneling distances and the junctions on CNTs

toc
%%
tic
[Nj, Gn, perco ] = A_ResNetwork_Soft_leafs_notouch_Nj(dCNT, ja, Ptrue, RVE); %find the network resistance matrix
toc

if perco ==0
    disp('// Warning: No percolation cluster exist!');
else
   bn = zeros(length(Gn),1); 
   bn(length(Gn))= 1;
   setup = [];
   setup.droptol = 1e-15;                       %tolerance 
   maxit = 3e6;
   setup.type = 'ict';
%    [L,U]=ilu(Gn,setup);
%    [xn,flag,relres,iter] = bicg(Gn,bn,1e-2,maxit,L,U);
%    disp(['flag ',num2str(flag),'%']);
%    Res1 = xn(length(Gn));
   
   L1=ichol(Gn,setup);
   xl = xu*ones(length(Gn),1);
   [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
   flag1
   [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
   flag2
   Res2 = xn(length(Gn));
   disp(['// Resistance2 = ',num2str(Res2),' Ohm']);

   Gtotal2=(RVE.size(3)*10^-6)/(Res2*RVE.size(1)*RVE.size(2)*10^-12);
   disp(['// Conductivity = ',num2str(Gtotal2),' S/m'])
end
Total (ii,:) = [Gtotal2 Nj flag1 flag2];
end
tstop = toc(tstart)
%save('Strain.mat','Result1')
