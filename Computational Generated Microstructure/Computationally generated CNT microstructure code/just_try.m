%%here we have the main algorithm
clear;clc
%parpool(16)
format long
% rng('default');
rng(1);
tstart = tic;

%% Define parameters of the RVE (Lx,Ly,Lz)
% rho_m=1.1;%g/cm^3
%     rho_f=2.1;%g/cm^3
%     Wt=3/100;
%     Vf=1/(1+(-1+1/Wt)*(rho_f/rho_m))
format long
%RandStream('mrg32k3a','seed',sum(200*clock)+ii);
RVE = [];
RVE.Vf = 0.0159;                               %Volume fraction
RVE.size = [2,2,2];                       %[Lx,Ly,Lz] mimcrometer
RVE.dvdw = 0.47e-3;                           %Van der Waals distance, micrometer
RVE.me = 9.1094e-31;                         %electron mass Kg
RVE.DE = 1.0*1.602e-19;                        %Joules 
RVE.M = 460;
RVE.hp = 6.62607004e-34;                     %plank constant meter^2 Kg/s
RVE.lambda = 2*1.6e-19;                    %J
RVE.e = 1.602e-19;                           %Coulomb
RVE.weibull = [21 2.4];                  %weibull distribution micrometer
RVE.D = 4e-3;                               %CNT diameter, micrometer
RVE.dcut = 1.8e-3;                           %cut off tunnelling distance %micrometer
RVE.sigma = 1.0e4;                             %intrinsic conductivity S/m
RVE.n=[0 0 1;0 0 1;1 0 0;0 1 0;1 0 0;0 1 0]; %vectors normal to RVE planes
RVE.Dir = 'z';                               %direction of strain and current, [x y z]
RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
RVE.pois = 0.3;                              %poisson's ratio
RVE.nBin = 15;                                %number of grids/bins in each direction (0.2*Lcnt)
xu = 0;

%% form CNT network
[ Resistor, Center, Length,Resistor1, Center1, L1] = A_HardCore_multiple_Grids_small_boxp_PBC(RVE);    %Generate the CNT
% rng(1);
% [ Ptrue2, Center2, Length2,Ptrue3, Center3, L3 ] = A_HardCore_BdB_PBC( RVE);
% pipi=Ptrue1';
% n=size(pipi,2);
% Linit = RVE.size;

