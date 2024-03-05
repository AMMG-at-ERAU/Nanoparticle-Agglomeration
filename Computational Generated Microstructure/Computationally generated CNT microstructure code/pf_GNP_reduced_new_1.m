%%%% We are using ellipse with end cap added
myclusterGA=parcluster('local');
myclusterGA.NumWorkers=36;
parpool(myclusterGA,35,'SpmdEnabled',true);
Number = 35;

Dum.int = 0;

% save('GNP_Cond_Vf15p_L1000n_c_2n_a50n_b50n_t5n_N60.mat','-struct','Dum');
% save('BackUp.mat','-struct','Dum');


%Total = zeros(Number,3);
tstart0 = tic;
parfor ii = 1:Number
format long
warning('off','all')
%tstart = tic;
%RandStream('mt19937ar','seed',sum(200*clock)+ii);
stream = RandStream('mrg32k3a','seed',sum(clock)+200*ii);
RandStream.setGlobalStream(stream); %restart from begining
%rng(1);
RVE = [];
RVE.Vf = 0.48/100;                               %Volume fraction
RVE.size = 0.3*[1 1 1];                       %[Lx,Ly,Lz] mimcrometer
RVE.dvdw = 0.47e-3;                           %Van der Waals distance, micrometer
RVE.me = 9.1094e-31;                         %electron mass Kg
RVE.DE = 1.0*1.602e-19;                        %Joules 
RVE.M = 460;
RVE.hp = 6.62607004e-34;                     %plank constant meter^2 Kg/s
RVE.lambda = 0.5*1.6e-19;                    %J
RVE.e = 1.602e-19;                           %Coulomb
RVE.weibull = [21 2.4];                  %weibull distribution micrometer
RVE.D = 0.345e-3;%4e-3;                               %CNT diameter, micrometer
RVE.dcut = 2e-3;                           %cut off tunnelling distance %micrometer
RVE.sigma = 1.0e4;                             %intrinsic conductivity S/m
RVE.Re=2.8*10^2;                             %graphene sheet resistance ohms/m^2 (in reality just ohms)
RVE.n=[0 0 1;0 0 1;1 0 0;0 1 0;1 0 0;0 1 0]; %vectors normal to RVE planes
RVE.Dir = 'z';                               %direction of strain and current, [x y z]
RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
RVE.pois = 0.3;                              %poisson's ratio
RVE.nBin = 16;                                %number of grids/bins in each direction
RVE.ae = 50e-3;                                  %semi major axis of ellipse
RVE.be = 50e-3;                                  %semi mino axis
RVE.th = (100/1400)*1e-3;%0.345e-3;                                 % ellipse thickness
RVE.procs = ii;                               %processor number
xu = 0;
RVE.Lit = [-1e-6 -1e-6 -1e-6 -RVE.size-[1e-6 1e-6 1e-6]]'; % to make sure that there is enough ellipse on
%either side of boundary, with length greater than e-6
RVE.gniko = 'BackUpGNP_AR1400_0.45p_L300_%d.mat';           % name of backup

poursave_bis( sprintf(RVE.gniko, ii),Dum ); %create the file Backupii.mat

[Vector, Center, Vector1, Center1, Resistor1, PA, K1, dCNT, ja, storE] = GNP_GenerateHardCore_reduced_PBC( RVE);%the best
poursave_new( Vector, Center, Vector1, Center1, Resistor1, PA, K1, dCNT, ja, storE,RVE,sprintf('GNP_Cond_Vf0.48p_L300n_c_2n_a50n_b50n_t100over1400n_N70_%d.mat', ii) ); %save Resistor1, dCNT and ja
[Nj, Nj1, Nj2, Lj, Lj1, Lj2, Gn, perco ] = ResNetwork_GNP_Soft_leafs_notouch_Nj_Nj1( dCNT, ja, Resistor1', RVE);

Vector = []; Center = []; Vector1 = []; Center1 = []; Resistor1 = []; PA = []; K1 = []; dCNT = []; ja = []; storE = [];


end

toc(tstart0)

delete(gcp('nocreate'))