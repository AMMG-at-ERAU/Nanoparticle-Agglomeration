% mycluster1=parcluster('local');
% mycluster1.NumWorkers=60;
% parpool(mycluster1,60,'SpmdEnabled',false);
% Number = 1000;
clear all; clc;
for ii = 1:1
format long
stream = RandStream('mrg32k3a','seed',sum(clock)+200*ii);
RandStream.setGlobalStream(stream); %restart from begining
%rng(2);
RVE = [];
RVE.Vf = 0.025;%0.0159;                               %Volume fraction
RVE.size = 10*[1,1,1];                       %[Lx,Ly,Lz] mimcrometer
RVE.dvdw = 0.47e-3;                           %Van der Waals distance, micrometer
RVE.me = 9.1094e-31;                         %electron mass Kg
RVE.DE = 1.0*1.602e-19;                        %Joules 
RVE.M = 460;
RVE.hp = 6.62607004e-34;                     %plank constant meter^2 Kg/s
RVE.lambda = 0.5*1.6e-19;                    %J
RVE.e = 1.602e-19;                           %Coulomb
RVE.weibull = [21 2.4];                  %weibull distribution micrometer
RVE.D = 10e-3;                               %CNT diameter, micrometer
RVE.li = 1;                               %CNT length, micrometer
RVE.dcut = 2e-3;                           %cut off tunnelling distance %micrometer
RVE.sigma = 1.0e4;                             %intrinsic conductivity S/m
RVE.n=[0 0 1;0 0 1;1 0 0;0 1 0;1 0 0;0 1 0]; %vectors normal to RVE planes
RVE.Dir = 'z';                               %direction of strain and current, [x y z]
RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
RVE.pois = 0.3;                              %poisson's ratio
RVE.nBin = 50;                                %number of grids/bins in each direction (0.2*Lcnt)


%% form CNT network
[ ~,~,~,~,Ptrue1, Center1, Length1] = A_HardCore_multiple_Grids_small_boxp_PBC(RVE);    %Generate the CNT
[dCNT, ja] = A_CNT_contact_multiple_Grids_boxp( Ptrue1, Length1, RVE); %Find all tunneling distances and the junctions on CNTs
poursave_CNT( dCNT, ja, Ptrue1, Length1, RVE,sprintf('CNT_Cond_Vf2p_L10000n_c_2n_li1000n_D10n_N60_%d.mat', ii) );
[Nj, Nj1, Nj2, Lj, Lj1, Lj2, Gn, perco ] = A_ResNetwork_Soft_leafs_notouch_Nj_Nj1(dCNT, ja, Ptrue1, RVE);
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
            [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
            Res2 = xn(length(Gn));
            %disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
            
            Gtotal2=(RVE.size(3)*10^-6)/(Res2*RVE.size(1)*RVE.size(2)*10^-12)
            %disp(['// Conductivity = ',num2str(Gtotal2),' S/m'])
        end
%poursave_CNT( dCNT, ja, Ptrue1, Length1, RVE,sprintf('CNT_Cond_Vf5p_L1000n_c_2n_li100n_D5n_N60_%d.mat', ii) );
%[Nj, Gn, perco ] = A_ResNetwork_Soft_leafs_notouch_Nj(dCNT, ja, Ptrue1, RVE); %find the network resistance matrix

end
delete(gcp('nocreate'))