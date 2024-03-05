mycluster1=parcluster('local');
mycluster1.NumWorkers=60;
parpool(mycluster1,60,'SpmdEnabled',false);
Number = 1000;

parfor ii = 1:Number
%% Define parameters of the RVE (Lx,Ly,Lz)
% rho_m=1.1;%g/cm^3
%     rho_f=2.1;%g/cm^3
%     Wt=3/100;
%     Vf=1/(1+(-1+1/Wt)*(rho_f/rho_m))
format long
stream = RandStream('mrg32k3a','seed',sum(clock)+200*ii);
RandStream.setGlobalStream(stream); %restart from begining
%rng(6);
RVE = [];
RVE.Vf = 0.004;%0.0159;                               %Volume fraction
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
RVE.pois = 0.4;                              %poisson's ratio
RVE.nBin = 15;                                %number of grids/bins in

%% form CNT network
[ Ptrue1, Center1, Length,Angle1, ~,~,~] = A_HardCore_multiple_Grids_small_boxp_PBC(RVE);    %Generate the CNT
pipi=Ptrue1';
n=size(pipi,2);
Linit = RVE.size;


    %RandStream('mrg32k3a','seed',sum(100*clock)+ii);
    %rng('shuffle'); 
p = 1;
for strain = [0 0.006]%-0.006:0.001:0.006
%tic
RVE.size = Linit .* [1-strain*RVE.pois,1-strain*RVE.pois,1+strain];
RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0;0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0;0,0;RVE.size(1),RVE.size(1);0,0;RVE.size(1),RVE.size(1)];
        Y=[0,RVE.size(2);RVE.size(2),RVE.size(2);RVE.size(2),0;0,0;0,RVE.size(2);RVE.size(2),RVE.size(2);RVE.size(2),0;0,0;0,0;0,0;RVE.size(2),RVE.size(2);RVE.size(2),RVE.size(2)];
        Z=[0,0;0,0;0,0;0,0;RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);0,RVE.size(3);0,RVE.size(3);0,RVE.size(3);0,RVE.size(3)];

%[ Ptrue, Center] = A_Reorientation( Ptrue1, Center1, Length, RVE, strain );       
[ Ptrue, Center] = A_Reorientation_tan_Z( Angle1, Ptrue1, Center1, Length, RVE, strain ); %Gives same result as affine transformation
PP = Ptrue'; 
n=size(PP,2);
%toc
%%

%  for i=1:n/2
%             Xx(i,:)=PP(1,2*i-1:2*i);
%             Yy(i,:)=PP(2,2*i-1:2*i);
%             Zz(i,:)=PP(3,2*i-1:2*i);
%         end
%         figure
%         plot3(X',Y',Z','r')
%         hold on
%         plot3(Xx',Yy',Zz','k')
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
%tic
 PiP=[];
 Pcenter=[];
 Plength=[];
 Pangle=[];
 CNT = [];
 for i=1:n/2
CNT.pt = Ptrue(2*i-1:2*i,:);
CNT.center = Center(i,:);
%CNT.angle = Angle(i,:);
CNT.length = Length(i);
[ P, C, L] = A_PBC1(CNT,RVE);
PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
Pcenter=[Pcenter C]; %%3 rows matrix with all centers
Plength=[Plength L]; %%row with all lengths
%Pangle=[Pangle A];   %%2 rows matrix with all angles
 end
 PP=PiP;
 n=size(PP,2)
%toc
%   for i=1:n/2
%             Xx(i,:)=PP(1,2*i-1:2*i);
%             Yy(i,:)=PP(2,2*i-1:2*i);
%             Zz(i,:)=PP(3,2*i-1:2*i);
%         end
%         figure
%         plot3(X',Y',Z','r')
%         hold on
%         plot3(Xx',Yy',Zz','k')
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
%tic
0.25*(pi*RVE.D^2)*sum(Plength)

[dCNT, ja] = A_CNT_contact_multiple_Grids_boxp( PP', Plength, RVE); %Find all tunneling distances and the junctions on CNTs
%toc
%tic
[Nj, Gn, perco ] = A_ResNetwork_Soft_leafs_notouch_Nj(dCNT, ja, PP', RVE); %find the network resistance matrix
%toc
if perco ==0
    %disp('// Warning: No percolation cluster exist!');
    break
else
   %tic
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
   [xn1,flag1]=pcg(Gn,bn,1e-4,maxit,L1,L1');
   flag1
   [xn,flag2]=pcg(Gn,bn,1e-4,maxit,L1,L1',xn1);
   flag2
   Res2 = xn(length(Gn));
   %disp(['// Resistance2 = ',num2str(Res2),' Ohm']);
   %Result(1,p)=Res1;
   result(1,p)=Res2;
   fLAG(1,p) = flag2;
   junct(1,p) = Nj;
      Gtotal1=(RVE.size(3)*10^-6)/(Res2*RVE.size(1)*RVE.size(2)*10^-12);
   Resistivity(1,p) = 1/Gtotal1;
   p = p+1;
   %toc

end

end
%Result = (Result-Result(1,1))./Result(1,1);
result = (result-result(1,1))./result(1,1);
Resistivity = (Resistivity-Resistivity(1,1))./Resistivity(1,1);
fLAG = fLAG*1;
junct = junct*1;
%Result1(ii,:) = Result(1:end);
result1(ii,:) = result(1:end);
Rvity(ii,:) = Resistivity(1:end);
FLAG(ii,:) = fLAG(1:end);
Junct(ii,:) = junct(1:end);
%disp(['// Change % = ',num2str(100*Result1(ii,:)),' %']);
% disp(['// Change % = ',num2str(100*result1(ii,:)),' %']);
% disp(['// Change % = ',num2str(100*Rvity(ii,:)),' %']);
% disp(['// Flag = ',num2str(FLAG(ii,:)),' %']);
% disp(['// Junct = ',num2str(Junct(ii,:)),' %']);
end
%tstop = toc(tstart)
save('SPIE_CNT_Strain_Vf4p_L300n_c_2n_li100n_D1n_miu_0.4_lambda_0.5.mat','result1','Rvity','FLAG','Junct')
delete(gcp('nocreate'))