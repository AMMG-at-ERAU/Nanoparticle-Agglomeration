%%here we have the main algorithm
clear;clc
%parpool(2)
format long
Number = 1;
rng('default');
rng(1);
for ii = 1:Number
%% Define parameters of the RVE (Lx,Ly,Lz)

format long
%RandStream('mrg32k3a','seed',sum(200*clock)+ii);
RVE = [];
RVE.Vf = 0.009;                               %Volume fraction
RVE.size = [1,1,1];                       %[Lx,Ly,Lz] mimcrometer
RVE.dvdw = 0.47e-3;                           %Van der Waals distance, micrometer
RVE.me = 9.1094e-31;                         %electron mass Kg
RVE.DE = 1.0*1.602e-19;                        %Joules 
RVE.M = 460;
RVE.hp = 6.62607004e-34;                     %plank constant meter^2 Kg/s
RVE.lambda = 1.5*1.6e-19;                    %J
RVE.e = 1.602e-19;                           %Coulomb
RVE.weibull = [21 2.4];                  %weibull distribution micrometer
RVE.D = 9e-4;                               %CNT diameter, micrometer
RVE.dcut = 3e-4;                           %cut off tunnelling distance %micrometer
RVE.sigma = 1.0e4;                             %intrinsic conductivity S/m
RVE.n=[0 0 1;0 0 1;1 0 0;0 1 0;1 0 0;0 1 0]; %vectors normal to RVE planes
RVE.Dir = 'z';                               %direction of strain and current, [x y z]
RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
RVE.pois = 0.3;                              %poisson's ratio
RVE.nBin = 50;                                %number of grids/bins in each direction


%% form CNT network
[ Ptrue1, Center1, Length,~,~,~] = workinG_on(RVE);    %Generate the CNT
pipi=Ptrue1';
n=size(pipi,2);
Linit = RVE.size;


    %RandStream('mrg32k3a','seed',sum(100*clock)+ii);
    %rng('shuffle'); 
p = 1;
for strain = 0:0.006:0.006
tic
RVE.size = Linit .* [1-strain*RVE.pois,1-strain*RVE.pois,1+strain];
RVE.V=[0 0 0;RVE.size];                      %Points common to the two "three adjacents" planes of RVE
        X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0;0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0;0,0;RVE.size(1),RVE.size(1);0,0;RVE.size(1),RVE.size(1)];
        Y=[0,RVE.size(2);RVE.size(2),RVE.size(2);RVE.size(2),0;0,0;0,RVE.size(2);RVE.size(2),RVE.size(2);RVE.size(2),0;0,0;0,0;0,0;RVE.size(2),RVE.size(2);RVE.size(2),RVE.size(2)];
        Z=[0,0;0,0;0,0;0,0;RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);0,RVE.size(3);0,RVE.size(3);0,RVE.size(3);0,RVE.size(3)];

[ Ptrue, Center] = A_Reorientation( Ptrue1, Center1, Length, RVE, strain );       
PP = Ptrue'; 
n=size(PP,2);
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
[ P, C, L] = A_PBC(CNT,RVE);
PiP=[PiP P];   %% Matrix with all the CNT generated in the box with BCs
Pcenter=[Pcenter C]; %%3 rows matrix with all centers
Plength=[Plength L]; %%row with all lengths
%Pangle=[Pangle A];   %%2 rows matrix with all angles
 end
 PP=PiP;
 n=size(PP,2)
 
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

[dCNT, ja] = A_CNT_contact_multiple( PP', Pcenter', Plength, RVE); %Find all tunneling distances and the junctions on CNTs
[Gn, perco ] = A_ResNetwork_Soft_leafs(dCNT, ja, PP', RVE); %find the network resistance matrix
%[Gn, perco] = A_ResNetwork_t(dCNT, PP', RVE);  %find the network resistance matrix
if perco ==0
    disp('// Warning: No percolation cluster exist!');
    Result1(ii,:) = 0;
    disp(['// Change % = ',num2str(100*Result1(ii,:)),' %']);
    break
else
   bn = zeros(length(Gn),1); 
   bn(length(Gn))= 1;
   setup = [];
   setup.droptol = 1e-15;                       %tolerance 
   maxit = 1e3;
   [L,U]=ilu(Gn,setup);
   [xn,flag,relres,iter] = bicg(Gn,bn,1e-4,maxit,L,U);
   Res1 = xn(length(Gn));
   
   L1=ichol(Gn,setup);
   [xn,~]=pcg(Gn,bn,1e-4,maxit,L1,L1');
   Res2 = xn(length(Gn));
   disp(['Strain = ',num2str(strain*100),'%']);
   disp(['// Resistance = ',num2str(Res1),' S/m']);
   disp(['// Resistance2 = ',num2str(Res2),' S/m']);
   toc
   Result(1,p)=Res1;
   p = p+1;
   %Gtotal1=(RVE.size(3)*10^-6)/(Res1*RVE.size(1)*RVE.size(2)*10^-12)
end

end
Result = (Result-Result(1,1))./Result(1,1);
Result1(ii,:) = Result(2:end);
disp(['// Change % = ',num2str(100*Result1(ii,:)),' %']);
end
save('Strain.mat','Result1')
disp(['// Change % = ',num2str(100*mean(Result1)),' %']);
