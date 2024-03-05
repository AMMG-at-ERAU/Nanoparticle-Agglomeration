function [ Resistor, Center, L, Vector, UU, Resistor1, Center1, L1, Vector1, PA, K1, Ncnt, Ncnt1, guarder ] = Hybrid_GNP_CNT_GenerateHardCore_Frac_PBC_save( RVE)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
% Here all GNPs are generated first before all CNTS
% ///////////////////////// This is working FASTEST///////////////////////////
% PP has all the CNT ends points as each row with 3 columns
% Pcenter is a 3 columns matrix with all centers
% Pangle is a 2 columns matrix with all angles
% Plength is a row vector with the CNT lengths
% dCNT is a N by N matrix with the tunneling disances
% ja is a N by N matrix with each row representing the different junction
% on the CNT of that row (their distance from the starting point of the
% CNT)
% N is the number of CNTs
% Nbin is number of bin in each axis
% ctobin is the number of bin per unit of length, assuming cubic bins(Nbin/Lx)
% ///////////////////// Grids are used //////////////////////////////
tic
%format long
proc = RVE.procs;
gniK = RVE.gniko;
choice = 1e6; %Maximum number of particles expected
Nom = sprintf(gniK, proc);


Vector1=load(Nom,'Vector1a');
if isfield(Vector1,'Vector1a')
    Vector1 = Vector1.Vector1a;
    Vector=load(Nom,'Vectora');
    Vector = Vector.Vectora;
    Resistor=load(Nom,'Resistora');
    Resistor = Resistor.Resistora;
    Center=load(Nom,'Centera');
    Center = Center.Centera;
    L = load(Nom,'La');
    L = L.La;
    Resistor1=load(Nom,'Resistor1a');
    Resistor1 = Resistor1.Resistor1a;
    Center1=load(Nom,'Center1a');
    Center1 = Center1.Center1a;
    L1 = load(Nom,'L1a');
    L1 = L1.L1a;
    MATRIX=load(Nom,'MATRIXa');
    MATRIX = MATRIX.MATRIXa;
    MATRIX_INV=load(Nom,'MATRIX_INVa');
    MATRIX_INV = MATRIX_INV.MATRIX_INVa;
    BDB=load(Nom,'BDBa');
    BDB = BDB.BDBa;
    K1=load(Nom,'K1a');
    K1 = K1.K1a;
    PA=load(Nom,'PAa');
    PA = PA.PAa;
    PA1=load(Nom,'PA1a');
    PA1 = PA1.PA1a;
    Grid=load(Nom,'Grida');
    Grid = Grid.Grida;
    guarder=load(Nom,'guardera');
    guarder = guarder.guardera;
    UU=load(Nom,'UUa');
    UU = UU.UUa;
    counT=load(Nom,'counTa');
    counT = counT.counTa;
    i = counT(1);
    ii = counT(2);
    Vpr = counT(3);
    overlap_gnp = counT(4);
    Vpr1 = counT(5);
    Overlap_cnt = counT(6);
    Ngnp = counT(7);
    Ngnp1 = counT(8);
    Tim=load(Nom,'Tima');
    Tim = Tim.Tima;
    rng(Tim) %load the random number generator settings
    
else
    
    i = 1;                                                          % useful CNT's index
    ii = 1;
    Vpr=0;
    overlap_gnp =0;
    Vpr1=0;
    Overlap_cnt = 0;
    Ngnp = 0;
    Ngnp1 = 0;
    Vector = zeros(3,3*choice);
    Resistor = zeros(3,2*choice); %Ptrue' for CNT
    Center = zeros(choice,3); %Same as Center for CNT
    L = zeros(1,choice); % Same as L for CNT
    Vector1 = zeros(3,3*choice);
    Resistor1 = zeros(3,2*choice);%Ptrue1' for CNT
    Center1 = zeros(choice,3);%Same as Center1 for CNT
    L1 = zeros(1,choice); % Same as L for CNT
    MATRIX = zeros(3,3*choice);
    MATRIX_INV = zeros(3,3*choice);
    BDB = zeros(choice,6);
    K1 = zeros(1,choice);
    PA = [];
    PA1 = [];
    Nbin = RVE.nBin;
    Grid = cell(Nbin,Nbin,Nbin); %Create cell of array with Nbin in each dimension
    guarder = zeros(choice,2);
    UU = zeros(choice,3);
end

%CNT PARAMETERS
%Weibul parameters in micrometer
aW = 5.6403;
bW = 2.4;
trio = [1 2 3];
Nbin = RVE.nBin
ctobin = Nbin./RVE.size %[ctobinx, ctobiny, ctobinz]
unit_CNT = [0 1 0;0 0 1;1 0 0]; %CNT axis is our Z, X is our Y and Y is our Z for CNT
unit = [1,0,0;0,1,0;0,0,1]; %major axis  is our X, minor axis is our Y and Z is our Z for GNP
RVE.VfCNT*prod(RVE.size)

%Plot RVE Boundaries
% X=[0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0;0,0;0,RVE.size(1);RVE.size(1),RVE.size(1);RVE.size(1),0;0,0;RVE.size(1),RVE.size(1);0,0;RVE.size(1),RVE.size(1)];
% Y=[0,RVE.size(2);RVE.size(2),RVE.size(2);RVE.size(2),0;0,0;0,RVE.size(2);RVE.size(2),RVE.size(2);RVE.size(2),0;0,0;0,0;0,0;RVE.size(2),RVE.size(2);RVE.size(2),RVE.size(2)];
% Z=[0,0;0,0;0,0;0,0;RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);RVE.size(3),RVE.size(3);0,RVE.size(3);0,RVE.size(3);0,RVE.size(3);0,RVE.size(3)];
% plot3(X',Y',Z','r')
% hold on
% xlabel('X')
% ylabel('Y')
% zlabel('Z')


%% GNP Hardcore
% Vpr=0;
% overlap_gnp =0;
while Vpr<RVE.VfGNP*prod(RVE.size)
    %Generate GNP
    r=rand(1,5);
    randP=r(1,1:3).*[1,1,1];
    
    u1 = 1.0-2.0*r(1,4);%cos theta
    v1 = sqrt(1.0-u1^2);%sin theta
    w1 = 2*pi*r(1,5);%phi
    center = RVE.size.*randP; %ellipse center point
    Matrix = [u1*cos(w1) -sin(w1) v1*cos(w1);u1*sin(w1) cos(w1) v1*sin(w1);-v1 0 u1]; %rotation matrix (Local to global)
    Matrix_inv = [u1*cos(w1) u1*sin(w1) -v1;-sin(w1) cos(w1) 0;v1*cos(w1) v1*sin(w1) u1]; %rotation matrix (global to local)
    Runit = Matrix*unit; %unit vector of ellipse after rotation [U V W]
    
    %Get the CNTs data before PBC
    Vector(:,3*ii-2:3*ii) = Runit; %ellipse base vectors
    Center(ii,:) = center; %ellipse center [x y z]
    UU(ii,:) = [u1 v1 w1];
    %     Center(1:i,:)
    %     i
    %Bounding box and PBC
    BdB = [center'-(((RVE.ae*Runit(:,1)).^2)+((RVE.be*Runit(:,2)).^2)).^0.5 ;-(center'+(((RVE.ae*Runit(:,1)).^2)+((RVE.be*Runit(:,2)).^2)).^0.5)];%[[xmin ymin zmin]';-[xmax ymax zmax]']
    temp = BdB>=RVE.Lit;
    if prod(double(temp),1) == 0 % Touching any RVE boundaries
        temp =(~temp);
        try
            [M, A, Pa, K] = PBC_elliptical( RVE, Matrix, Matrix_inv, Runit, center, temp );
        catch
            A = center;
            M = [center'-RVE.ae*Runit(:,1), center'+RVE.ae*Runit(:,1)];
            Pa = [0;0;-1/2;1];
            K = 1;
        end
    else
        A = center;
        M = [center'-RVE.ae*Runit(:,1), center'+RVE.ae*Runit(:,1)];
        Pa = [0;0;-1/2;1];
        K = 1;
    end
    sizA = size(A,1);
    
    for jj = 1:sizA
        kk =1;
        %Get the CNTs data after PBC
        Vector1(:,3*i-2:3*i) = Runit;
        Center1(i,:) = A(jj,:);
        MATRIX(:,3*i-2:3*i) = Matrix;
        MATRIX_INV(:,3*i-2:3*i) = Matrix_inv;
        vt = [A(jj,:)-((((RVE.ae*Runit(:,1)).^2)+((RVE.be*Runit(:,2)).^2)).^0.5)' (A(jj,:)+((((RVE.ae*Runit(:,1)).^2)+((RVE.be*Runit(:,2)).^2)).^0.5)')]+0.5*(RVE.th+2E-6)*[-1 -1 -1 1 1 1];%[xmin-th/2 ymin-th/2 zmin-th/2 xmax+th/2 ymax+th/2 zmax+th/2]
        BD(jj,:) = [0.5*(vt(1)+vt(4)) 0.5*(vt(2)+vt(5)) 0.5*(vt(3)+vt(6)) -vt(1)+vt(4) -vt(2)+vt(5) -vt(3)+vt(6)];%[vtcenter(x,y,z) vtlength(x,y,z)]
        
        if i> 1
            
            temp = (ones(i-1,1)*BD(jj,1:3)>=(BDB(1:i-1,1:3)-0.5*(ones(i-1,1)*BD(jj,4:6)+BDB(1:i-1,4:6)))) & (ones(i-1,1)*BD(jj,1:3)<=(BDB(1:i-1,1:3)+0.5*(ones(i-1,1)*BD(jj,4:6)+BDB(1:i-1,4:6))));
            temp = prod(double(temp),2);                                % find adjcent GNPs
            temp = temp > 0;
            if ~isempty(find(temp, 1))
                a = find(temp);
                Np = Runit(:,3);
                Ap = A(jj,:);
                for ju=1:length(a)
                    Nt = Vector1(:,3*a(ju));
                    At = Center1(a(ju),:);
                    [P,N] = A_Plane_Plane_intersect(Np',Ap,Nt',At);%two points on intersecting line
                    % bring points to local 1
                    Mp = [P' N'] - Ap'*ones(1,2);%vectors from ellipse center
                    Mp = Matrix_inv*Mp; %each point is a column
                    P1 = A_ellipse_Line( Mp(1:2,:), RVE.ae, RVE.be );
                    if ~isempty(P1) %line cut ellipse
                        % bring the 2 points to global
                        Mp = Matrix*P1; %each point is a column
                        Mp = Mp + Ap'*ones(1,2); %add vector from ellipse center
                        % bring points to local 2
                        Mt = [Mp P' N'] - At'*ones(1,4);%vectors from ellipse center
                        Mt = MATRIX_INV(:,3*a(ju)-2:3*a(ju))*Mt; %each point is a column
                        M0 = Mt(1:2,1:2);
                        crit = (M0(1,:)/RVE.ae).^2 + (M0(2,:)/RVE.be).^2 - 1 < 0;%inside ellpse
                        if ~isempty(find(crit, 1))%point is inside so there is overlap
                            i = i - 1;
                            ii = ii -1;
                            kk = 0;
                            overlap_gnp=overlap_gnp+1;
                            sizA = 1;
                            %Pt = [];
                            break;
                        end
                        
                        P1 = A_ellipse_Line( Mt(1:2,3:4), RVE.ae, RVE.be );
                        if ~isempty(P1)
                            % bring the 2 points to global
                            Mt = MATRIX(:,3*a(ju)-2:3*a(ju))*P1; %each point is a column
                            Mt = Mt + At'*ones(1,2); %add vector from ellipse center
                            % bring points to local 1
                            Mp = Mt - Ap'*ones(1,2);%vectors from ellipse center
                            Mp = Matrix_inv*Mp; %each point is a column
                            M0 = Mp(1:2,:);
                            crit = (M0(1,:)/RVE.ae).^2 + (M0(2,:)/RVE.be).^2 - 1 < 0;%inside ellpse
                            if ~isempty(find(crit, 1))%point is inside so there is overlap
                                i = i - 1;
                                ii = ii -1;
                                kk = 0;
                                overlap_gnp=overlap_gnp+1;
                                sizA = 1;
                                %Pt = [];
                                break;
                            end
                        end
                    end
                end
                if kk == 0
                    break;
                end
                
                for j=1:length(a)
                    %[ dt,pt1,pt2 ] = DistBtwEllipse( RVE.ae, RVE.be, (Center1(i,:))', (Center1(a(j),:))', Vector1(:,3*i-2:3*i), Vector1(:,3*a(j)-2:3*a(j)) );
                    [ dt,~,~,~,~ ] = DistBtwEllipse_Hybrid( RVE, (Center1(i,:))', (Center1(a(j),:))', Vector1(:,3*i-2:3*i), Vector1(:,3*a(j)-2:3*a(j)), M(:,2*jj-1:2*jj), Resistor1(:,2*a(j)-1:2*a(j)), Pa(:,sum(K(1:jj-1))+1:sum(K(1:jj))), PA(:,sum(K1(1:a(j)-1))+1:sum(K1(1:a(j)))) );
                    if dt < 2E-6 %error on distance
                        i = i - 1;
                        ii = ii -1;
                        kk = 0;
                        overlap_gnp=overlap_gnp+1;
                        sizA = 1;
                        break;
                    end
                end
                
                if kk == 0
                    break;
                end
            end
        end
    end
    
    for jj1 = 1:sizA
        if kk == 1
            Vector1(:,3*i-2:3*i) = Runit;
            Center1(i,:) = A(jj1,:);
            Resistor1(:,2*i-1:2*i) = M(:,2*jj1-1:2*jj1);
            PA = [PA Pa(:,sum(K(1:jj1-1))+1:sum(K(1:jj1)))];
            K1(i) = K(jj1);
            BDB(i,:) = BD(jj1,:);
            guarder(i,:)=[i ii]; %[i ii]
        end
        i = i + 1;
    end
    
    Vp = kk*pi*RVE.ae*RVE.be*RVE.th; %volume of one elliptical GNP
    Vpr = Vpr + Vp;
    ii = ii + 1;
    
    Ngnp = ii-1;
    Ngnp1 = i-1;
    
    if ((kk ==1) && mod(ii,100)==0)
%         ii
%         i
disp(['// where we are in the GNP part = ',num2str(ii),' second ']);
        T.Vectora = Vector;
        T.Resistora = Resistor;
        T.Centera = Center;
        T.La = L;
        T.Vector1a = Vector1;
        T.Resistor1a = Resistor1;
        T.Center1a = Center1;
        T.L1a = L1;
        T.MATRIXa = MATRIX;
        T.MATRIX_INVa = MATRIX_INV;
        T.BDBa = BDB;
        T.K1a = K1;
        T.PAa = PA;
        T.PA1a = PA1;
        T.Grida = Grid;
        T.guardera = guarder;
        T.UUa = UU;
        T.counTa = [i,ii,Vpr,overlap_gnp,Vpr1,Overlap_cnt,Ngnp,Ngnp1];
        T.Tima = rng;
        save(Nom,'-struct','T','-append');
    end
    %plot
    %     k=linspace(0,2*pi);
    %     XX = center'*ones(1,length(k)) + RVE.ae*Runit(:,1)*cos(k) + RVE.be*Runit(:,2)*sin(k);
    %     plot3(XX(1,:),XX(2,:),XX(3,:),'k','LineWidth',2)
    %     fill3(XX(1,:),XX(2,:),XX(3,:),'b')
    %     hold on
end

%% CNT Hardcore

while Vpr1<RVE.VfCNT*prod(RVE.size)
    % Generate CNT
    r=rand(1,5);
    %r=rand(1,6);
    randP=r(1,1:3).*[1,1,1];
    
    center = RVE.size.*randP;
    %li = 0.4; %micrometer
    %li=aW*(-log(1-r(1,6)))^(1/bW);
    
    u1 = 1.0-2.0*r(1,4);
    v1 = sqrt(1.0-u1^2);
    w1 = 2*pi*r(1,5);
    Matrix = [u1*cos(w1) -sin(w1) v1*cos(w1);u1*sin(w1) cos(w1) v1*sin(w1);-v1 0 u1];
    base = 0.5*[RVE.D RVE.D -RVE.D -RVE.D;RVE.D -RVE.D RVE.D -RVE.D;0 0 0 0]; %[X;Y;Z]four points around center with radius
    Runit = Matrix*unit_CNT; %unit vector of ellipse after rotation [U V W]
    
    P1=center-0.5*RVE.li*[v1*cos(w1),v1*sin(w1),u1];
    P2=center+0.5*RVE.li*[v1*cos(w1),v1*sin(w1),u1];
    %Get the CNTs data before PBC
    Vector(:,3*ii-2:3*ii) = Runit;
    iii = ii - Ngnp;
    Resistor(:,2*iii-1:2*iii) = [P1;P2]';
    Center(ii,:) = center;
    L(1,iii) = RVE.li;
    UU(ii,:) = [u1 v1 w1];
    
    list = cell(1,1);
    CNT.pt = Resistor(:,2*iii-1:2*iii)';
    CNT.center = Center(ii,:);
    CNT.length = L(1,iii);
    % PBC
    [ P, center1, li1] = A_PBC1(CNT,RVE);
    P = P';center1 = center1';
    D = [];
    K = [];
    Pa = [];
    sizC = size(center1,1);
    for jj =1:sizC
        k = 1;
        K = [K 3];
        Pa = [Pa [1 0 0.5*li1(jj) NaN]' [-1 0 0.5*li1(jj) 0]' [0 1 0 0]'];
        %Get the CNTs data after PBC
        Vector1(:,3*i-2:3*i) = Runit;
        Center1(i,:) = center1(jj,:);
        iv = i - Ngnp1;
        L1(1,iv) = li1(jj);
        Resistor1(:,2*i-1:2*i) = P(2*jj-1:2*jj,:)';
        vt = [min(P(2*jj-1:2*jj,:)) max(P(2*jj-1:2*jj,:))] + 0.5*RVE.D*[-1 -1 -1 1 1 1];
        BD(jj,:) = [0.5*(vt(1)+vt(4)) 0.5*(vt(2)+vt(5)) 0.5*(vt(3)+vt(6)) -vt(1)+vt(4) -vt(2)+vt(5) -vt(3)+vt(6)];%[vtcenter(x,yz) vtlength(x,y,z)]
        X2=(L1(1,iv))*[v1*cos(w1),v1*sin(w1),u1]'; %length of CNT for CNT vector
        Z = zeros(3,8);
        Z(:,1:2:7) = Matrix*base + (Resistor1(:,2*i-1))*ones(1,4);%four starting points of the 4 CNTs all around (Matrix*base find components in reference frame)
        Z(:,2:2:8) = Z(:,1:2:7) + X2*ones(1,4);%four end points of the 4 CNTs all around
        Z = Z'; %[X Y Z]
        
        
        Ptrue2 = (ones(8,1)*ctobin).*Z;%to scale everything for each bin with length 1 unit, to avoid rounding errors
        
        %Spatial patitioning
        Min = min(Ptrue2);%[xmin ymin zmin]
        Max = max(Ptrue2);%[xmax ymax zmax]
        %find the bins this filler may overlap
        Minn = floor(Min)+1;
        Minn = Minn + (Minn < 1);
        Maxn = floor(Max)+1;
        Maxn = Maxn - (Maxn > Nbin);
        
        B = floor(Ptrue2)+1;
        B = B + (B < 1);
        B = B - (B > Nbin);
        Ptrue3(:,:,1) = Ptrue2(1:2,:);
        Ptrue3(:,:,2) = Ptrue2(3:4,:);
        Ptrue3(:,:,3) = Ptrue2(5:6,:);
        Ptrue3(:,:,4) = Ptrue2(7:8,:);
        %%%%%% Find bins filler intersects
        for gi = trio
            n1 = unit(gi,:);
            V0 = (Minn(gi):Maxn(gi)-1)' * n1;
            % V0 = [1 0 0;2 0 0;3 0 0;4 0 0]
            I = A_bin_line_plane_multiple(n1,V0,Ptrue3(1,:,:),Ptrue3(2,:,:));
            tr = trio([1:gi-1 gi+1:end]);
            b = floor(I(:,tr,:))+1;
            A = zeros(2*size(b,1),3,4);
            at = (Minn(gi):Maxn(gi)-1)';
            at = at(:,:,ones(4,1));
            A(1:2:end,[gi tr],:)= [at b];
            at = (Minn(gi)+1:Maxn(gi))';
            at = at(:,:,ones(4,1));
            A(2:2:end,[gi tr],:)= [at b];
            B =[B;A(:,:,1);A(:,:,2);A(:,:,3);A(:,:,4)];
        end
        C = unique(B,'rows');
        C(~prod((C>0),2),:)=[];
        C(~prod((C<Nbin + 1),2),:)=[];
        D{jj} = C;
        
        for ij = 1:size(C,1)
            list{1} = [list{1} Grid{C(ij,1),C(ij,2),C(ij,3)}]; %list of all CNTs in overlapped bins
        end
        list2 = unique(cell2mat(list{1}));
        
        if ~isempty(list2)
            dt = A_DistBetween2Segment_multiple(Resistor1(:,2*i-1)',Resistor1(:,2*i)',Resistor1(:,2*list2-1)',Resistor1(:,2*list2)',RVE.dvdw,RVE.D);
            if prod(dt>0)==0%prod(dt>=RVE.dvdw)==0
                i = i-1;
                ii = ii - 1;
                Overlap_cnt = Overlap_cnt + 1;
                k = 0;
                sizC = 1;
                break;
            end
        end
        
        if Ngnp1>0 %iG = 1 at the begining withoug any GNP
            temp = (ones(Ngnp1,1)*BD(jj,1:3)>=(BDB(1:Ngnp1,1:3)-0.5*(ones(Ngnp1,1)*BD(jj,4:6)+BDB(1:Ngnp1,4:6)))) & (ones(Ngnp1,1)*BD(jj,1:3)<=(BDB(1:Ngnp1,1:3)+0.5*(ones(Ngnp1,1)*BD(jj,4:6)+BDB(1:Ngnp1,4:6))));
            temp = prod(double(temp),2);                                % find adjcent GNPs
            temp = temp > 0;
            if ~isempty(find(temp, 1))
                a = find(temp);
                for j=1:length(a)
                    [ dt,~,~,~,~ ] = DistBtwEllipse_Hybrid( RVE, (Center1(a(j),:))', (Center1(i,:))', Vector1(:,3*a(j)-2:3*a(j)), Vector1(:,3*i-2:3*i), Resistor1(:,2*a(j)-1:2*a(j)), P(2*jj-1:2*jj,:)', PA(:,sum(K1(1:a(j)-1))+1:sum(K1(1:a(j)))), Pa(:,sum(K(1:jj-1))+1:sum(K(1:jj))) );
                    if dt < 2E-6 %error on distance
                        i = i - 1;
                        ii = ii -1;
                        k = 0;
                        Overlap_cnt = Overlap_cnt + 1;
                        sizC = 1;
                        break;
                    end
                end
                
                if k == 0
                    break;
                end
            end
        end
        
    end
    
    for jj1 =1:sizC
        
        Vp = 0;
        if k == 1 %either distance between CNTs greatr than van der waal or iii>=100
            Vector1(:,3*i-2:3*i) = Runit;
            Center1(i,:) = center1(jj1,:);
            L1(1,iv) = li1(jj1);
            Resistor1(:,2*i-1:2*i) = P(2*jj1-1:2*jj1,:)';
            BDB(i,:) = BD(jj1,:);
            PA1 = [PA1 Pa(:,sum(K(1:jj1-1))+1:sum(K(1:jj1)))];
            K1(i) = K(jj1);
            C = D{jj1};
            guarder(i,:)=[i ii]; %[i ii]
            Vp=0.25*(pi*RVE.D^2)*L1(iv);
            for ij = 1:size(C,1)
                Grid{C(ij,1),C(ij,2),C(ij,3)} = [Grid{C(ij,1),C(ij,2),C(ij,3)} {i}];
            end
            
        end
        Vpr1=Vpr1+Vp;
        i = i + 1;
        iv = i - Ngnp1;
    end
    ii = ii + 1;
    
    if ((k ==1) && mod(ii,500)==0)
%         ii
%         i
disp(['// where we are in the CNT part = ',num2str(ii),' second '])
        T.Vectora = Vector;
        T.Resistora = Resistor;
        T.Centera = Center;
        T.La = L;
        T.Vector1a = Vector1;
        T.Resistor1a = Resistor1;
        T.Center1a = Center1;
        T.L1a = L1;
        T.MATRIXa = MATRIX;
        T.MATRIX_INVa = MATRIX_INV;
        T.BDBa = BDB;
        T.K1a = K1;
        T.PAa = PA;
        T.PA1a = PA1;
        T.Grida = Grid;
        T.guardera = guarder;
        T.UUa = UU;
        T.counTa = [i,ii,Vpr,overlap_gnp,Vpr1,Overlap_cnt,Ngnp,Ngnp1];
        T.Tima = rng;
        save(Nom,'-struct','T','-append');
    end
end
%since we preallocate choice slots of data points in the following
%variable, the non used slots will be zero. We need to take those away
iii = ii - Ngnp;
Resistor = Resistor(:,1:2*(iii-1));
L = L(1,1:iii-1); %Those are only for CNTs
iv = i - Ngnp1;
L1 = L1(1,1:iv-1);%Those are only for CNTs
Ncnt1 = length(L1)
Ncnt = length(L);
%Reorder so we have CNT first and GNP at the end
Vector = [Vector(:,3*Ngnp+1:3*(ii-1)) Vector(:,1:3*Ngnp)];
Center = [Center(Ngnp+1:ii-1,:);Center(1:Ngnp,:)];
UU = [UU(Ngnp+1:ii-1,:);UU(1:Ngnp,:)];
Vector1 = [Vector1(:,3*Ngnp1+1:3*(i-1)) Vector1(:,1:3*Ngnp1)];
Center1 = [Center1(Ngnp1+1:i-1,:);Center1(1:Ngnp1,:)];
Resistor1 = [Resistor1(:,2*Ngnp1+1:2*(i-1)) Resistor1(:,1:2*Ngnp1)];
K1 = [K1(1,Ngnp1+1:i-1) K1(1,1:Ngnp1)];
PA = [PA1 PA];
guarder = [guarder(Ngnp1+1:i-1,:)-ones(Ncnt1,1)*[Ngnp1 Ngnp];guarder(1:Ngnp1,:)+ones(Ngnp1,1)*[Ncnt1 Ncnt]];
Ncnt1 = size(Center1,1)-Ngnp1
% k=linspace(0,2*pi);
% for iii =Ncnt + 1:size(Center1,1)
%     XX = (Center1(iii,:))'*ones(1,length(k)) + RVE.ae*Vector1(:,3*iii-2)*cos(k) + RVE.be*Vector1(:,3*iii-1)*sin(k);
%     plot3(XX(1,:),XX(2,:),XX(3,:),'k','LineWidth',2)
%     fill3(XX(1,:),XX(2,:),XX(3,:),'r')
%     hold on
% end
%
% for iii =1:Ncnt
% plot3(Resistor1(1,2*iii-1:2*iii),Resistor1(2,2*iii-1:2*iii),Resistor1(3,2*iii-1:2*iii),'-k','LineWidth',2)
% end

Vpri = Vpr +Vpr1
Overlap_cnt
overlap_gnp
disp(['// Time for GenerateHardCore = ',num2str(toc),' second ']);
end

