clc
clear all
Nbin = 6;
Lx = 0.6;
Ly = 0.6;
Lz = 0.6;
trio = [1 2];
ctobin = Nbin/Lx;
% P0 = [0.5 0.5 0]/10;
% P1 = [4.5 0.5 0]/10;
% Nbin = 5;
% Lx = 0.5;
% Ly = 0.5;
% Lz = 0.5;
% P0 = [0.1 4.8 0]/10;
% P1 = [3 3 0]/10;
% Nbin = 5;
% Lx = 0.5;
% Ly = 0.5;
% Lz = 0.5;

% P0 = [1.5 4.5 0]/10;
% P1 = [4.8 2.8 0]/10;
Ptrue1 = [P0;P1];
Min = min(Ptrue1);%[xmin ymin zmin]
Max = max(Ptrue1);%[xmax ymax zmax]
        %find the bins this filler may overlap
        Minn = floor(Min*ctobin)+1
        Minn = Minn + (Minn < 1);
        Maxn = floor(Max*ctobin)+1
        Maxn = Maxn - (Maxn > Nbin);
n = [Lx/Nbin,0,0;0,Ly/Nbin,0;0,0,Lz/Nbin];
B =[];
for i = trio
n1 = n(i,:);
V0 = (Minn(i):Maxn(i)-1)' * n1
% V0 = [1 0 0;2 0 0;3 0 0;4 0 0]
I= A_bin_line_plane(n1,V0,P0,P1)
tr = trio([1:i-1 i+1:end])
b = floor(I(:,tr)*ctobin)+1
A = zeros(2*length(b),2)
[i tr]
A(1:2:end,[i tr])= [[Minn(i):Maxn(i)-1]' b];
A(2:2:end,[i tr])= [[Minn(i)+1:Maxn(i)]' b];
B =[B;A]
end
C = unique(B,'rows')