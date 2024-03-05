clc
clear all
format long
Nbin = 6;
Lx = 30;
Ly = 30;
Lz = 30;
trio = [1 2 3];
ctobin = Nbin/Lx; %if Lx = Ly = Lz
% P0 = ctobin*[0 0 0]; %to scale everything for each bin with length 1 unit
% P1 = ctobin*[30 30 30]; %to scale everything for each bin with length 1 unit, to avoid rounding errors
P0 = ctobin*[1 0 0]; %to scale everything for each bin with length 1 unit
P1 = ctobin*[30 0 0]; %to scale everything for each bin with length 1 unit, to avoid rounding errors

Ptrue1 = [P0;P1];
Min = min(Ptrue1);%[xmin ymin zmin]
Max = max(Ptrue1);%[xmax ymax zmax]
        %find the bins this filler may overlap
        Minn = floor(Min)+1;
        Minn = Minn + (Minn < 1)
        Maxn = floor(Max)+1;
        Maxn = Maxn - (Maxn > Nbin)
%n = ctobin*[Lx/Nbin,0,0;0,Ly/Nbin,0;0,0,Lz/Nbin]
n = [1,0,0;0,1,0;0,0,1]
B =floor(Ptrue1)+1;
B = B + (B < 1);
B = B - (B > Nbin);
% B =[];
for gi = trio
n1 = n(gi,:);
V0 = (Minn(gi):Maxn(gi)-1)' * n1
% V0 = [1 0 0;2 0 0;3 0 0;4 0 0]
I= A_bin_line_plane(n1,V0,P0,P1)
tr = trio([1:gi-1 gi+1:end])
lu= I(:,tr)
floor(lu)
b = floor(I(:,tr))+1
A = zeros(2*size(b,1),3)
[gi tr]
A(1:2:end,[gi tr])= [[Minn(gi):Maxn(gi)-1]' b];
A(2:2:end,[gi tr])= [[Minn(gi)+1:Maxn(gi)]' b];
B =[B;A]
end
C = unique(B,'rows')
size(C,1)