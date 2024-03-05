function [ Ptrue, CENTER ] = A_Reorientation( Ptrue, Center, Length, RVE, strain )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   Center 3 columns vectors
%   Ptrue 3 columns vectors

N = size(Center,1);
deform = ones(N,1) * [1-strain*RVE.pois,1-strain*RVE.pois,1+strain]; 
CENTER = Center.*deform;
deform = ones(2*N,1) * [1-strain*RVE.pois,1-strain*RVE.pois,1+strain];
PTRUE = Ptrue.*deform; 
LENGTH = PTRUE(1:2:2*N-1,:)-PTRUE(2:2:2*N,:);
LENGTH = sqrt(sum(LENGTH.*LENGTH,2));
Ptrue(1:2:2*N-1,:) = CENTER -(CENTER-PTRUE(1:2:2*N-1,:)).*((Length./LENGTH')'*ones(1,3));
Ptrue(2:2:2*N,:) = CENTER +(PTRUE(2:2:2*N,:)-CENTER).*((Length./LENGTH')'*ones(1,3));
end