function [ Ptrue, CENTER ] = A_Reorientation_tan_Z( U, Ptrue, Center, Length, RVE, strain )
%UNTITLED Summary of this function goes here
% Gives same result as affine transformation
%   Detailed explanation goes here
%   Center 3 columns vectors
%   Ptrue 3 columns vectors
Ud = U;
N = size(Center,1);
deform = ones(N,1) * [1-strain*RVE.pois,1-strain*RVE.pois,1+strain]; 
CENTER = Center.*deform;
Ud(:,3) = atan2(sin(U(:,3)),cos(U(:,3)));
Phi = atan2((1-strain*RVE.pois)*U(:,2).*cos(U(:,3))./cos(Ud(:,3)),(1+strain)*U(:,1));
Ud(:,1) = cos(Phi);
Ud(:,2) = sin(Phi);

Ptrue(1:2:2*N-1,:) = CENTER - 0.5*(Length'*ones(1,3)).*[Ud(:,2).*cos(Ud(:,3)),Ud(:,2).*sin(Ud(:,3)),Ud(:,1)];
Ptrue(2:2:2*N,:) = CENTER + 0.5*(Length'*ones(1,3)).*[Ud(:,2).*cos(Ud(:,3)),Ud(:,2).*sin(Ud(:,3)),Ud(:,1)];
end