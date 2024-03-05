function [ Rt ] = ResTun1( RVE, d )
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
% d               tunneling distance with CNT diameter substracted, micrometer
% RVE.me          electron mass Kg
% RVE.hp          planck constant meter^2 Kg/s
% RVE.D           CNT diameter, micrometer
% RVE.lambda      height of barrier for epoxy J
% RVE.e           electron charge Coulomb
% Rt              Tunelling resistance in Ohms

RVE.D = RVE.D*1e-6;      %from micrometer to meter
d = d*1e-6;        %from micrometer to meter
A = 0.25*pi*RVE.D^2;     %CNT cross sectional area (m^2)
Rt = ((d*RVE.hp^2)/(A*RVE.e^2*sqrt(2*RVE.me*RVE.lambda))).*exp((4*pi*d/RVE.hp)*sqrt(2*RVE.me*RVE.lambda));
end

