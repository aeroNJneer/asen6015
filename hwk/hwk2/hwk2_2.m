% ASEN-6015
% Hwk 1 
% Problem 2 Part a
close all
clear
clc
init;
mu = Gearth;
r0_hat = [1 0 0]';
r0 = 6578;
R0 = r0* [1 0 0]';

rf_hat = [cosd(5); sind(5); 0];
rf = 42164;
Rf = rf * rf_hat;  % assume circular orbit

Tof = 37864;

v0 = sqrt(mu/r0^3) * r0_hat;
vf = sqrt(mu/rf^3) * rf_hat;


[v1, v2] = lambert_univ(mu,R0,Rf,Tof,1)

% Problem 2 Part b
% simulate a case with simple lamber guidance routine where steering
% direction is determined by a_T = a_T v_hat_d

