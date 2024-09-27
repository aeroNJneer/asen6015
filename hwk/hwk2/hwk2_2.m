% ASEN-6015
% Hwk 1 Problem 2
close all
clear
clc
init;
mu = Gearth;
r0 = [6578 0 0];
rf = [-42164*cosd(5) 42164*sind(5) 0];

Tof = 37864;
TAf = 5;

v0 = [0 sqrt(mu/(norm(r0)^3)) 0];
vf = [sqrt(mu/abs(rf(1))) sqrt(mu/abs(rf(2))) 0];


[dv1, dv2] = lambert2(mu,r0,rf,v0,vf,Tof,"short");

