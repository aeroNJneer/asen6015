% hwk 4

clear
clc
close all
init;
set(0,'DefaultAxesXGrid','on', ...
    'DefaultAxesYGrid','on', ...
    'DefaultLineLineWidth',2, ...
    'DefaultAxesLineWidth',1.5);
set(groot, 'defaultAxesFontSize',12, ...
    'defaultTextFontSize',12, ...
    'defaultLegendFontSize',14, ...
    'defaultAxesTitleFontSizeMultiplier', 1.5, ...
    'defaultAxesTitleFontSizeMultiplier', 1.2);

% function 
theta = @(Ax,Az,Bx,Bz,t) atan2(Az+Bz*t,Ax+Bx*t);
alpha = @(Ax,Ay,Bx,By,t) atan2(Ay+By*t,Ax+Bx*t);
uhat = @(A,B,t) (A + B*t)/norm(A+B*t);

% initial conditions
A0 = [1;0;0];
B0 = [0;0;0];
r0 = [8276; 5612; 5];
v0 = [-3.142; 4.672; 0];
X0 = [ r0; v0];
t0 = 0;

% final conditions
af = 10000;
i = 0;
ecc = 0.0001
w = deg2rad(35);

athrust = 0.003; %   km/s^2
tfguess =  23;
rperi = Rp(af,ecc);
rapo = Ra(af,ecc)

[ai,ecci,ii,Raani,wi,true_anomaly] = ECI2classical(r0,v0,Gearth)

a_delta = af - ai;
ecc_delta = ecc - ecci;
w_delta = w - wi;

