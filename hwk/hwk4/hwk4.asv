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

% LVLH is aligned with desired orbital plane 
zhat = @(r) -r/norm(r);                          % ECI to LVLH
yhat = @(r,v) -cross(r,v)/norm(cross(r,v));      % ECI to LVLH @ tf
xhat = @(y,z) cross(y,z)/norm(cross(y,z));       % ECI to LVLH

lvlh2eci = @(z,y,x) [x y z];

mu = Gearth;
% initial conditions
A0 = [1;0;0];
B0 = [0;0;0];
r0 = [8276; 5612; 5];
v0 = [-3.142; 4.672; 0];
X0 = [ r0; v0];
t0 = 0;
hi = cross(r0,v0);

% final conditions
af = 10000;
incf = 0;
eccf = 0.0001;
wf = deg2rad(35);
hf = HfromAe(mu,af,eccf);

athrust = 0.03; %   km/s^2
tfguess =  23;

% LONf = LONfromH(hf);

rperf = Rp(af,eccf);                 % radius of periapsis
vperf = VfromRA(mu,rperf,af);       % velocity at periapsis
rapo = Ra(af,eccf);

[ai,ecci,inci,Raani,wi,true_anomaly] = ECI2classical(r0,v0,Gearth)
rperi = Rp(ai,ecci);
vperi = VfromRA(mu,rperi,ai);
a_delta = af - ai;
ecc_delta = eccf - ecci;
% w_delta = wf - wi;
h_delta = hf - norm(hi);
vp_delta = vperf-vperi;
rp_delta = rperf - rperi;

r = r0;
v = v0;
delu = 0.1;
ts = linspace(0,tfguess,24);
for i = 1:length(ts)
    A_eci = lvlh2eci(zhat(r), yhat(r,v), xhat(yhat(r,v),zhat(r)))  * A0;
    lamba_v = A_eci + B_eci*t;
    
end




%% functions

function r = finitediff(E,u,delu)
    
    dE = (E(u+delu) - E(u))/delu;

end
