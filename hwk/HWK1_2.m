clc
close all
clear

init
mu = Gearth;
a = 20000;
e = 0.01;
inc = deg2rad(0.1);
raan = 0;
w = 0;
n = MMfromMua(mu,a);
M = @(M0,n,t) M0 + n*t;
dM = @(a,n,t) -3*n*t/(2*a);
h = HfromAe(mu,a,e);
p = PfromAe(a,e);
t = 0:10:86400;
r0 = a;
r = r0;
TA = 0;
A = @(t,M0) [0 0 0 0 0 0
    0 0 0 0 0 0 
    0 0 0 0 0 0 
    0 0 0 0 0 0 
    0 0 0 0 0 0
    -dM(a,n,t) 0 0 0 0 0];

stm0 = expm(A(1000)*10)
xd = A(1000)*expm(A(1000)*10)
% create reference path
Mref = M(0,n,t);    % mean anomaly over time
% TAref = TAfrom;     % true anomaly over time
U = @(v_r,v_theta,phi_dot,v_z) [v_r; v_theta/r+sin(phi_dot); v_z];
Xref = zeros(6,length(t));
Xref(1,:) = a;
Xref(2,:) = e;
Xref(3,:) = inc;
Xref(4,:) = raan;
Xref(5,:) = w;
Xref(6,:) = Mref;

dB = [(2*a^2/h)*(e*sin(TA)) (2*a^2/h)*p/r 0
    p*sin(TA)/h (p+r)*cos(TA)+(r*e)   0
    0 0 (r*cos(TA+w))/h
    0 0 r*sin(TA+w)/(h*sin(inc))
    -p*cos(TA)/(h*e) (p+r)*sin(TA)/(h*e) -r*sin(TA+w)*cos(inc)/(h*sin(inc))
    n+a*sqrt(1-e^2)/(a*h*e)*(p*cos(TA)-(2*r*e)) a*sqrt(1-e^2)/(a*h*e)*(p+r)*sin(TA) 0];



% dMt = dM(a,n,t);
% U0 = U(0,0.1,0,0);
% dA0 = A(1,Mref(1));
% 
% C = [1 0 0 0 0 0; 0 0 0 0 0 1; 0 0 0 0 0 0 ];
% sys = ss(dA0,dB,C,[]);
% tsim = 0:0.1:100;
% Ut = zeros(length(tsim)-1,3);
% Uf = [U0' ; Ut];
% 
% lsim(ss,U0,tsim)