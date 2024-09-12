
close all
clear
init;
mu = 398600; % gr62)avitational constant
omegaE = 2*pi/86400;  % angular velocity of earth
rE = 6378; % radius of earth

Tof = 1.533641249191222e+04;
r = rE + 200;
rf = 35786;
t1 = PeriodfromA(Gearth,r)';
tf = PeriodfromA(Gearth,rf)';
dt = 1;
t = 0:dt:tf;

% initial conditions
x0 = 6678;
y0 = 0;
r0 = sqrt(x0^2+y0^2);
x0dot = 0;
y0dot = r0 * sqrt(mu/r0^3);
omega = sqrt(mu/r0^3); % angular velocity of satellite
X = [];
X(:,1) = [r0; 0; 0; omega];
dv1 = 2.3357;
dv2 = 1.4784;
u = zeros(length(t),2);

u(1,2) = dv1;
u(int32(Tof)/dt,2) = dv2;


Atilde = @(t,r) [0 1 0 0; 
    mu*(2*(cos(omega*t))^2-(sin(omega*t))^2)/r^3 0 3*mu*sin(omega*t)*cos(omega*t)/r^3 0;
    0 0 0 1;
    3*mu*sin(omega*t)*cos(omega*t)/r^3 0 mu*(2*(sin(omega*t)^2 -(cos(omega*t))^2))/r^3 0];
A0 = [0 1 0 0 ; 2*mu/r0^3 0 0 0 ; 0 0 0 1; 0 0 -mu/r0^3 0];
Btilde = @(r)[0 0; 1 0; 0 0; 0 1/r];
C = [0 1 0 0; 0 0 1 0];
D = zeros(2,2);
for i= 2:length(t)
    A = Atilde(t(i-1),X(1,i-1));
    B = Btilde(X(1,i-1));
    sys = ss(A,B,C,D);
    dx[i] = A*X(:,i-1) + B*u(i--1,:);
    lsim(sys,u(i-1,:),t)

end
