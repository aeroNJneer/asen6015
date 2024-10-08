clear
close all
init;
g = 32.174;         % gravitational accel

dt = 0.1;
nt = 3 * g;         % acceleration of target
Nvalues = [4];  % guidance gain
he = deg2rad(30);   % heading error
augmented = true;
dt = 0.01;
dtfine = 0.0001;

% Definition of state space system [y ydot nt]
% yddot = nt-nc     ntdot = 0
A = [0 1 0;
    0 0 1;
    0 0 0];
B = [0 -1 0]';

P = adjoint(A)
vm = 3000;      % velocity of missile ft/s
vt = 1000;      % velocity of target ft/s
vc = vm + vt;   % closing speed ft/s
R0 = 40000;     % ft
y0 = 0 ;        % no vertical component
tgo = R0/vc;    % time to go
tspan = 0:dt:tgo;   
% Simulation parameters
nt = 3 * g;         % acceleration of target
Nvalues = [3 4 5];  % guidance gain
he = deg2rad(30);   % heading error
augmented = true;
dt = 0.01;
dtfine = 0.0001;

C = [0 vm*sind(30) 0];

X0 = [0,3000*sind(30),4];
sys = ss(A,B,C,[],0.01)
lsim(sys,)

