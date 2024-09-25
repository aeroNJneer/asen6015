% HWK 1 problem 3.3
clc
close all
clear
init;

%% finite burn
r1 = 200 + Rearth;
dv1 = 2.4546;
dv2 = 1.4773;
Tof = 1.8932e+04;
v1 = sqrt(Gearth/r1);

% [r,r_dot, theta, theta_dot, h, h_dot]
accel = 0.03;                   % constant acceleration 30 m/s^2
dv = [0, dv1, 0]';        % first burn vector
btime1 = norm(dv)/norm(accel);
btime2 = norm(dv2)/norm(accel);
r0 = [r1, 0, 0]';
v0 = [0, v1, 0]';  % at periapsis

% solve ODEs
tspan = linspace(0,btime1,1000);
x1 = [r0;v0/r1];
[tfin,yfin] = ode45(@(tfin,yfin) finiteburn(tfin,yfin,dv,accel), tspan, [r0 ;v0]);

% propogate to apoapsis
% coasting after finite burn using previous ode
tspan = linspace(0,Tof,1000);
r2 =[yfin(1000,1); yfin(1000,2); yfin(1000,3)]; % initial conditions
v2 = [yfin(1000,4); yfin(1000,5); yfin(1000,6)]; 
[t2fin,y2fin] = ode45(@(t2fin,y2fin) finiteburn(t2fin,y2fin,0,0), tspan, [r2 ;v2]);

tfinite = [tfin; tfin(end)+t2fin];
yfinite = [yfin;y2fin];
X = [];
Q = [];

% Extract position components
for i = 1:length(yfinite)
    [raan inc theta] = RSW2Euler(yfinite(i,1:3),yfinite(i,4:6));
    C = Euler3132C([raan inc theta]);
    X = [X C*yfinite(i,1:3)'];
end

% Plot the 3D orbit
figure(1);
plot3(X(1,:),X(2,:),X(3,:));
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Orbit');

figure(2)
stackedplot(tfinite,X(1:3,:)')

%%
function dydt = finiteburn(t,y,deltaV,accel)
    % Constants
    init;
    mu = Gearth; % Gravitational parameter for Earth (km^3/s^2)
    r = y(1:3);
    v = y(4:6);
    DV = deltaV/norm(deltaV);   % direction of deltav
    Accel = accel*v/norm(v);           % acceleration along deltaV vector
    r_hat = r/norm(r);          % radial direction
    % EOM
    r_dot = v;
    v_dot = (-mu/(norm(r)^3)).*r_hat  + Accel;
    dydt = [r_dot; v_dot];
end
