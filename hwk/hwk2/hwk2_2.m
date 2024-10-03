% ASEN-6015
% Hwk 1 
% Problem 2 Part a
close all
clear
clc
init;
dt = 10;
mu = Gearth;
g = 9.08665;        % accel due to gravity

r0_hat = [1 0 0];
r0 = 6578;
R0 = r0* r0_hat';
V0 = VfromRA(mu,R0,R0);
rf_hat = [-cosd(5); sind(5); 0];
rf = 42164;
Rf = rf * rf_hat;  % final position

Tof = 37864;
tspan = [0 : dt : Tof];

[v1, v2,iter] = lambert_univ(mu,R0,Rf,Tof,1);
iter
% Problem 2 Part b
% simulate a case with simple lambert guidance routine where steering
% direction is determined by a_T = a_T v_hat_d

% iterate through the simulation, constantly updating r0 until we get to rf
ai = 0.03;      % constant acceleration 30 m/s
r = R0;
v = v1;
tspan = 0:dt:Tof;
tof = Tof;
burn = 0;
%v2_tol = v2+1e-2*v2;
for i = 1:length(tspan)
    vg = v2 - v1;       % velocity to go

    % Calculate gravitational force
    r_mag = norm(r);
    F_gravity = 0;     % -g * R0/r_mag; 

    % Apply transfer burns
    if abs(vg(1)) < 1e-4 ||  abs(vg(2)) < 1e-4 || abs(vg(3)) < 1e-4
       a_burn = -ai*vg/norm(vg);     % in direction of vg
       burn = burn + 1;             % count the burn time
    else
       a_burn = [0; 0; 0];
    end

    % Update velocity and position
    a_total = F_gravity + a_burn;
    v1 = v1 + a_total * dt;
    r = r + v1 * dt;
    tof = tof - dt;                 % decrement time remaining
    [v1, v2] = lambert_univ(mu,r,Rf,tof,1);
    % Store position and velocity
    R(:, i) = r;
    V(:, i) = v1;
end
burn
[x,y] = circularOrbit(mu,r0,length(tspan));

% Problem 2 Part 3
[v1, v2,~] = lambert_univ(mu,R0,Rf,Tof,1);

G = computeGravityGradient(rf, mu); % Compute gravity-gradient matrix
% Q = -(eye(3) + (tf-t)^2*G/2)/(tf-t);
% reset variables
r = R0;
v = v1;
tof = Tof;
% iterate through the simulation, constantly updating r0 until we get to rf
for i = 1:length(tspan)
    vg = v2 - v1;       % velocity to go 

    % Calculate gravitational force
    r_mag = norm(r);
    F_gravity = 0;     % -g * R0/r_mag; 

    % Apply transfer burns
    % when we get close to the target delta-v then no more burns until
    % error grows too big OR get sufficiently close to target 
    if abs(vg(1)) < 0.1 ||  abs(vg(2)) < 0.1 || abs(vg(3)) < 0.1

        a_burn = -ai*vg/norm(vg);     % in direction of vg
        Q = -1/(dt*i - tof) * (eye(3) + (dt*i - tof)^2 * G / 2); % Compute Q matrix
        a_burn = Q * a_burn; % Update acceleration vector using Q matrix
       burn = burn + 1;             % count the burn time
    else
       a_burn = [0; 0; 0];
    end

    % Update velocity and position
    a_total = F_gravity + a_burn;

    v1 = v1 + a_total * dt;
    r = r + v1 * dt;
    tof = tof - dt;                 % decrement time remaining
    [v1, v2] = lambert_univ(mu,r,Rf,tof,1);
    % Store position and velocity
    R1(:, i) = r;
    V1(:, i) = v1;
end

delR = R1(:,end) - Rf

figure(3);
plot(R(1, :), R(2, :), R1(1,:),R1(2,:));
xlabel('X (m)');
ylabel('Y (m)');
title('Orbits Solving Lamberts Problem') ;
grid on;
hold on;
plot(x,y)
legend('Lambert guidance','Q-guidance', 'initial orbit','Location','southwest')
hold off;
xlim([-60000 40000]);
ylim([-60000 40000]);

% plot from nominal

figure (4)
plot(tspan(1:end-1),vecnorm(V1(:,1:end-1))-vecnorm(V(:,1:end-1)))
title('Difference in \Deltav betweenn Lambert with Q-guidance')

% compute circular orbit data points
function [x,y] = circularOrbit(mu,r,pts)
init    
    T = PeriodfromA(mu,r);
    t = linspace(0,T,pts);
    n = MMfromP(T);
    x = r * sin(n.*t);
    y = r * cos(n.*t);
end

% compute gravity-gradient matrix G
function G = computeGravityGradient(R, mu)
    r = norm(R);
    G = (3 * mu / r^5) * (R' * R - (r^2) * eye(3));
end
