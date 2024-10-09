% ASEN-6015
% Hwk 1 
% Problem 2 Part a
close all
clear
clc
init;
dt = 10;
mu = Gearth;
g = 0.00908665;        % accel due to gravity

r0_hat = [1 0 0];
r0 = 6578;
R0 = r0* r0_hat';
V0 = VfromRA(mu,R0,R0)';
rf_hat = [-cosd(5); sind(5); 0];
rf = 42164;
Rf = rf * rf_hat;  % final position

Tof = 37864;
tspan = 0:dt:Tof;

[v1, v2,iter] = lambert_univ(mu,R0,Rf,Tof,1);
iter
% Problem 2 Part b
% simulate a case with simple lambert guidance routine with cross product
% guidance

% iterate through the simulation, constantly updating r0 until we get to rf
ai = 0.03;      % constant acceleration 30 m/s
r = R0;
v = v1;
tspan = 0:dt:Tof; 
tof = Tof; 
burn_a = 0;
for i = 1:length(tspan)
    vg = v2 - v1;       % velocity to go
    r_mag = norm(r);     
    % Calculate gravitational force
    F_gravity = 0; % -g * r/r_mag (ignore for now)

    % Apply transfer burns          
    if norm(vg) > 0.01 
       a_burn = -ai*vg/norm(vg);     % in direction of vg
       burn_a(i) = 1;                % track the burn time
    else
       a_burn = [0; 0; 0];
       burn_a(i) = 0;
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
burn_a(i) 
[x,y] = circularOrbit(mu,r0,length(tspan));

% Problem 2 Part 3
[v1, v2,~] = lambert_univ(mu,R0,Rf,Tof,1);

G = computeGravityGradient(Rf, mu); % Compute gravity-gradient matrix

% reset variables
r = R0;      
v = v1;
tof = Tof;
burn_b = 0;
dvg = (v1-V0)/dt;
gamma = 0.9;
% iterate through the simulation, constantly updating r0 until we get to rf
for i = 1:length(tspan)
    vg = v2 - v1;       % velocity to go 
    vghat = vg/norm(vg);
    r_mag = norm(r);    
    % Calculate gravitational force
    F_gravity = 0; % -g * r/r_mag (ignore for now)

    % Apply transfer burns
    % when we get close to the target delta-v then no more burns until
    % error grows too big OR get sufficiently close to target 
    if norm(vg) > 0.01
        Q = -(1/dt) * (eye(3) + (dt^2 * G / 2)); % Compute Q matrix
       % Q = -1/(dt*i - tof) * (eye(3) + (dt*i - tof)^2 * G / 2); % Compute Q matrix
        b = -Q *vg;
        C = norm(b - dvg);
        q = -dot(vghat,b)-sqrt(dot(vghat,b)^2 - norm(b)^2 + C^2);
        a_burn = ai*(gamma*b + (q-dot(gamma*b,vghat))*vghat);  % Update acceleration vector using Q matrix
        burn_b(i) = 1;             % count the burn time
    else
       a_burn = [0; 0; 0];
       burn_b(i) = 0;  
    end

    % Update velocity and position
    a_total = F_gravity + a_burn;

    v1 = v1 + a_total * dt;
    r = r + v1 * dt;
    tof = tof - dt;                 % decrement time remaining
    [v1, v2] = lambert_univ(mu,r,Rf,tof,1);
    dvg = (vg-(v2-v1))/dt;
    % Store position and velocity
    R1(:, i) = r;
    V1(:, i) = v1;
end
burn_b(i)
delR = R(:,end) -Rf
delR1 = R1(:,end) - Rf

figure(3);
plot(R(1, :), R(2, :), R1(1,:),R1(2,:));
xlabel('X (m)');
ylabel('Y (m)');
title('Orbits Solving Lamberts Problem') ;
grid on;
hold on;
plot(x,y)
legend('Lambert guidance','Combination steering', 'initial orbit','Location','southwest')
hold off;
% xlim([-60000 40000]);
% ylim([-60000 40000]);

% plot from nominal

figure (4)
plot(tspan(1:end-1),vecnorm(V1(:,1:end-1))-vecnorm(V(:,1:end-1)))
title('Difference in \Deltav betweenn Lambert with Q-guidance')

figure(5)
plot(tspan(1:end-1),vecnorm(V(:,1:end-1)),tspan(1:end-1),vecnorm(V1(:,1:end-1)))
title('vector')

figure(6)
plot(tspan(end-100:end),burn_a(end-100:end),tspan(end-100:end),burn_b(end-100:end))
title('Times of Burn')
legend('Lambert guidance','Cross-product guidance','Location','southwest')
ylim([0,2]);

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
x = R(1);
y = R(2);
z = R(3);
r = norm(R);

G = [mu/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - (3*mu*x*abs(x)*(x + conj(x)))/(2*(x*conj(x))^(1/2)*(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2)),   -(3*mu*x*abs(y)*(y + conj(y)))/(2*(y*conj(y))^(1/2)*(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2)),  0;
-(3*mu*y*abs(x)*(x + conj(x)))/(2*(x*conj(x))^(1/2)*(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2)), mu/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - (3*mu*y*abs(y)*(y + conj(y)))/(2*(y*conj(y))^(1/2)*(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2)),                                            0;
-(3*mu*z*abs(x)*(x + conj(x)))/(2*(x*conj(x))^(1/2)*(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2)),  -(3*mu*z*abs(y)*(y + conj(y)))/(2*(y*conj(y))^(1/2)*(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2)), mu/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2)] ;
end
