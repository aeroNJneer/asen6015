close all
clear
mu = 398600; % gr62)avitational constant
omegaE = 2*pi/86400;  % angular velocity of earth
rE = 6378; % radius of earth
dt = 10;
t = 0:dt:14000;
% initial conditions
x0 = 6678;
y0 = 0;
r0 = sqrt(x0^2+y0^2);
x0dot = 0;
y0dot = r0 * sqrt(mu/r0^3);
period = 2 * pi * sqrt(r0^3/mu); % orbital period
omega = 2 * pi / period; % angular velocity of satellite
Atilde = [0 1 0 0 ; 2*mu/r0^3 0 0 0 ; 0 0 0 1; 0 0 -mu/r0^3 0];
B = [0 0; 1 0; 0 0; 0 1];
F = eye(4) + Atilde*dt;
G = dt * B; 

% Ftilde = expm(Atilde*dt); 
% F = Ftilde(1:4,1:4); % linearized F matrix
% G = Ftilde(1:4,5:6); % linearized G matrix

% observability and stability
stable = isstable(ss(F,G,zeros(1,4), zeros(1,2)));

% Define nonlinear model as anonymous function
pert = [0, 0.075, 0, -0.021]'; % perturbation 
xdot = @(x, u) [x(2);
    -mu*x(1)/sqrt((x(1)^2+x(3)^2)^3) + u(1);
    x(4);
    -mu*x(3)/sqrt((x(1)^2+x(3)^2)^3) + u(2)];
% Define initial total state
X0 = [x0, x0dot, y0, y0dot]' +  pert; 
% 
% Define wrapper function for use in ODE45
xdotwrap = @(t, x) xdot(x, [0; 0]);
% Compute total state using ODE45
options = odeset('RelTol', 1e-12); % Set ODE45 options
[tout, xnl] = ode45(xdotwrap, t, X0, options);
xnl = xnl'; % Take transpose so that dimensions match with xnom

% Plot
figure;
subplot(4, 1, 1);
plot(t, xnl(1, :));
title('Total States vs. Time, Nonlinear Dynanmics - ODE45');
ylabel('x (km)');
grid on;
subplot(4, 1, 2);
plot(t, xnl(2, :));
ylabel('\deltax (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, xnl(3, :));
ylabel('y (km)');
grid on;
subplot(4, 1, 4);
plot(t, xnl(4, :));
ylabel('\deltay (km/s)');
xlabel('Time t (s)');
grid on;
% end of NL simulation

% linearized nominal simulation
t = 0:10:140000;
dx = zeros(4, length(t)); % Perturbation state vector
theta = omega * t;
xnom = r0*cos(theta(t));
ynom = r0*sin(theta(t));
xdotnom = ;
ydotnom = ;




% Define nominal state and control input over time
xnom = [r0*ones(1, length(t)); % x value (km)
    zeros(1, length(t)); % x velocity (km/s)
    zeros(1, length(t)); % y value (rad)
    y0dot*ones(1, length(t))]; % y velocity (rad/s)
u = zeros(2, length(t));
% Initialize perturbation state
dx(:, 1) = [0, 0.075, 0, -0.021]';
% Compute perturbation states over time
for i = 2:length(t)
    dx(:, i) = F*dx(:, i - 1) + G*u(:, i);
end

% Compute total states over time
x = xnom + dx;

% Plot
figure;
subplot(4, 1, 1);
plot(t, x(1, :));
title('Total States vs. Time, Linearized DT SS Model');
ylabel('x (km)'); ylim([-7000, 7000])
grid on;
subplot(4, 1, 2);
plot(t, x(2, :));
ylabel('\deltax (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, x(3, :));
ylabel('y (km)');
grid on;
subplot(4, 1, 4);
plot(t, x(4, :));
ylabel('\deltay (km/s)');
xlabel('Time t (s)');
grid on;

figure;
subplot(4, 1, 1);
plot(t, x(1, :) - xnl(1, 1:length(t)));
title('State Error vs. Time, Linearized DT SS Model vs. ODE45');
ylabel('x (km)');
grid on;
subplot(4, 1, 2);
plot(t, x(2, :) - xnl(2, 1:length(t)));
ylabel('dx (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, x(3, :) - xnl(3, 1:length(t)));
ylabel('y (km)');
grid on;
subplot(4, 1, 4);
plot(t, x(4, :) - xnl(4, 1:length(t)));
ylabel('dy (km/s)');
xlabel('Time t (s)');
grid on


