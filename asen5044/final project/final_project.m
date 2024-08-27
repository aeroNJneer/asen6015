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
omega = sqrt(mu/r0^3); % angular velocity of satellite
% Atilde is time varying so need to create paramaterized variable.... 
% -- pass in r and t (r0 for nominal)
% Atilde = [0 1 0 0; 
%     mu*(2*(cos(omega*t))^2-(sin(omega*t))^2)/r^3 0 3*mu*sin(omega*t)*cos(omega*t)/r^3 0;
%     0 0 0 1;
%     3*mu*sin(omega*t)*cos(omega*t)/r^3 0 mu*(2*(sin(omega*t)^2 -(cos(omega*t))^2))/r^3 0];
A0 = [0 1 0 0 ; 2*mu/r0^3 0 0 0 ; 0 0 0 1; 0 0 -mu/r0^3 0];
B = [0 0; 1 0; 0 0; 0 1];
F = eye(4) + Atilde*dt; % pass in r and t to Atilde... 
G = dt * B; 

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

% % linearized nominal simulation
% Define nominal state and control input over time
dx = zeros(4, length(t)); % Perturbation state vector
theta = omega * t;
xnom = r0*cos(omega * t);
ynom = r0*sin(omega * t);
xdotnom = -r0 * omega * sin(omega * t);
ydotnom = r0 * omega * cos(omega * t);

u = zeros(2, length(t));
% Initialize perturbation state
dx(:, 1) = [0, 0.075, 0, -0.021]';
% Compute perturbation states over time
for i = 2:length(t)
    dx(:, i) = F*dx(:, i - 1) + G*u(:, i); % need to pass in r and t to Atilde...
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


