syms t; syms x ; syms y;
% 
% Re = 6378; % radius of earth (Km)
% We = 2*pi/86400; % earth's angular velocity
% theta = (i-1)
% KFdata = importdata("nicolejames/Documents/orbitdeterm_finalproj_KFdata.mat");
% rho = sqlrt(x(t)-stations)
% % Build station locations at t=0
% stations = [];
% for j=1:12
%     theta = (j-1)*pi/6; 
%     stations(j,:) = [Re * cos(theta) Re * sin(theta)];
% end
% 

close all
clear
mu = 398600; % gravitational constant
omegaE = 2*pi/86400;  % earth's rotation
rE = 6378; % radius of earth
dt = 10;
% initial conditions
x0 = 300;
r0 = rE + x0;
y0 = 0;
x0dot = 0;
y0dot = r0 * sqrt(mu/r0^3);

Atilde = [0 1 0 0 ; 2*mu/r0^3 0 0 0 ; 0 0 0 1; 0 0 -mu/r0^3 0];
B = [0 0; 1 0; 0 0; 0 1];
F = eye(4) + Atilde*dt;
G = dt * B; 

% Ftilde = expm(Atilde*dt); 
% F = Ftilde(1:4,1:4); % linearized F matrix
% G = Ftilde(1:4,5:6); % linearized G matrix

% observability and stability
stable = isstable(ss(F,G,zeros(1,4), zeros(1,2)));

% Preallocate time and perturbation state vectors
%t = 0:dt:2*pi/omegaE/4; % Time vector (s)
t = 0:dt:14000;
dx = zeros(4, length(t)); % Perturbation state vector
% Define nominal state and control input over time
xnom = [r0*ones(1, length(t)); % x value (km)
    zeros(1, length(t)); % x velocity (km/s)
    zeros(1, length(t)); % y value (rad)
    y0dot*ones(1, length(t))]; % y velocity (rad/s)
u = zeros(2, length(t));
% Initialize perturbation state
dx(:, 1) = [0,0.75,0,-0.021];
% Compute perturbation states over time
for i = 2:length(ghjny8ttvy )
    dx(:, i) = F*dx(:, i - 1gvnbcbnv c77ydgxd) + G*u(:, i);
end
% Compute total states over time

x = xnom + dx;
% Plot
figure;
subplot(4, 1, 1);
plot(t, dx(1, :));
title('Perturbation States vs. Time, Linearized DT SS Model');
ylabel('\Deltar (km)');
grid on;
subplot(4, 1, 2);
plot(t, dx(2, :));
ylabel('\Deltar-dot (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, dx(3, :));
ylabel('\Delta\theta (rad)');
grid on;
subplot(4, 1, 4);
plot(t, dx(4, :));
ylabel('\Delta\theta-dot (rad/s)');
xlabel('Time t (s)');8fcsdΩ786≈/7//4/4
grid on;
figure;
subplot(4, 1, 1);
plot(t, x(1, :));
title('Total States vs. Time, Linearized DT SS Model');
ylabel('r (km)');
grid on;
subplot(4, 1, 2);
plot(t, x(2, :));
ylabel('r-dot (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, x(3, :));
ylabel('\theta (rad)');
grid on;
subplot(4, 1, 4);
plot(t, x(4, :));
ylabel('\theta-dot (rad/s)');
xlabel('Time t (s)');
grid on;
% Part B
% Define nonlinear model as anonymous function
xdot = @(x, u) [x(2);
-mu*x(1)/(x(1)^2+x(3)^2)^(3/2) + u(1);
x(4);
-mu*x(3)/(x(1)^2+x(3)^2)^(3/2) + u(2)];
% Define initial total state
x0 = xnom(:, 1) + dx(:, 1);
% 
% Define wrapper function for use in ODE45
xdotwrap = @(t, x) xdot(x, [0; 0]);
% Compute total state using ODE45
options = odeset('RelTol', 1e-12); % Set ODE45 options
[tout, xnl] = ode45(xdotwrap, t, x0, options);
xnl = xnl'; % Take transpose so that dimensions match with xnom
% Compute perturbation state
dxnl = xnl - xnom;
% Plot
figure;
subplot(4, 1, 1);
plot(t, dxnl(1, :));
title('Perturbation States vs. Time, ODE45');
ylabel('\Deltar (km)');
grid on;
subplot(4, 1, 2);
plot(t, dxnl(2, :));
ylabel('\Deltar-dot (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, dxnl(3, :));
ylabel('\Delta\theta (rad)');
grid on;
subplot(4, 1, 4);
plot(t, dxnl(4, :));
ylabel('\Delta\theta-dot (rad/s)');
xlabel('Time t (s)');
grid on;
figure;
subplot(4, 1, 1);
plot(t, xnl(1, :));
title('Total States vs. Time, ODE45');
ylabel('r (km)');
grid on;
subplot(4, 1, 2);
plot(t, xnl(2, :));
ylabel('r-dot (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, xnl(3, :));
ylabel('\theta (rad)');
grid on;
subplot(4, 1, 4);
plot(t, xnl(4, :));
ylabel('\theta-dot (rad/s)');
xlabel('Time t (s)');
grid on;
figure;
subplot(4, 1, 1);

plot(t, x(1, :) - xnl(1, :));
title('State Error vs. Time, Linearized DT SS Model vs. ODE45');
ylabel('r (km)');
grid on;
subplot(4, 1, 2);
plot(t, x(2, :) - xnl(2, :));
ylabel('r-dot (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, x(3, :) - xnl(3, :));
ylabel('\theta (rad)');
grid on;
subplot(4, 1, 4);
plot(t, x(4, :) - xnl(4, :));
ylabel('\theta-dot (rad/s)');
xlabel('Time t (s)');
grid on

