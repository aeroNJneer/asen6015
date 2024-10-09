clear
close all

% Constants
VM = 3000 * 0.3048; % Missile velocity (ft/s to m/s)
VT = 1000 * 0.3048; % Target velocity (ft/s to m/s)
R0 = 40000 * 0.3048; % Initial range (ft to m)
tgo = R0 / (VM + VT); % Time to go (s)
dt = 1; % Time step (s)
t_final = tgo; % Simulation time (s)

% Initial conditions
x0 = [0; VM * sin(deg2rad(30)); 0]; % Initial state vector [y, y_dot, nT]
theta_m = deg2rad(30); % Missile initial heading (30 degrees in radians)
theta_t = 0; % Target initial heading (rad)

% Time vector
t = 0:dt:t_final;

% Preallocate arrays for results
x = zeros(3, length(t));
lambda = zeros(3, length(t));
nc = zeros(1, length(t));
n1 = zeros(1, length(t));

% Initial state
x(:, 1) = x0;

% State transition matrix
A = [0 1 0; 0 0 1; 0 0 0];
B = [0; -1; 0];

% Performance index (final miss distance)
J = @(x) x(1, end)^2;

% Simulation loop
for i = 2:length(t)
    % Update target position
    x_t = R0 - VT * t(i);
    
    % Calculate line-of-sight (LOS) angle
    LOS_angle = atan2(0 - x(1, i-1), x_t - x(1, i-1));
    
    % Calculate closing velocity
    Vc = VM * cos(theta_m - LOS_angle) + VT * cos(theta_t - LOS_angle);
    
    % Calculate heading error
    heading_error = LOS_angle - theta_m;
    
    % Proportional Navigation (PN)
    nc(i) = 3 * Vc * heading_error; % Using N = 3 for simplicity
    
    % Update state
    x(:, i) = x(:, i-1) + (A * x(:, i-1) + B * nc(i)) * dt;
end

% Solve adjoint equations
for i = length(t):-1:2
    % Adjoint equations
    lambda(:, i-1) = lambda(:, i) - (A' * lambda(:, i)) * dt;
end

% Update control using adjoint variables
for i = 1:length(t)
    nc1(i) = -lambda(2, i) / B(2);
end

% Plot results
figure;
subplot(2, 1, 1);
plot(t, x(1, :), 'DisplayName', 'y');
xlabel('Time (s)');
ylabel('y (m)');
title('Miss Distance vs Time');
legend;
grid on;

subplot(2, 1, 2);
plot(t, nc, t,nc1,'DisplayName', 'n_c');
xlabel('Time (s)');
ylabel('n_c (m/s^2)');
title('Control Command vs Time');
legend;
grid on;
