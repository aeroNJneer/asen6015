%% MATLAB Script: Kluever's Simple Guidance Scheme for LEO to GEO Transfer
clear; clc;

%% Constants and Parameters
mu = 3.986e14; % Earth's gravitational parameter (m^3/s^2)
Re = 6371e3;   % Earth's radius (m)
g0 = 9.81;     % Gravity at Earth's surface (m/s^2)

% Thruster and spacecraft properties
efficiency = 0.65;       % Efficiency of the thruster
power_input = 10e3;      % Power input to the thruster (W)
Isp = 3300;              % Specific impulse (s)
initial_TW_ratio = 3.4e-5;
m0 = 1200;               % Initial mass (kg)
T = initial_TW_ratio * m0 * g0; % Thrust (N)

% Initial orbital elements
alt_LEO = 550e3;          % Initial altitude (m)
r_LEO = Re + alt_LEO;     % Initial orbit radius (m)
v_LEO = sqrt(mu / r_LEO); % Initial orbital velocity (m/s)
inclination_LEO = deg2rad(28.5); % Initial inclination (radians)
RAAN_LEO = deg2rad(0);    % Right Ascension of Ascending Node (radians)
omega_LEO = deg2rad(0);  % Argument of periapsis (radians)
true_anomaly_LEO = deg2rad(0); % True anomaly (radians)

% Target orbital elements
a_GEO = 42165e3; % Target semimajor axis (m)
e_target = 0;    % Target eccentricity
i_target = 0;    % Target inclination

% Initial state in ECI frame
r0 = [r_LEO; 0; 0];
v0 = [0; v_LEO; 0];

% Rotate initial state to account for inclination, RAAN, and argument of periapsis
rot_matrix = [cos(RAAN_LEO) * cos(omega_LEO) - sin(RAAN_LEO) * sin(omega_LEO) * cos(inclination_LEO), ...
              -cos(RAAN_LEO) * sin(omega_LEO) - sin(RAAN_LEO) * cos(omega_LEO) * cos(inclination_LEO), ...
              sin(RAAN_LEO) * sin(inclination_LEO);
              sin(RAAN_LEO) * cos(omega_LEO) + cos(RAAN_LEO) * sin(omega_LEO) * cos(inclination_LEO), ...
              -sin(RAAN_LEO) * sin(omega_LEO) + cos(RAAN_LEO) * cos(omega_LEO) * cos(inclination_LEO), ...
              -cos(RAAN_LEO) * sin(inclination_LEO);
              sin(omega_LEO) * sin(inclination_LEO), cos(omega_LEO) * sin(inclination_LEO), cos(inclination_LEO)];

r0 = rot_matrix * r0;
v0 = rot_matrix * v0;

%% Simulation Parameters
dt = 10; % Time step (s)
t_max = 220 * 86400; % Maximum simulation time (s)

%% Preallocate Storage for Results
r_ECI = r0;
v_ECI = v0;
m = m0;
time = 0;
trajectory = [];
mass_history = [];
time_history = [];
orbit_elements = [];
yaw_history = [];
pitch_history = [];

% Guidance parameters
Ki0 = 0.01; % Initial inclination control parameter
Ki1 = 1e-5; % Time-dependent inclination control parameter

%% Simulation Loop
while norm(r_ECI) < a_GEO && time < t_max
    % Compute RTN frame unit vectors
    r_norm = norm(r_ECI);
    r_hat = r_ECI / r_norm;
    qv_norm = norm(v_ECI);
    h_vec = cross(r_ECI, v_ECI);
    h_hat = h_vec / norm(h_vec);
    t_hat = cross(h_hat, r_hat);

    % Orbital parameters
    specific_angular_momentum = cross(r_ECI, v_ECI);
    h = norm(specific_angular_momentum);
    eccentricity_vector = cross(v_ECI, specific_angular_momentum) / mu - r_ECI / r_norm;
    e = norm(eccentricity_vector);
    inclination = acos(h_hat(3));
    true_anomaly = atan2(dot(r_ECI, t_hat), dot(r_ECI, r_hat));

    % Flight path angle
    gamma = atan2(dot(v_ECI, r_hat), dot(v_ECI, t_hat));

    % Guidance laws
    % In-plane steering
    alpha_a = 0; % For max da/dt
    alpha_e = atan2(r_norm * sin(gamma), 2 * a_GEO * (e_target + cos(gamma)));
    alpha_blended = (1 - e_target) * alpha_a + e_target * alpha_e;

    % Out-of-plane steering (yaw)
    Gi = Ki0 + Ki1 * time; % Time-dependent inclination control weighting
    beta = Gi * cos(true_anomaly + omega_LEO);

    % Thrust direction in RTN
    thrust_RTN = [sin(gamma + alpha_blended) * cos(beta); 
                  cos(gamma + alpha_blended) * cos(beta); 
                  sin(beta)];

    % Rotate thrust direction to ECI
    RTN_to_ECI = [r_hat, t_hat, h_hat];
    thrust_ECI = RTN_to_ECI * thrust_RTN;

    % Update dynamics
    a_thrust = T / m * thrust_ECI;
    a_grav = -mu / r_norm^3 * r_ECI;
    a_total = a_thrust + a_grav;

    v_ECI = v_ECI + a_total * dt;
    r_ECI = r_ECI + v_ECI * dt;

    % Mass depletion
    dm = T / (Isp * g0) * dt;
    m = m - dm;

    % Store results
    trajectory = [trajectory, r_ECI];
    mass_history = [mass_history, m];
    time_history = [time_history, time];
    yaw_history = [yaw_history, beta];
    pitch_history = [pitch_history, gamma + alpha_blended];
    orbit_elements = [orbit_elements; r_norm, e, inclination];

    % Update time
    time = time + dt;
end

%% Plot Results
figure;
plot3(trajectory(1, :) / 1e3, trajectory(2, :) / 1e3, trajectory(3, :) / 1e3);
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('LEO to GEO Transfer Trajectory');
grid on;

figure;
plot(time_history / 3600, mass_history);
xlabel('Time (hours)');
ylabel('Mass (kg)');
title('Spacecraft Mass vs Time');
grid on;

figure;
plot(time_history / 3600, yaw_history);
xlabel('Time (hours)');
ylabel('Yaw Angle (rad)');
title('Yaw Steering History');
grid on;

figure;
plot(time_history / 3600, pitch_history);
xlabel('Time (hours)');
ylabel('Pitch Angle (rad)');
title('Pitch Steering History');
grid on;

figure;
plot(time_history / 3600, orbit_elements(:, 3));
xlabel('Time (hours)');
ylabel('Inclination (rad)');
title('Inclination vs Time');
grid on;
