% LEO to GEO transfer using Kluever's simple guidance scheme
clear; clc;
close all;
%% Constants
mu = 3.986e14; % Earth's gravitational parameter (m^3/s^2)
Re = 6371e3;   % Earth's radius (m)
g0 = 9.81;     % Gravity at Earth's surface (m/s^2)

% Phase 1 and 2 guidance parameters
guidance_params.phase1.Ke0 = 0;
guidance_params.phase1.Ke1 = -1.5e3;
guidance_params.phase1.Ke2 = 0.0;
guidance_params.phase1.Ki0 = -0.33;
guidance_params.phase1.Ki1 = -4.3e-3;

guidance_params.phase2.Ke0 = -0.18;
guidance_params.phase2.Ke1 = -1.0e-2;
guidance_params.phase2.Ke2 = 7e-5;
guidance_params.phase2.Ki0 = -0.85;
guidance_params.phase2.Ki1 = -5.2e-3;

% Motor characteristics
efficiency = 0.65;       % Efficiency of the thruster
power_input = 10e3;      % Power input to the thruster (W)
Isp = 3300;              % Specific impulse (s)
initial_TW_ratio = 3.4e-5;

% Initial spacecraft mass
m0 = 1200;               % Spacecraft initial mass (kg)
T = initial_TW_ratio * m0 * g0; % Thrust (N)

% Orbital parameters
alt_LEO = 550e3;          % Initial altitude (m)
r_LEO = Re + alt_LEO;     % Initial orbit radius (m)
v_LEO = sqrt(mu / r_LEO); % Initial orbital velocity (m/s)
inclination_LEO = deg2rad(28.5); % Initial inclination (radians)
RAAN_LEO = deg2rad(0);    % Right Ascension of Ascending Node (radians)
a_LEO = -mu/(v_LEO^2 - (2*mu/r_LEO));
e_LEO = 0;

% Target orbit parameters
r_GEO = 42164e3; % Geostationary orbit radius (m)
e_target = 0;    % Target eccentricity
i_target = 0;    % Target inclination

% Initial state in ECI frame
r0 = [r_LEO; 0; 0]; % Initial position in ECI (m)
v0 = [0; v_LEO; 0]; % Initial velocity in ECI (m/s)

% Rotation to account for inclination and RAAN
rot_matrix = [cos(RAAN_LEO), -sin(RAAN_LEO), 0;
              sin(RAAN_LEO), cos(RAAN_LEO), 0;
              0, 0, 1] * ...
             [1, 0, 0;
              0, cos(inclination_LEO), -sin(inclination_LEO);
              0, sin(inclination_LEO), cos(inclination_LEO)];

r0 = rot_matrix * r0;
v0 = rot_matrix * v0;

%% Time settings
dt = 60; % Time step (s)
t_max = 200 * 86400; % Maximum simulation time (s)

%% Preallocation for variables
r_ECI = r0; % Position in ECI
v_ECI = v0; % Velocity in ECI
m = m0;     % Mass
time = 0;

% Storage for plotting
trajectory = [];
mass_history = [];
time_history = [];
a_history = [];
e_history = [];
i_history = [];
yaw_history = [];
pitch_history = [];

% initialize variables
a = a_LEO;
e = e_LEO;
inc = inclination_LEO;

%% Simulation Loop
while norm(r_ECI) < r_GEO && time < t_max
    % Compute RTN frame unit vectors
    r_norm = norm(r_ECI);
    r_hat = r_ECI / r_norm;
    v_norm = norm(v_ECI);
    v_hat = v_ECI / v_norm;
    h_vec = cross(r_ECI, v_ECI);
    h_norm = norm(h_vec);
    h_hat = h_vec / norm(h_vec);
    t_hat = cross(h_hat, r_hat);
    
    % Orbital parameters
    true_anomaly = atan2(dot(r_ECI, t_hat), dot(r_ECI, r_hat));
    argument_of_periapse = 0; % Assume for simplicity; refine as needed

    % Guidance law: thrust direction in RTN frame
    % Compute flight path angle
    gamma = atan2(dot(v_ECI, r_hat), dot(v_ECI, t_hat));

    % In-plane steering angles    
    c_a = [sin(gamma); cos(gamma); 0];  % for SMA
    % optimal in-plane thrust for eccentricity
    phi_e = atan2(r_norm * sin(true_anomaly) , (2 * a * (e + cos(true_anomaly))));
    c_e = [sin(gamma + phi_e); 
           cos(gamma + phi_e); 
           0];

    % Determine current phase and guidance parameters
    if time <= 120 * 24 * 3600
        phase_params = guidance_params.phase1;
    else
        phase_params = guidance_params.phase2;
    end
            
    % Guidance parameter schedule with
    % Quadratic weighting function for eccentricity
    Ge = phase_params.Ke0 + phase_params.Ke1 * time + phase_params.Ke2 * time^2;
    
    % Linear weighting for out-of-plane steering
    Gi = phase_params.Ki0 + phase_params.Ki1 * time ;

    % Normalize and blend control vectors
    c_blended = Ge * c_e + c_a;
    c_blended = c_blended / norm(c_blended);
    
    % Extract pitch angle
    pitch_angle = atan2(c_blended(1), c_blended(2));
    
    % Out-of-plane steering / yaw
    beta = Gi * (pi/2) * cos(true_anomaly + argument_of_periapse);  % yaw angle

    % Compute thrust direction in RTN frame with yaw
    thrust_unit_RTN = [
        sin(pitch_angle) * cos(beta);
        cos(pitch_angle) * cos(beta);
        sin(beta)
    ];

    % alpha_a = 0; % Align thrust with velocity for max da/dt
    % alpha_e = atan2(r_norm * sin(gamma), 2 * r_LEO * (e_target + cos(gamma)));
    % alpha_blended = (1 - e_target) * alpha_a + e_target * alpha_e;
    % 
    % 
    % 
    % % Thrust direction in RTN
    % thrust_RTN = [sin(gamma + alpha_blended) * cos(beta); 
    %               cos(gamma + alpha_blended) * cos(beta); 
    %               sin(beta)];
    % 

    % Rotate thrust direction to ECI
    RTN_to_ECI = [r_hat, t_hat, h_hat];
    thrust_ECI = RTN_to_ECI * thrust_unit_RTN;

    % Update dynamics
    a_thrust = T / m * thrust_ECI; % Acceleration from thrust (m/s^2)
    a_grav = -mu / r_norm^3 * r_ECI; % Gravitational acceleration (m/s^2)
    a_total = a_thrust + a_grav;

    v_ECI = v_ECI + a_total * dt; % Update velocity (m/s)
    r_ECI = r_ECI + v_ECI * dt;   % Update position (m)

    % Mass depletion
    dm = T / (Isp * g0) * dt;
    m = m - dm;
    
    % update orbital elements
    a = a + 2 * a^2 * v_norm * (norm(a_thrust)) * dt / mu;
    e = e + (norm(a_thrust)/v_norm) * (2*(e + cos(true_anomaly)* cos(phi_e)) + (r_norm/a)*sin(true_anomaly)*sin(phi_e)) * dt;
    inc = inc + (norm(a_thrust)) * (r_norm/h_norm) * cos(true_anomaly + argument_of_periapse) * sin(beta) *dt;
    
    % Record data
    trajectory = [trajectory, r_ECI];
    mass_history = [mass_history, m];
    time_history = [time_history, time];
    a_history = [a_history, a];
    e_history = [e_history, e];
    i_history = [i_history, inc];
    yaw_history = [yaw_history, beta];
    pitch_history = [pitch_history, pitch_angle];


    % Update time
    time = time + dt;
end

%% Plot Results
figure;
plot3(trajectory(1, :) / 1e3, trajectory(2, :) / 1e3, trajectory(3, :) / 1e3);
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('LEO to GEO Transfer');
grid on;

figure;
plot(time_history / 3600, mass_history);
xlabel('Time (hours)');
ylabel('Mass (kg)');
title('Spacecraft Mass vs Time');
grid on;
