clear;
close all;
clc;
init;
% Define constants

mu = Gearth;    % Gravitational parameter (m^3/s^2)
Re = Rearth;    % Earth's radius (m)
alt_init = 200;             % initial altitude (km)
alt_ref = 35786;            % reference altitude (km)

%% Define reference Orbit
a_ref = alt_ref + Rearth;   % reference semi-major axis (km)
e_ref = 0;                  % eccentricity of the reference orbit (circular)

% Reference state vector (a, e, i, Omega, omega, nu)
ref_state = [a_ref; e_ref; 0; 0; 0; 0]; % Assuming circular orbit for simplicity

%% Define initial  
a_init = alt_init + Rearth; % initial semi-major axis
e_init = 0;                 % eccentricity of initial orbit (circular)
% initial state vector (a, e, i, Omega, omega, nu)
initial_state = [a_init; e_init; 0; 0; 0; 0]; % Initial state

% Define the perturbing accelerations (example values)
a_r = 0; % Radial component
a_theta = 0; % Transverse component
a_h = 0; % Normal component


%% Define the gain matrix (example values)
K = diag([1, 1, 1, 1, 1, 1]); % Gain matrix

% Numerical integration
tspan = [0 86400]; % Time span (1 day)
state = initial_state; % Initial orbital elements

% Integrate the variational equations with control
[t, state] = ode45(@(t, state) gauss_variational_equations_with_control(t, state, ref_state, K, mu), tspan, initial_state);

%% functions

% Control law
function u = control_law(state, ref_state, K)
    error = ref_state - state;
    u = K * error;
end


% Gauss' variational equations
function dxdt = gauss_variational_equations(t, state)
    mu = Gearth;
    a = state(1);
    e = state(2);
    i = state(3);
    Omega = state(4);
    omega = state(5);
    nu = state(6);
    
    % Perturbing accelerations
    a_r = 0; % Example value
    a_theta = 2.5/a; % Example value
    a_h = 0; % Example value
    
    % Specific angular momentum
    h = sqrt(mu * a * (1 - e^2));
    
    % Gauss' variational equations
    da_dt = (2 * a^2 / h) * (e * sin(nu) * a_r + (1 + e * cos(nu)) * a_theta);
    de_dt = (1 / h) * (a * sin(nu) * a_r + ((a + Re) * cos(nu) + Re * e) * a_theta);
    di_dt = (Re * cos(nu + omega) / h) * a_h;
    dOmega_dt = (Re * sin(nu + omega) / (h * sin(i))) * a_h;
    domega_dt = (1 / (h * e)) * (-a * sin(nu) * a_r + ((a + Re) * cos(nu) + Re * e) * a_theta) - (Re * sin(nu + omega) * cos(i) / (h * sin(i))) * a_h;
    dnu_dt = h / Re^2 + (1 / (h * e)) * (a * sin(nu) * a_r - ((a + Re) * cos(nu) + Re * e) * a_theta);
    
    dxdt = [da_dt; de_dt; di_dt; dOmega_dt; domega_dt; dnu_dt];
end




% Gauss' variational equations with control
function dxdt = gauss_variational_equations_with_control(t, state, ref_state, K, mu)
    % Compute the control input
    Re = 6378.14;

    u = control_law(state, ref_state, K);
    
    % Perturbing accelerations (example values)
    a_r = u(1); % Radial component
    a_theta = 2.5/state(1); % Transverse component
    a_h = u(3); % Normal component
    
    % Specific angular momentum
    h = sqrt(mu * state(1) * (1 - state(2)^2));
    
    % Gauss' variational equations
    da_dt = (2 * state(1)^2 / h) * (state(2) * sin(state(6)) * a_r + (1 + state(2) * cos(state(6))) * a_theta);
    de_dt = (1 / h) * (state(1) * sin(state(6)) * a_r + ((state(1) + Re) * cos(state(6)) + Re * state(2)) * a_theta);
    di_dt = (Re * cos(state(6) + state(5)) / h) * a_h;
    dOmega_dt = (Re * sin(state(6) + state(5)) / (h * sin(state(3)))) * a_h;
    domega_dt = (1 / (h * state(2))) * (-state(1) * sin(state(6)) * a_r + ((state(1) + Re) * cos(state(6)) + Re * state(2)) * a_theta) - (Re * sin(state(6) + state(5)) * cos(state(3)) / (h * sin(state(3)))) * a_h;
    dnu_dt = h / Re^2 + (1 / (h * state(2))) * (state(1) * sin(state(6)) * a_r - ((state(1) + Re) * cos(state(6)) + Re * state(2)) * a_theta);
    
    dxdt = [da_dt; de_dt; di_dt; dOmega_dt; domega_dt; dnu_dt];
end


