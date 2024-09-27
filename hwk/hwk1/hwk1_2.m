clear
clc
close all
init;
% Define constants
mu = Gearth;
Re = Rearth;
t_final = 3600*4;

% state parameters
a_initial = 200 + Re;
a_ref = 1000 + Re ;
e_initial = 0.1 ;
e_ref = 0.11;
i_initial = 0.22;
i_ref = 0.25;
Omega_initial = 0.5;
Omega_ref = 0.4;
omega_initial = 0.22;
omega_ref = 0.23;
M_initial = 0 ;

mm_ref = MMfromMua(mu,a_ref);

M_ref = mm_ref*t_final;

% Define the gain matrix (example values)
K = diag([1, 500, 1, 1, 1, 1]); % Gain matrix

% Define the reference state vector (reference COEs)
xref = [a_ref; e_ref; i_ref; Omega_ref; omega_ref; M_ref];

% Numerical integration
tspan = linspace(0,t_final); % Time span 
initial_state = [a_initial; e_initial; i_initial; Omega_initial; omega_initial; M_initial]; % Initial orbital elements

% Integrate the variational equations with control
[t, state] = ode45(@(t, state) gauss_variational_equations_with_control(t, state, xref, K), tspan, initial_state);

% Plot the results
figure;
s  = stackedplot(t,state);
xlabel('Time t (s)');
s.DisplayLabels={'Semi-major axis (km)', 'Eccentricity', 'Inclination (rad)', 'RAAN (rad)','Argument of Perigee (rad)','Mean Anomaly (rad)'};
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Orbital Element State Change with Linear Control');

%% functions 
% Control law with linear feedback
function u = control_law(state, xref, K)
    error = state -xref;
    u = -K * error;
end

% Gauss' variational equations with control
function dxdt = gauss_variational_equations_with_control(t, state, xref, K)
    mu = evalin('base', 'mu');
    % Compute the control input
    u = control_law(state, xref, K);

    [r,v] = perifocal_sv_from_coe(state);
    [x y  z] = eci2rdl(r,v);
    A = [x ;y ;z]' * v';
    % Perturbing accelerations (example values)
    a_r = A(1); % Radial component
    a_theta = A(2)/norm(r); % Transverse component
    a_h = A(3); % Normal component
    

    % Specific angular momentum
    h = sqrt(mu * state(1) * (1 - state(2)^2));
    
    % Gauss' variational equations
    da_dt = (2 * state(1)^2 / h) * (state(2) * sin(state(6)) * a_r + (1 + state(2) * cos(state(6))) * a_theta);
    de_dt = (1 / h) * (sin(state(6)) * a_r + (cos(state(6)) + state(2)) * a_theta);
    di_dt = (state(1) * cos(state(6) + state(5)) / h) * a_h;
    dOmega_dt = (state(1) * sin(state(6) + state(5)) / (h * sin(state(3)))) * a_h;
    domega_dt = (1 / (h * state(2))) * (-cos(state(6)) * a_r + (2 + state(2) * cos(state(6))) * a_theta) - (state(1) * sin(state(6) + state(5)) * cos(state(3)) / (h * sin(state(3)))) * a_h;
    dM_dt = sqrt(mu / state(1)^3) - (2 / h) * (state(1) / (1 - state(2)^2)) * (state(2) * sin(state(6)) * a_r + (1 + state(2) * cos(state(6))) * a_theta);
    
    % Return the derivatives
    dxdt = [da_dt; de_dt; di_dt; dOmega_dt; domega_dt; dM_dt];
end

