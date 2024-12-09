
% constants
constants.mu = 3.986e14;         % Earth's gravitational parameter (m^3/s^2)
constants.g0 = 9.81;             % Gravitational acceleration (m/s^2)
constants.tgo = 220;                % time to go (simulation time in days)
constants.dt = 600;           % time between iterations 
constants.tmax = 250;           % max time of simulation

%  WATT = kg m^2/s^3  Isp = sec  
% Simulation parameters
thrust_params.power = 10e3;           % 10 kW
thrust_params.Isp = 3300;             % Specific impulse
thrust_params.efficiency = 0.65;      % Engine efficiency
thrust_params.initial_mass = 1200;    % kg
thrust_params.initial_tw = 3.4e-5;    % thrust in mN (?)          
thrust_params.thrust_acceleration = thrust_params.initial_tw * constants.g0 * thrust_params.initial_mass;
thrust_params.mass_flow_rate = 2*thrust_params.efficiency* thrust_params.power/(constants.g0*thrust_params.Isp)^2; 

% Phase 1 and 2 guidance parameters  
guidance_params.phase1.Ke0 = 0;
guidance_params.phase1.Ke1 = -1.5e-3;
guidance_params.phase1.Ke2 = 0.0;
guidance_params.phase1.Ki0 = -0.33;
guidance_params.phase1.Ki1 = -4.3e-3;

guidance_params.phase2.Ke0 = -0.18;
guidance_params.phase2.Ke1 = -1.0e-2;
guidance_params.phase2.Ke2 = 7e-5;
guidance_params.phase2.Ki0 = -0.85;
guidance_params.phase2.Ki1 = -5.2e-3;

% INITIAL CONDITIONS
initial_altitude = 550e3;  
initial_inclination = deg2rad(28.5);
r_earth = 6378e3;  % Earth radius

% state vectors in orbit frame
r_initial = [r_earth + initial_altitude; 0; 0 ];
v_initial = [0; sqrt(constants.mu / norm(r_initial)); 0];

% Initial state in ECI - assume simulation starts at RAAN
rtn2eci = Euler313(0, initial_inclination, 0, 0);
r0 = rtn2eci * r_initial;
v0 = rtn2eci * v_initial;

% Target GEO conditions
target_semimajor_axis = 42164e3;
vf_mag = sqrt(constants.mu/target_semimajor_axis);

% initialize variables
r = r0;     % in ECI
v = v0;     % in ECI
t = 0;
m = thrust_params.initial_mass;
i = 1;

a_calc = norm(r0);
e_calc = 0;
inc_calc = initial_inclination;
results.state = [];
results.COEs = [];
results.time = [];

while (norm(r) < target_semimajor_axis && t < constants.tmax*24*3600)
    % get orbital and flight parameters
    COEs = eci2orb1(constants.mu, r,v);
    theta = COEs.omega + COEs.true_anomaly;
    gamma = getFlightPathAngle(r,v);
    period = 2*pi*sqrt(COEs.a^3/constants.mu);
    dt = period/100;

    % calculate optimal angles
    [phi_a, phi_e, beta] = getOptimalAngles(r, COEs.a, COEs.e, COEs.omega, COEs.true_anomaly);

    % Determine current phase and guidance parameters
    if t <= 120 * 24 * 3600
        phase_params = guidance_params.phase1;
    else
        phase_params = guidance_params.phase2;
    end

    % compute weighting functions
    G_a = 1;
    [G_e,G_i] = computeWeightingFunctions(t, phase_params);

    [c_a, c_e, c] = getBlendedCvector(G_a, G_e, gamma, phi_a, phi_e);

    % determine pitch and yaw steering angles
    [alpha, beta] = getPitchYaw(c, G_i, theta);

    thrust_vector_rtn = getThrustVector(alpha,beta);
    accel_thrust = thrust_params.thrust_acceleration / m;
    % variational orbital elements
    da_dt = 2*a_calc^2*norm(v)*accel_thrust*cos(phi_a)/constants.mu;
    a_calc = a_calc + da_dt*dt;
    de_dt = (accel_thrust/norm(v)) *(2*(e_calc+cos(COEs.true_anomaly))*cos(phi_e) + (r/a_calc)*sin(COEs.true_anomaly)*sin(phi_e));
    e_calc = e_calc + de_dt*dt;
    di_dt = accel_thrust*cos(theta)*sin(beta)/COEs.h;
    inc_calc = inc_calc + di_dt*dt;
  
    % simulate each orbit
    dt = period/100;
    tspan = t:dt:period;
            
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    [t_sim, states] = ode45(@(t,y) systemDynamics(t, y, thrust_vector_rtn, constants.mu, COEs, thrust_params), tspan, [r; v; m], options);
   
    results.states(i) = [results.state; states(end,1:3)];
  %  results.COEs = [results.COEs, COEs'];
    results.time(i) = [results.time, t_sim(end)];
    r = states(end,1:3)';
    v = states(end,4:6)';
    m = states(end,7);
    t = t+dt;
    i = i+1;
end



function gamma = getFlightPathAngle(r,v)
    % Compute the magnitudes of the position and velocity vectors
    r_mag = norm(r);
    v_mag = norm(v);

    % Compute the radial velocity component
    vr = dot(r, v) / r_mag;

    % Compute the flight path angle in radians
    gamma = acos(vr / v_mag);
end

function [phi_a, phi_e, beta] = getOptimalAngles(r, a, e, omega, nu)
    phi_a = 0;
    phi_e = atan2(norm(r)*sin(nu),2*a*(e+cos(nu)));
    beta = pi*cos(omega + nu)/2;
end

function [Ge, Gi] = computeWeightingFunctions(t, phase_params)
    % Quadratic weighting function for eccentricity
    Ge = phase_params.Ke0 + phase_params.Ke1 * t + phase_params.Ke2 * t^2;
    
    % Linear weighting for out-of-plane steering
    Gi = phase_params.Ki0 + phase_params.Ki1 * t;
end


function [c_a, c_e, c] = getBlendedCvector(G_a, G_e, gamma, phi_a, phi_e)
    % gamma = fpa, beta = yaw angle
    c_a = [sin(gamma + phi_a); cos(gamma + phi_a); 0];
    c_e = [sin(gamma + phi_e); cos(gamma +  phi_e); 0];
   
    c = (G_a*c_a + G_e*c_e)/(norm(G_a*c_a)+norm(G_e*c_e));
end 

function [alpha, beta] = getPitchYaw(c, G_i, theta)
    alpha = atan2(c(1),c(2));   
    beta = G_i*cos(theta);
end

function u_hat_rtn = getThrustVector(alpha, beta)
    u_hat_rtn = [sin(alpha)*cos(beta); cos(alpha)*cos(beta); sin(beta)];
end

function state_dot = systemDynamics(t, state, thrust_vector_RTN, mu, coes, thrust_params)
    % State vector: [x, y, z, vx, vy, vz, mass]
    r = state(1:3);
    v = state(4:6);
    m = state(7);

    % Gravitational acceleration
    r_mag = norm(r);
    a_grav = -mu * r / (r_mag^3);

    % Transform thrust direction to ECI frame
    R_RTN2ECI = compute_RTN_to_ECI_matrix(r, v);
    a_rtn2eci = Euler313(coes.raan, coes.inc, coes.omega, coes.true_anomaly);

    thrust_unit_ECI = R_RTN2ECI * thrust_vector_RTN;
    a_thrust = thrust_params.thrust_acceleration * thrust_unit_ECI / m;

    % Mass rate change
    mdot = -thrust_params.mass_flow_rate;
    
    % State derivative
    state_dot = [v; a_grav + a_thrust; mdot];    
end

% Compute RTN to ECI transformation matrix
function R_RTN2ECI = compute_RTN_to_ECI_matrix(r, v)
    % Radial unit vector (R direction)
    r_unit = r / norm(r);
    
    % Velocity unit vector for transverse direction
    v_unit = v / norm(v);
    
    % Normal unit vector (cross product of r and v)
    n_unit = cross(r, v) / norm(cross(r, v));
    
    % Transverse unit vector (completes right-handed system)
    t_unit = cross(n_unit, r_unit);
    
    % Transformation matrix from RTN to ECI
    R_RTN2ECI = [r_unit, t_unit, n_unit]';
end