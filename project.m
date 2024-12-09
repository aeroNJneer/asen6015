
% constants
constants.mu = 3.986e14;         % Earth's gravitational parameter (m^3/s^2)
constants.g0 = 9.81;             % Gravitational acceleration (m/s^2)
constants.tgo = 216;                % time to go (simulation time)
constants.dt = 3600;           % time between iterations 

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
v_initial = [0; constants.mu / norm(r_initial); 0];

% Initial state in ECI - assume simulation starts at RAAN
rtn2eci = Euler313(0, initial_inclination, 0, 0);
r0 = rtn2eci * r_initial;
v0 = rtn2eci * v_initial;

% Target GEO conditions
target_semimajor_axis = 42164e3;
vf_mag = sqrt(constants.mu/target_semimajor_axis);

while (1)




end






da_dt = 2*a^2*v*accel_thrust*cos(phi_a)/mu;
de_dt = (accel_thrust/v) *(2*(e+cus(nu))*cos(phi) + (r/a)*sin(nu)*sin(phi_e));
di_dt = accel_thrust*cos(theta)*sin(beta)/h;

function [gamma] = getFlightPathAngle(r,v)
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
    phi_e = atan2(r*sin(nu),2*a*(e+cos(nu)));
    beta = pi*cos(omega + nu)/2;
end

function [c_a, c_e, c] = getBlendedCvector(gamma, phi_a, phi_e)
    % gamma = fpa, beta = yaw angle
    c_a = [sin(gamma + phi_a); cos(gamma + phi_a); 0];
    c_e = [sin(gamma + phi_e); cos(gamma +  phi_e); 0];
    
    c = (G_a*c_a + G_e*c_e)/(norm(G_a*c_a)+norm(G_e*c_e));
end 

function [alpha, beta] = getPitchYaw(c, theta)
    alpha = atan2(c(1),c(2));   
    beta = G_i*cos(theta);
end

function u_hat_rtn = getThrustVector(alpha, beta)
    u_hat_rtn = [sin(alpha)*cos(beta)]; cos(alpha)*cos(beta); sin(beta);
end

function [t, state_dot] = systemDynamics(t, state)
    % State vector: [x, y, z, vx, vy, vz, mass]
    r = state(1:3);
    v = state(4:6);
    m = state(7);

 % Gravitational acceleration
    r_mag = norm(r);
    a_grav = -thrust_params.mu * r / (r_mag^3);
    
    % variational orbital elements
    da_dt = 2*a^2*v*accel_thrust*cos(phi_a)/mu;
    de_dt = (accel_thrust/v) *(2*(e+cus(nu))*cos(phi) + (r/a)*sin(nu)*sin(phi_e));
    di_dt = accel_thrust*cos(theta)*sin(beta)/h;

    % Thrust acceleration
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