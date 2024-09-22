% Constants
init;
mu = Gearth; % Earth's gravitational parameter, m^3/s^2
burn_acceleration = 30; % Constant acceleration during burn, m/s^2
dv1 = 2.23;

a = 7000e3; % Semi-major axis, meters
e = 0.01; % Eccentricity
inc = deg2rad(13); % Inclination, radians
RAAN = 0; % Right Ascension of Ascending Node, radians
w = 0; % Argument of Perigee, radians
M0 = 0; % Mean Anomaly at Epoch, radians
burn_duration = dv1/burn_acceleration; % Duration of the burn, seconds

% Time vector
tspan = [0 burn_duration];

% Initial state vector [rx ry rz vx vy vz]
[r0, v0] = keplerian_to_cartesian(a, e, inc, RAAN, w, M0, mu);
X0 = [r0; v0];

% ODE options
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

% Solve the ODE
[t, state] = ode45(@(t, y) satellite_dynamics(t, y, burn_acceleration, mu), tspan, X0, options);

% Extract position and velocity
r = state(:, 1:3);
v = state(:, 4:6);

% Plot the orbit
figure;
plot3(r(:,1), r(:,2), r(:,3));
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Satellite Orbit with Finite Burn');
grid on;

% Function to convert Keplerian elements to Cartesian coordinates
function [r, v] = keplerian_to_cartesian(a, e, i, RAAN, w, M0, mu)
    % Mean anomaly to eccentric anomaly
    init;
    E = FindEccAnomFromMeanAnom(M0,e);
    % True anomaly
    nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    theta = nu + w;
    % Distance
    r_mag = RfromEA(a,e,E);
    % Position in orbital plane
    r_orbital = r_mag * [cos(nu); sin(nu); 0];
    % Velocity in orbital plane
    v_orbital = sqrt(mu*a)/r_mag * [-sin(E); sqrt(1-e^2)*cos(E); 0];
    % Rotation matrices
    R3_W = [cos(-w) sin(-w) 0; -sin(-w) cos(-w) 0; 0 0 1];
    R1_i = [1 0 0; 0 cos(-i) sin(-i); 0 -sin(-i) cos(-i)];
    R3_RAAN = [cos(-RAAN) sin(-RAAN) 0; -sin(-RAAN) cos(-RAAN) 0; 0 0 1];
    % Combined rotation matrix
    Q = R3_RAAN * R1_i * R3_W;

    Q = DCM(RAAN,i,theta);
    % Position and velocity in inertial frame
    r = Q * r_orbital;
    v = Q * v_orbital;
end

% Function to compute satellite dynamics
function dydt = satellite_dynamics(t, y, burn_acceleration, mu)
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);
    a_gravity = -mu / r_norm^3 * r;
    a_burn = burn_acceleration * v / norm(v);
    dydt = [v; a_gravity + a_burn];
end
