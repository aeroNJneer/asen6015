function [t, y] = earth_shadow_orbital_dynamics(t_span, initial_state, params)
% EARTH_SHADOW_ORBITAL_DYNAMICS Simulate satellite orbit with shadow effects
%   Implements Keplerian dynamics with thrust during sunlight and coast during shadow
%
%   Inputs:
%   - t_span: Time span for integration [t_start, t_end]
%   - initial_state: Initial state vector [a, e, i, Omega, omega, M]
%   - params: Struct containing mission parameters
%
%   Outputs:
%   - t: Time vector
%   - y: State history matrix

    % Constants
    mu = 3.986e14;  % Earth gravitational parameter (m^3/s^2)
    R_earth = 6371e3;  % Earth radius (m)
    R_sun = 696340e3;  % Sun radius (m)
    dist_earth_sun = 149.6e9;  % Earth-Sun distance (m)

    % Unpack parameters
    Isp = params.Isp;  % Specific impulse (s)
    g0 = 9.81;  % Standard gravity (m/s^2)
    max_thrust = params.max_thrust;  % Maximum thrust (N)
    mass = params.initial_mass;  % Initial spacecraft mass (kg)

    % ODE options
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @shadow_event);

    % Initialize storage
    t = [];
    y = [];
    current_mass = mass;
    current_state = initial_state;

    % Time integration loop
    while t(end) < t_span(2)
        % Propagate orbit
        [t_segment, y_segment, te, ye, ie] = ode45(@orbit_dynamics, t_span, current_state, options);

        % Store results
        t = [t; t_segment];
        y = [y; y_segment];

        % Check if reached end of simulation
        if isempty(te)
            break;
        end

        % Determine if in shadow or sunlight
        switch ie
            case 1  % Entering shadow
                % Coast in shadow
                current_state = y_segment(end, :);
            case 2  % Exiting shadow
                % Apply thrust when exiting shadow
                current_state = y_segment(end, :);
                
                % Calculate thrust
                [~, thrust_vec] = orbit_dynamics(te, current_state);
                
                % Mass update
                dm = max_thrust / (Isp * g0);
                current_mass = current_mass - dm * (te - t_segment(end));
        end

        % Update time span for next segment
        t_span = [te, t_span(2)];
    end

    % Orbital dynamics function
    function [dx, thrust_vec] = orbit_dynamics(t, x)
        % Extract orbital elements
        a = x(1);
        e = x(2);
        i = x(3);
        Omega = x(4);
        omega = x(5);
        M = x(6);

        % Compute additional orbital parameters
        n = sqrt(mu / a^3);  % Mean motion
        
        % Compute thrust vector (example simple thrust direction)
        thrust_mag = max_thrust;
        thrust_vec = [thrust_mag * cos(omega), thrust_mag * sin(omega), 0];

        % Basic Keplerian dynamics without perturbations
        dx = zeros(6, 1);
        dx(1) = 2 * a^2 / (mu * sqrt(1 - e^2)) * thrust_vec(1);  % da/dt
        dx(2) = sqrt(1 - e^2) / (mu * a) * (thrust_vec(2) * sin(M) + thrust_vec(1) * (1 + e * cos(M)));  % de/dt
        dx(3) = 0;  % Simple case: no inclination change
        dx(4) = 0;  % No right ascension change
        dx(5) = 0;  % No argument of periapsis change
        dx(6) = n;  % Mean anomaly progression
    end

    % Shadow detection event function
    function [value, isterminal, direction] = shadow_event(t, x)
        % Simplified shadow detection 
        % Assumes circular orbit and simplified geometry
        
        % Compute eclipse condition
        umbra_angle = asin((R_earth) / (dist_earth_sun));
        
        % You would replace this with a more accurate shadow calculation
        % This is a simplified approximation
        value = [1, -1];  % Two events: entering and exiting shadow
        isterminal = [1, 1];  % Stop integration at both events
        direction = [0, 0];  % Detect crossings in both directions
    end
end