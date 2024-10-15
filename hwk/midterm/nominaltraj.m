%initial conditions
x0 = -1.0;
y0 = -0.25
xdot0 = 0.003;
ydot0 = -0.0005;
mu = 3.986004418e5;

% final conditions
xf = -0.25 ;
yf = 0. ;
xdotf = 0.0002 ;
ydotf = 0. ;
Xf = [xf; yf; xdotf; ydotf];

r0 = [x0;y0;0];
v0 = [xdot0;ydot0];
rf = [xf;yf;0];
vf = [xdotf;ydotf];


    % Call the Superior Lambert Solver
    [v1, v2] = superior_lambert(r0, rf, 1200, mu)
    


function [v1, v2] = superior_lambert(r1, r2, tof, mu)
    % Implement the Superior Lambert Solver here
    
    % Calculate the distance between the initial and final positions
    R1 = norm(r1);
    R2 = norm(r2);
    dtheta = acos(dot(r1, r2) / (R1 * R2));
    
    % Guess initial values
    a_min = (R1 + R2) / 2; % Minimum semi-major axis
    beta = dtheta / 2; % Half the transfer angle
    
    % Iterate to find the semi-major axis and velocities
    a = a_min; % Initial guess
    max_iter = 100;
    tol = 1e-6;
    
    for iter = 1:max_iter
        % Compute the current time of flight
        [~, tof_current] = lambert_time_of_flight(a, mu, R1, R2, dtheta);
        
        % Check for convergence
        if abs(tof_current - tof) < tol
            break;
        end
        
        % Adjust semi-major axis
        a = a * (tof / tof_current);
    end
    
    % Calculate the velocities
    [v1, v2] = lambert_velocities(a, mu, r1, r2, dtheta);
end

function [E, tof] = lambert_time_of_flight(a, mu, R1, R2, dtheta)
    % Calculate the eccentric anomaly and time of flight for a given semi-major axis
    
    % Compute the semi-latus rectum and eccentricity
    p = (2 * a * R1 * R2 * cos(dtheta)) / (R1 + R2);
    e = sqrt(1 - p / a);
    
    % Calculate the eccentric anomaly
    E = acos((1 - R1 / a) / e);
    
    % Calculate the time of flight
    tof = sqrt(a^3 / mu) * (E - e * sin(E));
end

function [v1, v2] = lambert_velocities(a, mu, r1, r2, dtheta)
    % Calculate the velocities at the initial and final positions
    
    % Compute the velocity magnitudes
    V = sqrt(mu * (2 / norm(r1) - 1 / a));
    V2 = sqrt(mu * (2 / norm(r2) - 1 / a));
    
    % Compute the velocity vectors
    v1 = V * (r2 - r1) / norm(r2 - r1);
    v2 = V2 * (r1 - r2) / norm(r1 - r2);
end
