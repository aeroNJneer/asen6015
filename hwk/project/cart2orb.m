% Utility function to convert state vectors to orbital elements
function [a, e, i, omega, RAAN, true_anomaly, mu] = cart2orb(r, v)
    mu = 3.986004418e14;  % Earth gravitational constant, m^3/s^2
    
    % Angular momentum
    h_vec = cross(r, v);
    h = norm(h_vec);
    
    % Eccentricity vector
    e_vec = cross(v, h_vec)/mu - r/norm(r);
    e = norm(e_vec);
    
    % Semimajor axis
    specific_energy = norm(v)^2/2 - mu/norm(r);
    a = -mu / (2*specific_energy);
    
    % Inclination
    i = acos(h_vec(3) / h);
    
    % Node vector
    n_vec = cross([0;0;1], h_vec);
    
    % RAAN
    if norm(n_vec) > 0
        RAAN = atan2(n_vec(2), n_vec(1));
    else
        RAAN = 0;
    end
    
    % Argument of periapsis
    if norm(n_vec) > 0
        omega = atan2(dot(e_vec, cross(n_vec, r)), ...
                      dot(n_vec, e_vec) * norm(r));
    else
        omega = 0;
    end
    % True anomaly
    % Compute the angle between the eccentricity vector and position vector
    true_anomaly = atan2(dot(h_vec, cross(e_vec, r)) / (h * e), ...
          dot(e_vec, r) / (e * norm(r)));
end
