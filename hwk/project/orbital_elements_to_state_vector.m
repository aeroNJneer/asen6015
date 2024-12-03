function [r, v] = orbital_elements_to_state_vector(a, e, i, RAAN, omega, nu, mu)
    % Convert orbital elements to state vector
    % 
    % Inputs:
    %   a       - semimajor axis (km)
    %   e       - eccentricity (dimensionless)
    %   i       - inclination (radians)
    %   RAAN    - Right Ascension of Ascending Node (radians)
    %   omega   - Argument of periapsis (radians)
    %   nu      - True anomaly (radians)
    %   mu      - Gravitational parameter (km^3/s^2), default is Earth's if not provided
    %
    % Outputs:
    %   r       - Position vector in ECI frame (km)
    %   v       - Velocity vector in ECI frame (km/s)

    % Default to Earth's gravitational parameter if not specified
    if nargin < 7
        mu = 3.986004418e5;  % Earth's standard gravitational parameter (km^3/s^2)
    end

    % Calculate p (semi-latus rectum)
    p = a * (1 - e^2);

    % Calculate position in orbital (perifocal) frame
    r_pqw = [
        p * cos(nu) / (1 + e * cos(nu));
        p * sin(nu) / (1 + e * cos(nu));
        0
    ];

    % Calculate velocity in orbital (perifocal) frame
    v_pqw = [
        -sqrt(mu/p) * sin(nu);
        sqrt(mu/p) * (e + cos(nu));
        0
    ];

    % Rotation matrices
    % 1. Rotation about z-axis by -RAAN
    R_z_RAAN = [
        cos(RAAN)   -sin(RAAN)   0;
        sin(RAAN)    cos(RAAN)   0;
        0            0            1
    ];

    % 2. Rotation about x-axis by -i
    R_x_i = [
        1    0            0;
        0    cos(i)   -sin(i);
        0    sin(i)    cos(i)
    ];

    % 3. Rotation about z-axis by -(omega)
    R_z_omega = [
        cos(omega)   -sin(omega)   0;
        sin(omega)    cos(omega)   0;
        0            0             1
    ];

    % Combined rotation matrix (PQW to ECI)
    R_pqw_to_eci = (R_z_RAAN * R_x_i * R_z_omega)';

    % Transform position and velocity to ECI frame
    r = R_pqw_to_eci * r_pqw;
    v = R_pqw_to_eci * v_pqw;
end

% Example usage and validation function
function validate_orbital_elements_conversion()
    % Test case: GEO satellite
    mu = 3.986004418e5;  % Earth's gravitational parameter
    
    % Orbital elements
    a = 42164;           % semimajor axis (km)
    e = 0.0;             % eccentricity
    i = deg2rad(0);      % inclination
    RAAN = deg2rad(0);   % Right Ascension of Ascending Node
    omega = 0;           % Argument of periapsis
    nu = 0;              % True anomaly
    
    % Convert to state vector
    [r, v] = orbital_elements_to_state_vector(a, e, i, RAAN, omega, nu, mu);
    
    % Validate results
    disp('Position vector (km):');
    disp(r);
    disp('Velocity vector (km/s):');
    disp(v);
    
    % Compute orbital parameters from state vector for verification
    [a_check, e_check, i_check, RAAN_check, omega_check, nu_check] = state_vector_to_orbital_elements(r, v, mu);
    
    disp('Verification of orbital elements:');
    fprintf('Semimajor axis:     %.4f km (Original: %.4f)\n', a_check, a);
    fprintf('Eccentricity:       %.4f (Original: %.4f)\n', e_check, e);
    fprintf('Inclination:        %.4f deg (Original: %.4f)\n', rad2deg(i_check), rad2deg(i));
    fprintf('RAAN:               %.4f deg (Original: %.4f)\n', rad2deg(RAAN_check), rad2deg(RAAN));
    fprintf('Argument of periapsis: %.4f deg (Original: %.4f)\n', rad2deg(omega_check), rad2deg(omega));
    fprintf('True anomaly:       %.4f deg (Original: %.4f)\n', rad2deg(nu_check), rad2deg(nu));
end

% Complementary function to convert state vector back to orbital elements
function [a, e, i, RAAN, omega, nu] = state_vector_to_orbital_elements(r, v, mu)
    % Convert state vector to orbital elements
    
    % Default to Earth's gravitational parameter if not specified
    if nargin < 3
        mu = 3.986004418e5;  % Earth's standard gravitational parameter (km^3/s^2)
    end

    % Angular momentum vector
    h = cross(r, v);

    % Eccentricity vector
    e_vec = cross(v, h)/mu - r/norm(r);

    % Eccentricity magnitude
    e = norm(e_vec);

    % Node vector (perpendicular to orbital plane)
    n = cross([0, 0, 1], h);

    % Semimajor axis
    r_mag = norm(r);
    v_mag = norm(v);
    a = 1 / (2/r_mag - v_mag^2 / mu);

    % Inclination
    i = acos(h(3) / norm(h));

    % Right Ascension of Ascending Node
    if norm(n) ~= 0
        RAAN = atan2(n(2), n(1));
        if RAAN < 0
            RAAN = RAAN + 2*pi;
        end
    else
        RAAN = 0;
    end

    % Argument of periapsis
    if norm(n) ~= 0
        omega = atan2(dot(e_vec, cross(n, r/r_mag)), dot(n, e_vec));
        if omega < 0
            omega = omega + 2*pi;
        end
    else
        omega = 0;
    end

    % True anomaly
    nu = atan2(dot(h, cross(r, v)), dot(r, v));
    if nu < 0
        nu = nu + 2*pi;
    end
end