function B = calculate_B_matrix_gauss(a, e, i, RAAN, omega, M, mu)
    % This function calculates the B matrix for the dynamic system of a spacecraft in orbit
    % using Gauss's planetary equations and classical orbital elements.
    %
    % Inputs:
    %   a     - Semi-major axis (km)
    %   e     - Eccentricity
    %   i     - Inclination (rad)
    %   RAAN  - Right Ascension of the Ascending Node (rad)
    %   omega - Argument of Periapsis (rad)
    %   M     - Mean Anomaly (rad)
    %   mu    - Gravitational parameter (km^3/s^2)
    %
    % Output:
    %   B - B matrix (6x3)

    % Calculate the true anomaly from the mean anomaly
    E = mean_to_eccentric_anomaly(M, e);
    nu = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

    % Calculate the orbital radius
    r = a * (1 - e^2) / (1 + e * cos(nu));

    % Calculate the specific angular momentum
    h = sqrt(mu * a * (1 - e^2));

    % Rotation matrix from perifocal to ECI frame
    cos_RAAN = cos(RAAN);
    sin_RAAN = sin(RAAN);
    cos_omega = cos(omega);
    sin_omega = sin(omega);
    cos_i = cos(i);
    sin_i = sin(i);

    R = [cos_RAAN*cos_omega - sin_RAAN*sin_omega*cos_i, -cos_RAAN*sin_omega - sin_RAAN*cos_omega*cos_i, sin_RAAN*sin_i;
         sin_RAAN*cos_omega + cos_RAAN*sin_omega*cos_i, -sin_RAAN*sin_omega + cos_RAAN*cos_omega*cos_i, -cos_RAAN*sin_i;
         sin_omega*sin_i, cos_omega*sin_i, cos_i];

    % Calculate the partial derivatives of the orbital elements with respect to the perturbing accelerations
    p = a * (1 - e^2);
    f_r = [2*a^2*e*sin(nu)/h; p*sin(nu)/h; 0];
    f_theta = [2*a^2*(1 + e*cos(nu))/h; (p + r)*cos(nu)/h; 0];
    f_h = [0; 0; r*cos(nu)/h];

    % Transform the partial derivatives to the ECI frame
    f_r_ECI = R * f_r;
    f_theta_ECI = R * f_theta;
    f_h_ECI = R * f_h;

    % Construct the B matrix
    B = [f_r_ECI, f_theta_ECI, f_h_ECI];

    % Display the B matrix
    disp('B matrix:');
    disp(B);function B = calculate_B_matrix_gauss(a, e, i, RAAN, omega, M, mu)
    % This function calculates the B matrix for the dynamic system of a spacecraft in orbit
    % using Gauss's planetary equations and classical orbital elements.
    %
    % Inputs:
    %   a     - Semi-major axis (km)
    %   e     - Eccentricity
    %   i     - Inclination (rad)
    %   RAAN  - Right Ascension of the Ascending Node (rad)
    %   omega - Argument of Periapsis (rad)
    %   M     - Mean Anomaly (rad)
    %   mu    - Gravitational parameter (km^3/s^2)
    %
    % Output:
    %   B - B matrix (6x3)

    % Calculate the true anomaly from the mean anomaly
    E = mean_to_eccentric_anomaly(M, e);
    nu = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

    % Calculate the orbital radius
    r = a * (1 - e^2) / (1 + e * cos(nu));

    % Calculate the specific angular momentum
    h = sqrt(mu * a * (1 - e^2));

    % Rotation matrix from perifocal to ECI frame
    cos_RAAN = cos(RAAN);
    sin_RAAN = sin(RAAN);
    cos_omega = cos(omega);
    sin_omega = sin(omega);
    cos_i = cos(i);
    sin_i = sin(i);

    R = [cos_RAAN*cos_omega - sin_RAAN*sin_omega*cos_i, -cos_RAAN*sin_omega - sin_RAAN*cos_omega*cos_i, sin_RAAN*sin_i;
         sin_RAAN*cos_omega + cos_RAAN*sin_omega*cos_i, -sin_RAAN*sin_omega + cos_RAAN*cos_omega*cos_i, -cos_RAAN*sin_i;
         sin_omega*sin_i, cos_omega*sin_i, cos_i];

    % Calculate the partial derivatives of the orbital elements with respect to the perturbing accelerations
    % p = a * (1 - e^2);
    % f_r = [2*a^2*e*sin(nu)/h; p*sin(nu)/h; 0];
    % f_theta = [2*a^2*(1 + e*cos(nu))/h; (p + r)*cos(nu)/h; 0];
    % f_h = [0; 0; r*cos(nu)/h];

    p = a * (1 - e^2);
    f_r = [(2*a^2/h)*(e*sin(TA)); p*sin(TA)/h; 0; 0; ];
    f_theta = [2*a^2*(1 + e*cos(nu))/h; (p + r)*cos(nu)/h; 0];
    f_h = [0; 0; r*cos(nu)/h];


dB = [(2*a^2/h)*(e*sin(TA)) (2*a^2/h)*p/r 0
    p*sin(TA)/h (p+r)*cos(TA)+(r*e)   0
    0 0 (r*cos(TA+w))/h
    0 0 r*sin(TA+w)/(h*sin(inc))
    -p*cos(TA)/(h*e) (p+r)*sin(TA)/(h*e) -r*sin(TA+w)*cos(inc)/(h*sin(inc))
    n+a*sqrt(1-e^2)/(a*h*e)*(p*cos(TA)-(2*r*e)) a*sqrt(1-e^2)/(a*h*e)*(p+r)*sin(TA) 0];

F_vec = dB'


dMt = dM(a,n,t);


    % Transform the partial derivatives to the ECI frame
    f_r_ECI = R * f_r;
    f_theta_ECI = R * f_theta;
    f_h_ECI = R * f_h;

    % Construct the B matrix
    B = [f_r_ECI, f_theta_ECI, f_h_ECI];

    % Transform the B matrix for the orbital element state vector
    B_transformed = zeros(6, 3);
    B_transformed(1, :) = B(1, :); % Semi-major axis
    B_transformed(2, :) = B(2, :); % Eccentricity
    B_transformed(3, :) = B(3, :); % Inclination
    B_transformed(4, :) = B(4, :); % RAAN
    B_transformed(5, :) = B(5, :); % Argument of Periapsis
    B_transformed(6, :) = B(6, :); % Mean Anomaly

    % Display the transformed B matrix
    disp('Transformed B matrix:');
    disp(B_transformed);
    B = B_transformed;
end

function E = mean_to_eccentric_anomaly(M, e)
    % This function calculates the eccentric anomaly from the mean anomaly using
    % an iterative method (Newton-Raphson).
    %
    % Inputs:
    %   M - Mean Anomaly (rad)
    %   e - Eccentricity
    %
    % Output:
    %   E - Eccentric Anomaly (rad)

    % Initial guess
    E = M;
    tol = 1e-10;
    max_iter = 1000;
    iter = 0;

    % Newton-Raphson iteration
    while iter < max_iter
        E_new = E - (E - e*sin(E) - M) / (1 - e*cos(E));
        if abs(E_new - E) < tol
            break;
        end
        E = E_new;
        iter = iter + 1;
    end
end

end

function E = mean_to_eccentric_anomaly(M, e)
    % This function calculates the eccentric anomaly from the mean anomaly using
    % an iterative method (Newton-Raphson).
    %
    % Inputs:
    %   M - Mean Anomaly (rad)
    %   e - Eccentricity
    %
    % Output:
    %   E - Eccentric Anomaly (rad)

    % Initial guess
    E = M;
    tol = 1e-10;
    max_iter = 1000;
    iter = 0;

    % Newton-Raphson iteration
    while iter < max_iter
        E_new = E - (E - e*sin(E) - M) / (1 - e*cos(E));
        if abs(E_new - E) < tol
            break;
        end
        E = E_new;
        iter = iter + 1;
    end
end
