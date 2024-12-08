function RTN2ECI = Euler313(Raan, i, omega, true_anomaly)
    % Rotation matrix for RAAN
    theta = omega + true_anomaly;
    R_RAAN = [cos(Raan), -sin(Raan), 0;
              sin(Raan), cos(Raan), 0;
              0, 0, 1];
    
    % Rotation matrix for inclination
    R_i = [1, 0, 0;
           0, cos(i), -sin(i);
           0, sin(i), cos(i)];
    
    % Rotation matrix for argument of perigee
    R_omega = [cos(theta), -sin(theta), 0;
               sin(theta), cos(theta), 0;
               0, 0, 1];
    
    % Combined rotation matrix
    RTN2ECI = R_RAAN * R_i * R_omega;
end