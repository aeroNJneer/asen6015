function RAAN = calculate_RAAN(a, e, i)
    % Constants
    mu = 398600; % Earth's gravitational parameter, km^3/s^2

    % Calculate the mean motion (n)
    n = sqrt(mu / a^3);

    % Calculate the specific angular momentum (h)
    h = sqrt(mu * a * (1 - e^2));

    % Calculate the node line vector (N)
    N = [-h * sin(i); h * cos(i); 0];

    % Calculate the RAAN
    RAAN = atan2(N(2), N(1));
    
  
    % Ensure RAAN is in the range [0, 360]
    if RAAN < 0
        RAAN = RAAN + 2*pi;
    end
end