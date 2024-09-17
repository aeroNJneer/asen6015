function compare_burns(deltaV,duration,acceleration)
    % Constants
    mu = 398600; % Gravitational parameter for Earth (km^3/s^2)
    r0 = 7000; % Initial orbit radius (km)
    v0 = sqrt(mu / r0); % Initial orbit velocity (km/s)
    % burn_duration  Duration of the finite burn (s)
    % acceleration  Constant acceleration (km/s^2)

    % Normalize the burn vector
    burn_vector = burn_vector / norm(burn_vector);

    % Time settings
    t_impulsive = 0; % Impulsive burn time
    t_finite = linspace(0, burn_duration, 1000); % Finite burn time array

    % Impulsive burn
    v_impulsive = v0 + deltaV * burn_vector(2); % Only tangential component affects velocity
    r_impulsive = r0 + deltaV * burn_vector(1); % Only radial component affects position

    % Finite burn
    v_finite = v0 + acceleration * t_finite * burn_vector(2); % Tangential velocity change
    r_finite = r0 + 0.5 * acceleration * t_finite.^2 * burn_vector(1); % Radial position change

    % Plotting
    figure;
    hold on;
    plot_orbit(r_impulsive, v_impulsive, mu, 'Impulsive Burn');
    plot_orbit(r_finite, v_finite, mu, 'Finite Burn');
    legend('Impulsive Burn', 'Finite Burn');
    title('Comparison of Impulsive and Finite Burns');
    xlabel('X (km)');
    ylabel('Y (km)');
    axis equal;
    grid on;
    hold off;
end

function plot_orbit(r, v, mu, label)
    % Calculate orbital elements
    a = r / (2 - (r * v^2 / mu));
    e = sqrt(1 - (r * v^2 / mu));
    theta = linspace(0, 2 * pi, 1000);
    r_orbit = a * (1 - e^2) ./ (1 + e * cos(theta));

    % Convert to Cartesian coordinates
    x = r_orbit .* cos(theta);
    y = r_orbit .* sin(theta);

    % Plot orbit
    plot(x, y, 'DisplayName', label);
end
