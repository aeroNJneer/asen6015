% Define gravitational parameter (mu) for Earth in km^3/s^2
mu = 398600;

% Initial and final position vectors (in km)
r1 = [6578, 0, 0]; % Example initial position
r2 = [-42004, 3674.8, 0]; % Example final position

% Time of flight (in seconds)
tof = 37864; % Example time of flight

% Solve Lambert's problem
[v1, v2] = lambert_univ(mu,r1, r2, tof, 1);

% Initial velocity vector (in km/s)
v_initial = [0, 7.7483, 0]; % Example initial velocity

% Calculate the delta-v for the first burn
delta_v1 = v1 - v_initial;

% Calculate the delta-v for the second burn
delta_v2 = v2 - v_initial;

% Display the results
fprintf('Delta-v for the first burn: [%f, %f, %f] km/s\n', delta_v1);
fprintf('Delta-v for the second burn: [%f, %f, %f] km/s\n', delta_v2);

% Time vector for simulation
t = linspace(0, tof, 1000);

% Simulate the orbits
orbit1 = simulate_orbit(r1, v_initial, t, mu);
transfer_orbit = simulate_orbit(r1, v1, t, mu);
orbit2 = simulate_orbit(r2, v2, linspace(t(end),t(end)+5000,1000), mu);

% Plot the orbits
figure;
plot(orbit1(:,1), orbit1(:,2), 'b', 'DisplayName', 'Initial Orbit');
hold on;
plot(transfer_orbit(:,1), transfer_orbit(:,2), 'r--', 'DisplayName', 'Transfer Orbit');
 plot(orbit2(:,1), orbit2(:,2), 'g', 'DisplayName', 'Final Orbit');
legend;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Orbit Simulation with Lambert Transfer');
grid on;
axis equal;

% Function to solve Lambert's problem
function [v1, v2] = lambert(r1, r2, tof, mu)
    % This function uses the universal variable formulation to solve Lambert's problem
    % For simplicity, we'll use a built-in function from MATLAB Central File Exchange
    % You can replace this with your own implementation if needed
    [v1, v2] = lambert_universal_variable(r1, r2, tof, mu);
end

% Function to simulate the orbit
function orbit = simulate_orbit(r0, v0, t, mu)
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    [~, state] = ode45(@(t, y) two_body(t, y, mu), t, [r0, v0], options);
    orbit = state(:, 1:3);
end

% Two-body equations of motion
function dydt = two_body(~, y, mu)
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);
    a = -mu / r_norm^3 * r;
    dydt = [v; a];
end
