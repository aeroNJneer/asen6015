% Constants
mu = 3.986e14; % Earth's gravitational parameter (m^3/s^2)
r_LEO = 6678e3; % Radius of LEO (m)
r_GEO = 42164e3; % Radius of GEO (m)
Isp = 3000; % Specific impulse of the ion engine (s)
g0 = 9.81; % Standard gravity (m/s^2)
T_max = 0.5; % Maximum thrust (N)
m0 = 500; % Initial mass of the spacecraft (kg)

% Time settings
dt = 60; % Time step (s)
t_max = 1e7; % Maximum simulation time (s)

% Initial conditions
r = r_LEO; % Initial radius (m)
v = sqrt(mu/r_LEO); % Initial velocity (m/s)
m = m0; % Initial mass (kg)
t = 0; % Initial time (s)

% Arrays to store results
time = [];
radius = [];
velocity = [];
mass = [];

% Simulation loop
while r < r_GEO && t < t_max
    % Calculate acceleration
    a = T_max / m;
    
    % Update velocity and radius
    v = v + a * dt;
    r = r + v * dt;
    
    % Update mass
    m_dot = T_max / (Isp * g0);
    m = m - m_dot * dt;
    
    % Update time
    t = t + dt;
    
    % Store results
    time = [time; t];
    radius = [radius; r];
    velocity = [velocity; v];
    mass = [mass; m];
end

% Plot results
figure;
subplot(3,1,1);
plot(time/3600, radius/1e3);
xlabel('Time (hours)');
ylabel('Radius (km)');
title('Orbital Radius vs Time');

subplot(3,1,2);
plot(time/3600, velocity/1e3);
xlabel('Time (hours)');
ylabel('Velocity (km/s)');
title('Orbital Velocity vs Time');

subplot(3,1,3);
plot(time/3600, mass);
xlabel('Time (hours)');
ylabel('Mass (kg)');
title('Spacecraft Mass vs Time');
