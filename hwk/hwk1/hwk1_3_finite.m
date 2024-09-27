clear
close all
init;
% Constants
mu = Gearth;
dt = 10; % Time step (s)

% Initial conditions
r = [6571 ; 0; 0];              % Initial position at 200 km altitude (km)
v = [0; sqrt(mu / norm(r)); 0]; % Initial velocity for circular orbit (km/s)
a_i = 0.03;                     % Constant acceleration during burn (km/s^2)

% Hohmann transfer parameters
r1 = 6571 ;                     % Initial orbit radius (m)
r2 = 42157 ;                    % Final orbit radius (m)
v1 = sqrt(mu / r1);             % Initial orbit velocity (m/s)
v2 = sqrt(mu / r2);             % Final orbit velocity (m/s)
a_t = (r1+r2)/2;                % semi major axis of transfer orbit
v_transfer1 = VfromRA(mu,r1,a_t);   % Initial velocity of transfer orbit (km/s)
v_transfer2 = VfromRA(mu,r2,a_t);   % Final velocity of transfer orbit (km/s)

t_transfer = PeriodfromA(mu,a_t)/2 ; % Time of transfer
t_orbit1 = PeriodfromA(mu,r1);    % Period of initial orbit
t_orbit2 = PeriodfromA(mu,r2);    % Period of final orbit

total_time = t_orbit2 + t_transfer + t_orbit1; % Total simulation time (s)

% calculate max burn time by dividing delta-v by max acceleration
burn1 = (v_transfer1-v1)/a_i;       
burn2 = (v2-v_transfer2)/a_i;

% Preallocate arrays for position and velocity
num_steps = floor(total_time / dt);
positions = zeros(3, num_steps);
velocities = zeros(3, num_steps);
tspan = linspace(0,total_time,num_steps);
j = 1;
% Simulation loop
for i = 1:length(tspan)
    % Calculate gravitational force
    r_mag = norm(r);
    F_gravity = -mu / r_mag^3 * r;
    
    % Apply Hohmann transfer burns
    if tspan(i) >= floor(t_orbit1 )  &&  tspan(i) < floor(t_orbit1+ 0.39*(burn1) )  % start first burn at 200 km altitude for 1/3 max burn time
       a_burn = a_i*v/norm(v) + [0.001; 0; 0];             % acceleration in the direction of the velocity vector (+ thrust vector error)
    elseif norm(r) >= 0.62*r2 && norm(v) < v2          % when at nearing target and velocity less than requiredi) 
       a_burn = a_i*v/norm(v) + [0.001; 0; 0];
    else
       a_burn = [0; 0; 0];
    end

    % Update velocity and position
    a_total = F_gravity + a_burn;
    v = v + a_total * dt;
    r = r + v * dt;
    
    % Store position and velocity
    positions(:, i) = r;
    velocities(:, i) = v;

    % Update velocity and position
    a_total = F_gravity + a_burn;
    v = v + a_total * dt;
    r = r + v * dt;
    
    % Store position and velocity
    positions(:, i) = r;
    velocities(:, i) = v;
end

figure(1);
plot3(positions(1, :), positions(2, :), positions(3,:));
xlabel('X (m)');
ylabel('Y (m)');
title('3D Orbit Hohmann Transfer (radial thrust vector error)');
grid on;
