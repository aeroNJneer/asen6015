% ASEN 6015 
% Hwk 4

clear;
close all;
init;

% Constants
mu = Gearth; 
thrust_acc = 0.030; % Thrust acceleration in km/s^2
tf = 23; % Initial guess for burn-out time in seconds

% Function to transform ECI to LVLH frame
eci_to_lvlh = @(r, v) [r/norm(r), cross(v, r)/norm(cross(v, r)), cross(cross(v, r), r)/norm(cross(cross(v, r), r))];

% Initial state vectors
r0 = [8276; 5612; 5]; % Position in km
v0 = [-3.142; 4.672; 0]; % Velocity in km/s
x0 = [r0;v0];

% Initial COEs
a_init = AfromVR(mu,norm(v0),norm(r0));
h_init = cross(r0,v0);
E_init = EccfromHRV(mu,r0,v0,h_init);
e_init = norm(E_init);
i_init = IfromH(h_init);
LON = cross([0 0 1],h_init);
n_mag = norm(LON);  
Raan_init = RaanfromLON(LON) * raan_sign(LON);
w_init = WfromLONe(LON,E_init)*w_sign(E_init);
TA_init = TAfromRe(r0,E_init);

% Target orbital elements
a_target = 10000; % Semi-major axis in km
i_target = 0.0; % Inclination in degrees
e_target = 0.0001; % Eccentricity
w_target = 35; % Argument of perigee in degrees
h_target = HfromAe(mu,a_target,e_target);
Raan_target = 0;
rp_target = Rp(a_target,e_target);
C = DCM(Raan_target,i_target,w_target);     % DCM
Rp_target = C * [rp_target; 0; 0];
Vp_target = C * [0; VfromRA(mu,rp_target,a_target); 0];

% Initial guesses for guidance parameters
A = [1; 0; 0];
B = [0; 0; 0];
lambda_r = -B;
flamba_v = @(t) A+B*t;

X0 = [x0; lambda_r; flamba_v(0)];

% plot initial and final orbits
orbit_init = plot_orbit(a_init,e_init,i_init,Raan_init,w_init,TA_init);
orbit_target = plot_orbit(a_target,e_target,i_target,Raan_target,w_target,0);

%     figure(1) 
%     % Plot the orbits
%     plot3(orbit_init(1, :), orbit_init(2, :), orbit_init(3, :), 'b','DisplayName','Initial Orbit');
%      hold on;
%     plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r','DisplayName','Earth'); % Central body
%     plot3(r0(1),r0(2),r0(3),'bo','DisplayName','s/c at t_0')
%     plot3(orbit_target(1, :), orbit_target(2, :), orbit_target(3, :), 'g','DisplayName','Final Orbit');
%     xlabel('X (km)');
%     ylabel('Y (km)');
%     zlabel('Z (km)');
%     title('Orbit Plot');
%     legend('Location','bestoutside')
%     zlim([-20 20]);
%     hold off;
%% Question 1

% Simulate the burn without correction
% [r, v, t] = simulate_burn(r0, v0, A, B, thrust_acc, tf, mu);
[r, v, t, lambda_r, lambda_v] = propagate_with_costate_LVLH_ECI(X0,A,B,C,thrust_acc,tf,mu);
% coasting
Xc = [r(end,:)'; v(end,:)'; lambda_r(end,:)'; lambda_v(end,:)'];
[rc, vc, tc, lambda_rc, lambda_vc] = propagate_with_costate_LVLH_ECI(Xc,A,B,C,0,60*45,mu);

% Plot the resulting trajectory
    figure(2);
    plot3(r(:,1), r(:,2), r(:,3));
    hold on;
    plot3(rc(:,1), rc(:,2), rc(:,3));
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Trajectory of the Spacecraft with no Correcction');
    hold off
    
%% Question 2

% constraints

    % Implement predictor-corrector guidance system
    [A_corr, B_corr, tf_corr] = predictor_corrector(r0, v0, A, B, thrust_acc, tf, mu, a_target, e_target, i_target, w_target);

    % Monte Carlo simulation with thrust dispersion
    monte_carlo_simulation(r0, v0, A_corr, B_corr, thrust_acc, tf_corr, mu, a_target, e_target, i_target, w_target);



function [r, v, t] = simulate_burn(r0, v0, A, B, thrust_acc, tf, mu)
    % Time span for integration
    tspan = linspace(0,tf,100);
    % Initial state vector
    y0 = [r0; v0];
    
    % Options for ODE solver
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    
    % Integrate the equations of motion
    [t, y] = ode45(@(t, y) satellite_dynamics(t, y, A, B, thrust_acc, mu), tspan, y0, options);
    
    % Extract position and velocity vectors
    r = y(:, 1:3);
    v = y(:, 4:6);
end

function dydt = satellite_dynamics(t, y, A, B, thrust_acc, mu)
    % Extract position and velocity vectors
    r = y(1:3);
    v = y(4:6);
    
    % Calculate thrust direction for each time step
    u = guidancelaw(A,B,t);     % LVLH

    % Calculate the acceleration due to thrust
    a_thrust = thrust_acc * u;  % LVLH
    
    % Calculate the gravitational acceleration
    a_gravity = -mu * r / norm(r)^3;    % ECI
    
    % Equations of motion
    drdt = v;       %ECI
    dvdt = a_gravity + a_thrust;   % ECI
    
    % Combine derivatives into a single column vector
    dydt = [drdt; dvdt];
end

function u = guidancelaw(A, B, t)
    % Calculate the thrust direction vector u
    At = A + B * t;
    u = At / norm(At);
end

function plot_trajectory(r, v, t)
    % Plot the trajectory in 3D space
    figure;
    plot3(r(:,1), r(:,2), r(:,3));
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Trajectory of the Spacecraft');
    grid on;
end

function [A_corr, B_corr, tf_corr] = predictor_corrector(r0, v0, A, B, thrust_acc, tf, mu, a_target, e_target, i_target, w_target)
    % Implement the predictor-corrector algorithm
    % Adjust the guidance parameters to minimize the errors

    % line search ... 

    % Return the corrected guidance parameters and burn-out time
end

function monte_carlo_simulation(r0, v0, A, B, thrust_acc, tf, mu, a_target, e_target, i_target, w_target)
    % Implement random thrust dispersion
    % Run multiple simulations and perform statistical analysis
    % Plot histograms of the results
end


function orbit_3D = plot_orbit(a, e, i, RAAN, w, theta)
    % a: semi-major axis (km)
    % e: eccentricity
    % i: inclination (degrees)
    % RAAN: right ascension of ascending node (degrees)
    % w: argument of periapsis (degrees)
    % theta: true anomaly (degrees)

    % Convert angles from degrees to radians
    i = deg2rad(i);
    RAAN = deg2rad(RAAN);
    w = deg2rad(w);
    theta = deg2rad(theta);

    % Define the number of points in the orbit
    num_points = 1000;
    theta_vals = linspace(0, 2*pi, num_points);

    % Calculate the radius for each point in the orbit
    r = a * (1 - e^2) ./ (1 + e * cos(theta_vals));

    % Calculate the position in the orbital plane
    x_orbit = r .* cos(theta_vals);
    y_orbit = r .* sin(theta_vals);

    % Rotation matrices
    R3_W = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
    R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_w = [cos(w) -sin(w) 0; sin(w) cos(w) 0; 0 0 1];

    % Combined rotation matrix
    R = R3_W * R1_i * R3_w;

    % Initialize the orbit in 3D space
    orbit_3D = zeros(3, num_points);

    for j = 1:num_points
        % Position in the orbital plane
        pos_orbit = [x_orbit(j); y_orbit(j); 0];
        % Rotate to 3D space
        pos_3D = R * pos_orbit;
        orbit_3D(:, j) = pos_3D;
    end
end

function dlambda = costate_dynamics(lambdav)
    dlamdav = B;
    dlambdar = lambdav


end


function [r, v, t, lambda_r, lambda_v] = propagate_with_costate_LVLH_ECI(x0, A, B, C, thrust_acc, tf, mu)
    % Time span for integration
    tspan = [0 tf];
    
    % Options for ODE solver
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    
    % Integrate the equations of motion and costate equations
    [t, y] = ode45(@(t, y) equations_of_motion_LVLH_ECI(t, y, A, B, C, thrust_acc, mu), tspan, x0, options);
    
    % Extract position, velocity, and costate vectors
    r = y(:, 1:3);
    v = y(:, 4:6);
    lambda_r = y(:, 7:9);
    lambda_v = y(:, 10:12);
end

function dydt = equations_of_motion_LVLH_ECI(t, y, A, B, C, thrust_acc, mu)
    % Extract position, velocity, and costate vectors
    r = y(1:3);
    v = y(4:6);
    lambda_r = y(7:9);
    lambda_v = y(10:12);
    
    % Calculate the LVLH frame transformation matrix
    R_LVLH_ECI = calculate_LVLH_to_ECI(r, v);
    
    % Transform position and velocity to LVLH frame
    r_LVLH = R_LVLH_ECI * r;
    v_LVLH = R_LVLH_ECI * v;
    
    % Calculate the thrust direction using bilinear tangent guidance law in LVLH frame
    u_LVLH = bilinear_tangent_guidance(A, B, t);
    
    % Transform thrust direction to ECI frame
    u_ECI = R_LVLH_ECI' * u_LVLH;
    
    % Calculate the acceleration due to thrust in ECI frame
    a_thrust = thrust_acc * u_ECI;
    
    % Calculate the gravitational acceleration in ECI frame
    r_norm = norm(r);
    a_gravity = -mu * r / r_norm^3;
    
    % Equations of motion in ECI frame
    drdt = v;
    dvdt = a_gravity + a_thrust;
    
    % Costate equations in ECI frame
    dlambda_rdt = A+B*t;    % lambda_v
    dlambda_vdt = B;
    
    % Combine derivatives into a single column vector
    dydt = [drdt; dvdt; dlambda_rdt; dlambda_vdt];
end

function R_LVLH_ECI = calculate_LVLH_to_ECI(r, v)
    % Calculate the unit vectors of the LVLH frame
    r_norm = norm(r);
    v_norm = norm(v);
    
    z_LVLH = r / r_norm;
    y_LVLH = cross(r, v) / norm(cross(r, v));
    x_LVLH = cross(y_LVLH, z_LVLH);
    
    % Transformation matrix from LVLH to ECI
    R_LVLH_ECI = [z_LVLH'; y_LVLH'; y_LVLH'];
end

function u = bilinear_tangent_guidance(A, B, t)
    % Calculate the thrust direction vector u
    At = A + B * t;
    u = At / norm(At);
end

function E = calculate_errors(x)
    rp = evalin('base','Rp_target');
    vp = evalin('base','Vp_target');

    nu = atan2


    h = cross(x(1:3),x(4:6));
    E(5) = dot(h,x(1:3));

end


% % Example usage:
% A = [1; 0; 0];
% B = [0; 0; 0];
% t = 0:0.1:23; % Time vector from 0 to 23 seconds
% 
% % Preallocate thrust direction array
% u = zeros(3, length(t));
% 
% % Calculate thrust direction for each time step
% for i = 1:length(t)
%     u(:, i) = bilinear_tangent_guidance(A, B, t(i));
% end
% 
% % Plot the thrust direction components over time
% figure;
% plot(t, u(1, :), 'r', t, u(2, :), 'g', t, u(3, :), 'b');
% xlabel('Time (s)');
% ylabel('Thrust Direction Components');
% legend('u_x', 'u_y', 'u_z');
% title('Thrust Direction Components Over Time');
% grid on;
