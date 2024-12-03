
% MATLAB script for low-thrust orbit transfer from LEO to GEO

% Constants
mu = 3.986e5; % Earth's gravitational parameter (km^3/s^2)
P = 10e3; % Input power to the engine (W)
Isp = 3300; % Specific impulse (s)
g = 9.81; % Earth's gravitational acceleration (m/s^2)-- convert
eta = 0.65; % Efficiency of the low-thrust engine
T = 2 * eta * P / (Isp * g); % Thrust (N)
m0 = 1200; % Initial mass of the spacecraft (kg)
r0 = 6378 + 550; % Initial radius (km)
v0 = sqrt(mu / r0); % Initial velocity (km/s)
a0 = r0; % Initial semi-major axis (km)
e0 = 0; % Initial eccentricity
i0 = 28.5; % Initial inclination (degrees)
aGEO = 42164; % Semi-major axis of GEO (km)

% determine state from oes


% create transfer object
transfer = LowThrustTransfer([a0, e0, i0], [aGEO, 0, 0], P, Isp, eta);

% Time parameters
dt = 60; % Time step (s)
t_max = 216.3 * 24 * 3600; % Maximum time (s)
t = 0:dt:t_max; % Time vector

% Initialize arrays
a = zeros(size(t));
e = zeros(size(t));
i = zeros(size(t));
m = zeros(size(t));
a(1) = a0;
e(1) = e0;
i(1) = i0;
m(1) = m0;

% Simulation loop
for k = 1:length(t)-1
    % Compute thrust acceleration
    aT = T / m(k);
    transfer.compute_thrust_direction(X,t)
%     % Compute optimal control laws
%     alpha_a = transfer.compute_semimajor_thrust_angle(transfer,X);
% 
%     alpha_e = transfer.compute_eccentricity_thrust_angle(X);
% 
%     beta = transfer.compute_yaw_steering_angle(X);
%     
%     % Blend control laws
%     Ga = 1;
%     [Ge, Gi] = transfer.compute_guidance_weights(t,phase);
%     [Ca, Ce] = transfer.
%     alpha = atan2(Ga * sin(alpha_a) + Ge * sin(alpha_e), Ga * cos(alpha_a) + Ge * cos(alpha_e));
%     beta = Gi * cos(i(k));
    
    % Update orbital elements
    a(k+1) = a(k) + 2 * a(k)^2 * aT * cos(alpha) * dt / v0;
    e(k+1) = e(k) + aT * (e(k) + cos(e(k))) * cos(alpha) * dt / v0;
    i(k+1) = i(k) + aT * cos(i(k)) * sin(beta) * dt / v0;
    
    % Update mass
    m(k+1) = m(k) - T * dt / (Isp * g);
end

% Plot results
figure;
subplot(3,1,1);
plot(t/3600/24, a);
xlabel('Time (days)');
ylabel('Semi-major axis (km)');
title('Semi-major axis vs Time');

subplot(3,1,2);
plot(t/3600/24, e);
xlabel('Time (days)');
ylabel('Eccentricity');
title('Eccentricity vs Time');

subplot(3,1,3);
plot(t/3600/24, i);
xlabel('Time (days)');
ylabel('Inclination (degrees)');
title('Inclination vs Time');

function dxdt = dynamics(t, state, acc)
    global mu
    % Unpack state
    r = state(1:3);
    v = state(4:6);
    mdot = -2*eta*P/(g*Isp)^2;

%     1) Compute the unit vectors ca and ce that define the in-plane optimal
% controls for maximum da=dt and de=dt , respectively. The unit
% vectors are expressed in a local rotating radial-transverse-normal
% (RTN) coordinate frame where the R axis is along the radial direction,
% the T axis is in the orbit plane along the transverse direction,
% and the N axis is normal to the orbit plane (however, the in-plane
% steering unit vectors ca and ce will only have R and T components).
% A general expression for the unit vector in the RTN frame is
% ck D sin ° C Á¤k ; cos ° C Á¤k ; 0
% T
% ; k D a; e (8)
% where Á¤a and Á¤e are computed using the optimal controls (4) and
% (5). The angle ° is the  ight-path angle and is measured from the
% local horizon to the velocity vector.
   

% Calculate LVLH frame
    R_eci2lvlh = calculateLVLHFrame(r, v);
    R_lvlh2eci = R_eci2lvlh';
    
    % Current control direction in LVLH
    u_lvlh = (A + B*t)/norm(A + B*t);
    
    % Convert control to ECI
    u_eci = R_lvlh2eci * u_lvlh;
    
    % Gravitational acceleration
    g_eci = -mu*r/norm(r)^3;
    
    % State derivatives
    dr = v;
    dv = g_eci + acc*u_eci;
    
    dxdt = [dr; dv; mdot];
end
% 
% function R_eci2rtn = calculateRTNFrame(r, v)
%     % Calculate LVLH frame vectors
%     N = cross(r, v);
%     
%     % Unit vectors
%     r_hat = r / norm(r); % z
%     n_hat = N / norm(N); % y
%     t_hat = cross(h_hat, r_hat); % x
%     
%     % Rotation matrix from ECI to LVLH
%     % gamma = fpa
%     c_a = [sin(gamma + phi_star_a); cos(gamma + phi_star_a); 0];
%     c_e = [sin(gamma + phi_star_e); cos(gamma + phi_star_e); 0]
%     R_eci2RTN = [r_hat, t_hat, n_hat];
% end
% 
