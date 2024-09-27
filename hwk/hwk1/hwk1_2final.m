clear
%close all
clc
init;
% Constants
mu = Gearth; 

% Initial and desired orbit elements
a0 = 7000; % Semi-major axis (m)
e0 = 0.01; % Eccentricity
i0 = deg2rad(30); % Inclination (rad)
raan0 = deg2rad(40); % Right ascension of ascending node (rad)
w0 = deg2rad(60); % Argument of perigee (rad)
M0 = deg2rad(0); % Mean anomaly (rad)

af = 10000; % Target semi-major axis (m)
ef = 0.005; % Target eccentricity
incf = deg2rad(25); % Target inclination (rad)
raanf = deg2rad(35); % Target right ascension of ascending node (rad)
wf = deg2rad(55); % Target argument of perigee (rad)
Mf = deg2rad(0); % Target mean anomaly (rad)

% Gain matrix
K = diag([1e-7, 1e-6, 1e-5, 1e-7, 1e-7, 1e-7]);

% Time settings
dt = 10; % Time step (s)
total_time = 60*60; % Total simulation time (s)
num_steps = total_time / dt;

% Preallocate arrays for orbit elements
a = zeros(1, num_steps);
e = zeros(1, num_steps);
i = zeros(1, num_steps);
raan = zeros(1, num_steps);
w = zeros(1, num_steps);
M = zeros(1, num_steps);

% Initial orbit elements
a(1) = a0;
e(1) = e0;
i(1) = i0;
raan(1) = raan0;
w(1) = w0;
M(1) = M0;

% Desired orbit elements
xd = [af; ef; incf; raanf; wf; Mf];

% Simulation loop
for t = 2:num_steps
    % Current orbit elements
    x = [a(t-1); e(t-1); i(t-1); raan(t-1); w(t-1); M(t-1)];
    EA = FindEccAnomFromMeanAnom(x(6),x(2));
    TA = TAfromEA(x(2),EA);
    % Control law
    u = -K * (x - xd);
    p = PfromAe(a(t-1),e(t-1));
    h = HfromP(mu,p);
    n = MMfromP(p);
    r = RfromEA(a(t-1),e(t-1),EA);

    [~, v] = perifocal_sv_from_coe(x,mu);
    % Gauss' Variational equations
    da_dt = (2*a(t-1)^2/h) * (e(t-1)*sin(TA)*u(1) + p*u(2)/r);
    de_dt = ((p*sin(TA)*u(1) + ((p+r)*cos(TA) + r*e(t-1))*u(2)))/h;
    di_dt = r * cos(TA+w(t-1)) * u(3)/h;
    dRaan_dt = r* sin(TA+ w(t-1)) *u(3)/(h* sin(i(t-1)));
    dw_dt = (-p*cos(TA)*u(1) + (p+r)*sin(TA)*u(2))/(h*e(t-1)) - r*sin(TA+2)*cos(i(t-1))*u(3)/(h*sin(i(t-1)));
    dM_dt = n + sqrt(1-e(t-1)^2) * ((p*cos(TA) - 2*r*e(t-1))*u(1) - (p+r)*sin(TA)*u(2))/(h*e(t-1));

    % Update orbit elements
    a(t) = a(t-1) + da_dt * dt;
    e(t) = e(t-1) + de_dt * dt;
    i(t) = i(t-1) + di_dt * dt;
    raan(t) = raan(t-1) + dRaan_dt * dt;
    w(t) = w(t-1) + dw_dt * dt;
    M(t) = M(t-1) + dM_dt * dt;
end

% Plotting the results
figure;
subplot(3,2,1);
plot((0:num_steps-1)*dt, a);
xlabel('Time (s)');
ylabel('Semi-major axis (m)');
title('Semi-major axis');

subplot(3,2,2);
plot((0:num_steps-1)*dt, e);
xlabel('Time (s)');
ylabel('Eccentricity');
title('Eccentricity');

subplot(3,2,3);
plot((0:num_steps-1)*dt, rad2deg(i));
xlabel('Time (s)');
ylabel('Inclination (deg)');
title('Inclination');

subplot(3,2,4);
plot((0:num_steps-1)*dt, rad2deg(raan));
xlabel('Time (s)');
ylabel('RAAN (deg)');
title('RAAN');

subplot(3,2,5);
plot((0:num_steps-1)*dt, rad2deg(w));
xlabel('Time (s)');
ylabel('Argument of Perigee (deg)');
title('Argument of Periapsis');

subplot(3,2,6);
plot((0:num_steps-1)*dt, rad2deg(M));
xlabel('Time (s)');
ylabel('Mean Anomaly (deg)');
title('Mean Anomaly');