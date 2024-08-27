close all
clear
% constants
    ymax = 10;     % y-limit for NEES/NIS plots
    alpha = 0.05;  % 95% confidence
    mu = 398600; % gravitational constant
    omegaE = 2*pi/86400;  % angular velocity of earth
    rE = 6378; % radius of earth
    dt = 10;
    t = 0:dt:14000;
    stationNums = 1:12;
% initial conditions
    x0 = 6678;
    y0 = 0;
    r0 = sqrt(x0^2+y0^2);
    x0dot = 0;
    y0dot = r0 * sqrt(mu/r0^3);
    X0 = [x0, x0dot, y0, y0dot]';
    omega = sqrt(mu/r0^3);          % angular velocity of satellite
    A0 = getAtilde(X0);
    B = [0 0; 1 0; 0 0; 0 1]; 
    G = dt * B; 
% end initial conditions

% nonlinear Model
    xdot = "@(x, u, w) [x(2); -mu*x(1)/sqrt((x(1)^2+x(3)^2)^3) + u(1) + w(1); x(4); -mu*x(3)/sqrt((x(1)^2+x(3)^2)^3) + u(2) + w(2)]";
% end nonlinear Model
    kf = KalmanFilter();
    kf.modelState(t, xdot, X0, mu, [0 0]);
