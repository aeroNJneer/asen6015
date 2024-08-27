function [tout, xnl] = nlmodel(t, X0, mu, w, u) 
% Define nonlinear model as anonymous function and set options
    xdot = @(x, u, w) [x(2);
        -mu*x(1)/sqrt((x(1)^2+x(3)^2)^3) + u(1) + w(1);
        x(4);
        -mu*x(3)/sqrt((x(1)^2+x(3)^2)^3) + u(2) + w(2)];
    options = odeset('RelTol', 1e-12);

    % Define wrapper function for use in ODE45
    xdotwrap = @(t, x) xdot(x, [0; 0], w);
    % Compute total state using ODE45
    [tout, xnl] = ode45(xdotwrap,t, X0, options);
end