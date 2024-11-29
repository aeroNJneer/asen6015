function [t, state, control] = optimizeLowThrustTransfer()
    % Optimize low-thrust transfer from LEO to GEO
    % Uses indirect optimization method with costate variables
    
    % Constants
    mu = 398600.4418;        % Earth's gravitational parameter [km^3/s^2]
    g0 = 9.81;              % Standard gravity [m/s^2]
    Isp = 3000;             % Specific impulse [s]
    T = 0.5;                % Thrust [N]
    m0 = 1500;              % Initial mass [kg]
    mdot = T/(g0*Isp);      % Mass flow rate [kg/s]
    
    % Initial and final orbits
    r_LEO = 6678;           % LEO radius [km]
    r_GEO = 42164;          % GEO radius [km]
    
    % Initial guess for time of flight
    tf_guess = 30*24*3600;  % 30 days in seconds
    
    % Initial state vector [r, theta, vr, vtheta, m]
    x0 = [r_LEO; 0; 0; sqrt(mu/r_LEO); m0];
    
    % Target state (excluding theta)
    xf_target = [r_GEO; 0; sqrt(mu/r_GEO); NaN];
    
    % Initial guess for costate variables
    lambda0_guess = zeros(5,1);
    
    % Optimization options
    options = optimset('Display', 'iter', 'MaxFunEvals', 1000);
    
    % Optimize costate initial conditions and final time
    vars_guess = [lambda0_guess; tf_guess];
    [vars_opt, fval] = fmincon(@(vars) costFunction(vars, x0, xf_target, mu, T, mdot), ...
        vars_guess, [], [], [], [], [], [], ...
        @(vars) constraints(vars, x0, xf_target, mu, T, mdot), options);
    
    % Extract optimal results
    lambda0_opt = vars_opt(1:5);
    tf_opt = vars_opt(end);
    
    % Integrate optimal trajectory
    tspan = linspace(0, tf_opt, 1000);
    [t, state] = ode45(@(t,x) dynamics(t, x, mu, T, mdot), tspan, [x0; lambda0_opt]);
    
    % Calculate optimal control
    control = calculateControl(state, T);
end

function dxdt = dynamics(t, x, mu, T, mdot)
    % State and costate equations of motion
    r = x(1);
    theta = x(2);
    vr = x(3);
    vtheta = x(4);
    m = x(5);
    lambda = x(6:10);
    
    % Control law from Pontryagin's Maximum Principle
    [u_r, u_theta] = optimalControl(x, lambda, T);
    
    % State derivatives
    dr = vr;
    dtheta = vtheta/r;
    dvr = vtheta^2/r - mu/r^2 + T*u_r/m;
    dvtheta = -vr*vtheta/r + T*u_theta/m;
    dm = -mdot;
    
    % Costate derivatives
    dlambda = costateEquations(x, lambda, mu, T, u_r, u_theta);
    
    dxdt = [dr; dtheta; dvr; dvtheta; dm; dlambda];
end

function [u_r, u_theta] = optimalControl(x, lambda, T)
    % Calculate optimal control direction from primer vector
    lambda_v = [lambda(3); lambda(4)];
    norm_lambda_v = norm(lambda_v);
    
    if norm_lambda_v < 1e-10
        u_r = 0;
        u_theta = 0;
    else
        u_r = lambda_v(1)/norm_lambda_v;
        u_theta = lambda_v(2)/norm_lambda_v;
    end
end

function dlambda = costateEquations(x, lambda, mu, T, u_r, u_theta)
    % Costate differential equations
    r = x(1);
    vr = x(3);
    vtheta = x(4);
    m = x(5);
    
    dlambda = zeros(5,1);
    
    % Detailed costate equations (derived from Hamiltonian)
    dlambda(1) = lambda(2)*vtheta/r^2 - lambda(3)*vtheta^2/r^2 + ...
                 2*lambda(3)*mu/r^3 + lambda(4)*vr*vtheta/r^2;
    dlambda(2) = 0;
    dlambda(3) = -lambda(1) + lambda(4)*vtheta/r;
    dlambda(4) = -lambda(2)/r + 2*lambda(3)*vtheta/r - lambda(4)*vr/r;
    dlambda(5) = T*(u_r*lambda(3) + u_theta*lambda(4))/m^2;
    
end

function [c, ceq] = constraints(vars, x0, xf_target, mu, T, mdot)
    % Boundary constraints for optimization
    lambda0 = vars(1:5);
    tf = vars(end);
    
    % Integrate trajectory
    [~, state] = ode45(@(t,x) dynamics(t, x, mu, T, mdot), [0 tf], [x0; lambda0]);
    xf = state(end,1:5)';
    
    % Terminal constraints (matching target conditions)
    ceq = [xf(1) - xf_target(1);   % Final radius
           xf(3);                   % Final radial velocity
           xf(4) - xf_target(3)];   % Final tangential velocity
    
    c = [];  % No inequality constraints
end

function cost = costFunction(vars, x0, xf_target, mu, T, mdot)
    % Cost function (time of flight)
    tf = vars(end);
    cost = tf;
end

function control = calculateControl(state, T)
    % Calculate control history from state and costate
    n = size(state,1);
    control = zeros(n,2);
    
    for i = 1:n
        [u_r, u_theta] = optimalControl(state(i,:)', state(i,6:10)', T);
        control(i,:) = [u_r, u_theta];
    end
end