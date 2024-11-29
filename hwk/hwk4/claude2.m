% Main script
function main_guidance
    % Constants
    global mu = 398600; % Earth's gravitational parameter [km^3/s^2]
    acc = 0.03; % Thrust acceleration [km/s^2]
    
    % Initial state
    x0 = [8276; 5612; 5; -3.142; 4.672; 0];
    
    % Target orbital elements
    target.a = 10000; % km
    target.i = 0; % deg
    target.e = 0.0001;
    target.omega = 35; % deg
    
    % Initial guidance parameters
    u0 = [1; 0; 0; 0; 0; 0; 23]; % [Ax Ay Az Bx By Bz tf]
    
    % Guidance loop parameters
    dt_guidance = 1.0; % Guidance update interval [s]
    max_iterations = 20;
    tolerance = 1e-6;
    
    % Run predictor-corrector guidance
    [final_state, guidance_history] = run_guidance(x0, u0, target, acc, dt_guidance, max_iterations, tolerance);
end

function [final_state, guidance_history] = run_guidance(x0, u0, target, acc, dt_guidance, max_iterations, tolerance)
    t = 0;
    x = x0;
    u = u0;
    guidance_history.t = [];
    guidance_history.u = [];
    guidance_history.errors = [];
    
    while t < u(end) % Continue until burn time
        % Predict trajectory and compute errors
        [errors, predicted_state] = compute_errors(x, u, target, acc);
        
        % Check if converged
        if norm(errors) < tolerance
            break;
        end
        
        % Compute Jacobian using forward finite differencing
        J = compute_jacobian(x, u, target, acc);
        
        % Update guidance parameters using Newton-Raphson
        delta_u = -pinv(J) * errors;
        
        % Apply step size control for better convergence
        alpha = 1.0;
        while norm(compute_errors(x, u + alpha*delta_u, target, acc)) > norm(errors) && alpha > 1e-3
            alpha = alpha * 0.5;
        end
        
        u = u + alpha*delta_u;
        
        % Store history
        guidance_history.t = [guidance_history.t; t];
        guidance_history.u = [guidance_history.u; u'];
        guidance_history.errors = [guidance_history.errors; errors'];
        
        % Propagate state
        [~, x_prop] = ode45(@(t,x) dynamics(t, x, u, acc), [t t+dt_guidance], x);
        x = x_prop(end,:)';
        t = t + dt_guidance;
    end
    
    final_state = x;
end

function J = compute_jacobian(x, u, target, acc)
    % Compute Jacobian using forward finite differencing
    epsilon = 1e-6;
    n_errors = 7;
    n_controls = length(u);
    J = zeros(n_errors, n_controls);
    
    base_errors = compute_errors(x, u, target, acc);
    
    for i = 1:n_controls
        u_perturbed = u;
        u_perturbed(i) = u_perturbed(i) + epsilon;
        perturbed_errors = compute_errors(x, u_perturbed, target, acc);
        J(:,i) = (perturbed_errors - base_errors) / epsilon;
    end
end

function [errors, final_state] = compute_errors(x0, u, target, acc)
    % Simulate trajectory to final time
    [~, x] = ode45(@(t,x) dynamics(t, x, u, acc), [0 u(end)], x0);
    final_state = x(end,:)';
    
    % Extract final position and velocity
    rf = final_state(1:3);
    vf = final_state(4:6);
    
    % Compute target orbit parameters
    [rp, vp] = get_target_vectors(target);
    
    % Compute true anomaly
    nu = atan2(dot(rf, vp), dot(rf, rp));
    
    % Compute required position and velocity
    p = norm(vp)^2 * norm(rp)^2 / 398600;
    r_req = p/(1 + target.e * cos(nu));
    v_req = compute_target_velocity(nu, rp, vp, p, 398600);
    
    % Compute angular momentum direction of target orbit
    h_req = cross(rp, vp)/norm(cross(rp, vp));
    
    % Compute lambda
    lambda = u(1:3) + u(4:6)*u(end);
    
    % Compute gravity at final point
    g = -398600*rf/norm(rf)^3;
    
    % Compute errors
    errors = zeros(7,1);
    errors(1:3) = v_req - vf;
    errors(4) = r_req - norm(rf);
    errors(5) = dot(h_req, rf);
    errors(6) = dot(lambda, g) - dot(u(4:6), vf);
    errors(7) = sqrt(dot(lambda, lambda)) - 1;
end

function dxdt = dynamics(t, x, u, acc)
    % Extract position and velocity
    r = x(1:3);
    v = x(4:6);
    
    % Compute thrust direction using bilinear tangent law
    A = u(1:3);
    B = u(4:6);
    lambda = A + B*t;
    u_thrust = lambda/norm(lambda);
    
    % Compute acceleration
    g = -398600*r/norm(r)^3; % Gravity
    a_thrust = acc*u_thrust; % Thrust
    
    % State derivative
    dxdt = [v; g + a_thrust];
end

function [rp, vp] = get_target_vectors(target)
    % Compute perigee position and velocity vectors for target orbit
    p = target.a * (1 - target.e^2);
    rp = [target.a*(1-target.e); 0; 0];
    vp = [0; sqrt(398600/p)*(1+target.e); 0];
    
    % Apply rotation for argument of perigee
    omega_rad = deg2rad(target.omega);
    R = [cos(omega_rad) -sin(omega_rad) 0;
         sin(omega_rad) cos(omega_rad) 0;
         0 0 1];
    rp = R*rp;
    vp = R*vp;
end

function v_req = compute_target_velocity(nu, rp, vp, p, mu)
    % Compute required velocity at true anomaly nu
    sin_nu = sin(nu);
    cos_nu = cos(nu);
    v_req = -sqrt(mu/p)*sin_nu*rp + (1 - norm(rp)/p*(1-cos_nu))*vp;
end