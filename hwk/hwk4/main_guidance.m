% Main script
function main_guidance
    % Constants
    init;
    global mu acc
    mu = 3.986004415e5; % Earth's gravitational parameter [km^3/s^2]
    acc = 0.03; % Thrust acceleration [km/s^2]
    
    % Initial state
    x0 = [8276; 5612; 5; -3.142; 4.672; 0];
    tf = 23;
    % Target orbital elements
    target.a = 10000; % km
    target.i = 0; % deg
    target.e = 0.0001;
    target.omega = 35; % deg
    target.p = PfromAe(target.a,target.e); % semi-latus rectum
    % target.v_peri = [VfromRA(mu,norm(target.r_peri),target.a); 0; 0] ; % velocity at periapsis
    % target.r_peri = [0; 0; Rp(target.a,target.e)]; % radius at periapsis in lvlh frame
    target.h = [0; sqrt(mu*target.p); 0];      % angular momentum vector

    % compute periapsis information
    % First compute vectors in perifocal frame (PQW)
    target.rp_pqw = [target.a*(1-target.e); 0; 0];
    target.vp_pqw = [0; sqrt(mu/target.p)*(1+target.e); 0];
    
    % Create rotation matrix from perifocal to ECI
    % For zero inclination, this is just the argument of perigee rotation
    R_pqw2eci = [cosd(target.omega) -sind(target.omega) 0;
                 sind(target.omega) cosd(target.omega)  0;
                 0              0               1];
    
    % Transform to ECI frame
    target.rp_eci = R_pqw2eci * target.rp_pqw;
    target.vp_eci = R_pqw2eci * target.vp_pqw;
    
    % Initial guidance parameters (in LVLH frame)
    A = [1; 0; 0];
    B = [0; 0; 0];
    u0 = [A; B; tf]; % [Ax Ay Az Bx By Bz tf]
    lambda0 = [B; A+B*tf; tf]

    % Guidance loop parameters
    dt_guidance = 1.0; % Guidance update interval [s]
    max_iterations = 20;
    tolerance = 1e-3;
    
    % Run predictor-corrector guidance
    [final_state, guidance_history] = run_guidance(x0, u0, target, acc, dt_guidance, max_iterations, tolerance, lambda0);
end

function [final_state, guidance_history] = run_guidance(x0, u0, target, acc, dt_guidance, max_iterations, tolerance, lambda)
    t = 0;
    x = x0;
    u = u0;
    guidance_history.t = [];
    guidance_history.u = [];
    guidance_history.errors = [];
    
    while t < 23 % Continue until burn time
        for i = 1:max_iterations
            % Predict trajectory and compute errors
            [errors, predicted_state] = compute_errors(x, u, target, acc);
            
            % Check if converged
            if norm(errors) < tolerance
                break;
            end
            
            % Compute Jacobian using forward finite differencing
            J = compute_jacobian(x, u, errors, target, acc);
            
            % Update guidance parameters using Newton-Raphson
            delta_lambda = -pinv(J) * errors;
            deltaA = delta_lambda(4:6) - delta_lambda(1:3)*delta_lambda(7);
            delta_u = [deltaA; delta_lambda(1:3); delta_lambda(7)];
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
        end
        t = t + dt_guidance;
        [~, x_prop] = ode45(@(t,x) dynamics(t, x, u, acc), [t t+dt_guidance], x);
        x = x_prop(end,:)';
    end
    
    final_state = x;
end

function T_lvlh2eci = get_lvlh2eci_matrix(r, v)
    % Compute LVLH to ECI transformation matrix
    % z-axis points along -r (radial)
    % y-axis points along -h (orbit angular momentum)
    % x-axis completes right-handed system
    
    h = cross(r, v);
    z = r/norm(r);
    y = h/norm(h);
    x = cross(y, z);
    
    % Ensure orthonormal basis
    x = x/norm(x);
    y = y/norm(y);
    z = z/norm(z);
    
    % Form transformation matrix
    T_lvlh2eci = [x y z];
end

function dxdt = dynamics(t, x, u, acc)
    % Extract position and velocity
    global mu
    r = x(1:3);
    v = x(4:6);
    
    % Get transformation matrix from LVLH to ECI
    T_lvlh2eci = get_lvlh2eci_matrix(r, v);
    
    % Compute thrust direction in LVLH frame using bilinear tangent law
    A = u(1:3);
    B = u(4:6);
    lambda_lvlh = A + B*t;
    u_thrust_lvlh = lambda_lvlh/norm(lambda_lvlh);
    
    % Transform thrust direction to ECI frame
    u_thrust_eci = T_lvlh2eci * u_thrust_lvlh;
    
    % Compute accelerations in ECI frame
    g = -mu*r/norm(r)^3; % Gravity
    a_thrust = acc*u_thrust_eci; % Thrust
    
    % State derivative in ECI frame
    dxdt = [v; g + a_thrust];
end

function [errors, final_state] = compute_errors(x0, u, target, acc)
    % Simulate trajectory to final time
    global mu
    [~, x] = ode45(@(t,x) dynamics(t, x, u, acc), [0 u(end)], x0);
    final_state = x(end,:)';
    
    % Extract final position and velocity in ECI
    rf = final_state(1:3);
    vf = final_state(4:6);
    vp_eci = target.vp_eci;
    rp_eci = target.rp_eci;
    % Get transformation matrices at final state
    T_lvlh2eci = get_lvlh2eci_matrix(rf, vf);
    T_eci2lvlh = T_lvlh2eci';
    
    % Transform target vectors to LVLH for true anomaly computation

    rp_lvlh = T_eci2lvlh * target.rp_eci;
    rf_lvlh = T_eci2lvlh * rf;
    vp_lvlh = T_eci2lvlh * vp_eci;

    nu = atan(dot(rf,vp_eci/norm(vp_eci))/dot(rf,rp_eci/norm(rp_eci)));
    % Compute required position and velocity
   % p = norm(vp_eci)^2 * norm(rp_eci)^2 / mu;
    r_req = target.p/(1 + target.e * cos(nu));
    
    % Compute required velocity in ECI frame
    v_req_eci = compute_target_velocity(nu, rp_eci, vp_eci, target.p, mu);
    v_req_lvlh = T_eci2lvlh * v_req_eci;

    vf_lvlh = T_eci2lvlh * vf;
    % Compute angular momentum direction of target orbit
    h_req = cross(rp_eci, vp_eci)/norm(cross(rp_eci, vp_eci));
    
    % Compute lambda in LVLH frame
    lambda_lvlh = u(1:3) + u(4:6)*u(end);
    
    % Transform lambda to ECI for gravity alignment check
    lambda_eci = T_lvlh2eci * lambda_lvlh;
    
    % Compute gravity at final point in ECI frame
    g = -mu*rf/norm(rf)^3;
    g_lvlh = T_eci2lvlh * g;
    % Transform final velocity to LVLH for B·v check
    vf_lvlh = T_eci2lvlh * vf;
    
    % Compute errors
    errors = zeros(7,1);
    errors(1:3) = v_req_lvlh - vf;                          % Velocity error in ECI
    errors(4) = r_req - norm(rf);                          % Scalar radius error
    errors(5) = dot(h_req, rf);                            % Orbit plane alignment in ECI
    errors(6) = dot(lambda_lvlh, g_lvlh) - dot(u(4:6), vf_lvlh); % Primer vector condition
    errors(7) = sqrt(dot(lambda_lvlh, lambda_lvlh)) - 1;   % Lambda magnitude in LVLH
end

function J = compute_jacobian(x0, u, base_error, target, acc)
    % Extracts parameters
    A = u(1:3);
    B = u(4:6);
    tf = u(7);
    
    % Baseline control calculations
    u_baseline = u;
    
    % Initialize Jacobian matrix
    n_params = length(u);  % 7 parameters: A(3), B(3), tf
    n_errors = length(base_error);
    J = zeros(n_errors, n_params);
    
    % Perturbation size
    epsilonA = 0.1;
    epsilonB = 0.01;
    epsilonT = 0.1;
    
    % Compute Jacobian columns for B parameters
    for i = 1:3
        B_perturbed = B;
        B_perturbed(i) = B(i) + epsilonB;
        
        % Compute perturbed control and lambda
        u_perturbed = [A; B_perturbed; tf];
        
        % Simulate with perturbed parameters
        [~, x_final] = ode45(@(t, x) dynamics(t, x, u_perturbed, acc), [0 tf], x0);

        error_perturbed = compute_errors(x_final(end,:)', u_perturbed, target, acc);
        
        % Compute Jacobian column
        J(:,i) = (error_perturbed - base_error) / epsilonB;  
    end

    % Compute Jacobian columns for A parameters
    for i = 1:3
        A_perturbed = A;
        A_perturbed(i) = A(i) + epsilonA;
        
        % Compute perturbed control and lambda
        u_perturbed = [A_perturbed; B; tf];
        
        % Simulate with perturbed parameters
        [~, x_final] = ode45(@(t, x) dynamics(t, x, u_perturbed, acc), [0 tf], x0);

        error_perturbed = compute_errors(x_final(end,:)', u_perturbed, target, acc);
        
        % Compute Jacobian column
        J(:,i+3) = (error_perturbed - base_error) / epsilonA;
    end
    
    % Compute Jacobian column for tf
    tf_perturbed = tf + epsilonT;
    u_perturbed = [A; B; tf_perturbed];
    lambda_perturbed = compute_lambda(A, B, tf_perturbed);
    
    % Simulate with perturbed parameters
    [~, x_final] = ode45(@(t, x) dynamics(t, x, u_perturbed, acc), [0 tf], x0);

    error_perturbed = compute_errors(x_final(end,:)', u_perturbed, target, acc);
    
    J(:,7) = (error_perturbed - base_error) / epsilonT;
end

% function [rp_eci, vp_eci] = get_target_vectors(target)
%     % Compute perigee position and velocity vectors for target orbit in ECI frame
%     global mu
% 
%     % First compute vectors in perifocal frame (PQW)
%     rp_pqw = [target.a*(1-target.e); 0; 0];
%     vp_pqw = [0; sqrt(mu/target.p)*(1+target.e); 0];
% 
%     % Create rotation matrix from perifocal to ECI
%     % For zero inclination, this is just the argument of perigee rotation
%     omega_rad = deg2rad(target.omega);
%     R_pqw2eci = [cos(omega_rad) -sin(omega_rad) 0;
%                  sin(omega_rad) cos(omega_rad)  0;
%                  0              0               1];
% 
%     % Transform to ECI frame
%     rp_eci = R_pqw2eci * rp_pqw;
%     vp_eci = R_pqw2eci * vp_pqw;
% end
% 

function v_req_eci = compute_target_velocity(nu, rp_eci, vp_eci, p, mu)
    % Compute required velocity at true anomaly nu in ECI frame
    sin_nu = sin(nu);
    cos_nu = cos(nu);
    v_req_eci = -sqrt(mu/p)*sin_nu*rp_eci/norm(rp_eci) + (1 - norm(rp_eci)*(1-cos_nu)/p)*vp_eci;
end

function lambda = compute_lambda(A, B, t)
    lambda = A + B*t;
end