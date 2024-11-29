% Main script for orbital transfer using predictor-corrector guidance
function main_orbital_transfer()
init;    
% Constants
    global mu thrust_acc
    mu = 398600.4418; % Earth's gravitational parameter (km^3/s^2)
    thrust_acc = 0.03; % 30 m/s^2 converted to km/s^2
    
    % Initial state
    x0 = [8276; 5612; 5; -3.142; 4.672; 0];
    
    % Target orbital elements
    target.a = 10000;  % km
    target.i = 0;      % deg
    target.e = 0.0001;
    target.omega = 35; % deg
    target.r_peri = [0; 0; Rp(target.a,target.e)]; % radius at periapsis in lvlh frame
    target.p = PfromAe(target.a,target.e); % semi-latus rectum
    target.v_peri = [VfromRA(mu,norm(target.r_peri),target.a); 0; 0] ; % velocity at periapsis
    target.h = 
    % Initial guidance parameters
    A = [1; 0; 0];
    B = [0; 0; 0];
    tf = 23; % Initial burn time guess
    
    % Initialize control unit vector
    u = [A; B; tf];
    lambda = [B; A+B*tf; tf];
    primer = lambda(4:6)/norm(4:6);
    % Guidance loop parameters
    max_iterations = 40;
    tolerance = 1e-6;
    
    % Main guidance loop
    for iter = 1:max_iterations
        % Predict trajectory with current parameters
        [t, tf_state] = ode45(@(t, x) dynamics(x, t, primer), [0 tf], x0);
      
        % Compute error vector
        error = compute_error_vector(tf_state(end,:)', target, lambda);
        
        if norm(error) < tolerance
            break;
        end
        
        % Compute Jacobian using forward differences
        J = compute_jacobian(x0, u, error, target);
        
        % Update control vector
        delta_lambda = -pinv(J) * error;
        
        % Apply step size control for better convergence
        alpha = 1.0;
        while norm(delta_lambda) > 0.5
            alpha = alpha * 0.5;
            delta_lambda = alpha * delta_lambda;
        end
        lambda = lambda + delta_lambda;
        a = lambda(4:6) - lambda(1:3)*lambda(7);     % A = B*t - lambda
        u = [a; lambda(4:6); lambda(7)];
    end
    
    % Final simulation with converged parameters
   % [states, times] = simulate_trajectory(x0, u);
    [times, states] = ode45(@(t, x) dynamics(x, t, u), [0 tf], x0);

    plot_results(states, times, target);
end

% function [states, times] = simulate_trajectory(x0, u)
%     % Simulate trajectory using RK4 integration
%     t0 = 0;
%     tf = u(end);
%     dt = 0.1; % Time step
% 
%     times = t0:dt:tf;
%     states = zeros(6, length(times));
%     states(:,1) = x0;
% 
%     for i = 1:length(times)-1
%         k1 = dynamics(states(:,i), times(i), u);
%         k2 = dynamics(states(:,i) + dt/2*k1, times(i) + dt/2, u);
%         k3 = dynamics(states(:,i) + dt/2*k2, times(i) + dt/2, u);
%         k4 = dynamics(states(:,i) + dt*k3, times(i) + dt, u);
% 
%         states(:,i+1) = states(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
%     end
% end

function dxdt = dynamics(x, t, u)
    global mu thrust_acc
    
    % Extract position and velocity
    r = x(1:3);
    v = x(4:6);
    % A = u(1:3);
    % B = u(4:6);
    tf = u(7);
    % Compute gravity in ECI
    r_mag = norm(r);
    g_eci = -mu * r / r_mag^3;
    
    % Get transformation matrix
    ECI_to_LVLH = get_ECI_to_LVLH(r, v);
    
    % Compute normalized control in LVLH
   % u_lvlh = (A + B*t)/norm(A + B*t);
    
    % Transform normalized control from LVLH to ECI
    u_eci = ECI_to_LVLH' * u;
    
    % Compute state derivatives
    dxdt = zeros(6,1);
    dxdt(1:3) = v;
    dxdt(4:6) = g_eci + thrust_acc * u_eci;
end

function error = compute_error_vector(final_state, target, u)
    global mu
    
    % Extract final position and velocity in ECI
    rf = final_state(1:3);
    vf = final_state(4:6);
    
    % Get transformation to LVLH frame at final state
    ECI_to_LVLH = get_ECI_to_LVLH(rf, vf);
    
    % Transform final velocity to LVLH
    vf_lvlh = ECI_to_LVLH * vf;
    
    % Compute current orbital elements
    [a, e, i, omega, ~, nu] = rv2oe(rf, vf, mu, target);
    
    % Compute target orbit periapsis state
    rp = norm(target.r_peri);
    vp = norm(target.v_peri);
    
    % Compute required velocity at current position (in LVLH)
    vR_lvlh = compute_required_velocity(rf, target, nu);
    
    % Compute angular momentum direction of target orbit in ECI frame
    h_target = [sind(target.i)*sind(target.omega);
                -sind(target.i)*cosd(target.omega);
                cosd(target.i)];
    
    % Compute lambda in LVLH
    lambda_lvlh = u;
    
    % Compute gravity at final state and transform to LVLH
    gf_eci = -mu * rf / norm(rf)^3;
    gf_lvlh = ECI_to_LVLH * gf_eci;
    
    % Assemble error vector according to equations 10-16
    error = zeros(7,1);
    error(1:3) = vR_lvlh - vf_lvlh;  % Velocity errors in LVLH frame
    error(4) = target.a * (1-target.e^2)/(1+target.e*cosd(nu)) - norm(rf);
    error(5) = dot(h_target, rf);
    error(6) = dot(lambda_lvlh, gf_lvlh) - dot(u(4:6), vf_lvlh);
    error(7) = sqrt(dot(lambda_lvlh, lambda_lvlh)) - 1;
end


function [v_req] = compute_required_velocity(r, target, nu)
    global mu
    nu = deg2rad(nu);
    % Compute velocity required for target orbit at current position
    h = sqrt(mu * target.p);
    r_mag = norm(r);

    % compute required position
    r_req = target.p/(1-target.e*cos(nu))

    % compute required velocity
   % v_req = -target.r_peri*sqrt(mu/target.p) * sin(nu)/norm(target.r_peri) + (1 - norm(target.r_peri)/target.p * (1-cos(nu)) ) * target.v_peri;
    
    % % Compute velocity components in orbital frame
    v_r = mu/h * target.e * sin(nu);
    v_theta = h/r_mag;
    V = [v_r v_theta 0];
    % Transform to inertial frame (simplified for this implementation)
    vR1 = -v_r * sin(nu) + v_theta * cos(nu);
    vR2 = v_r * cos(nu) + v_theta * sin(nu);
    vR3 = 0; % Assuming planar transfer
    v_req = [vR1; vR2; vR3];
end

function [a, e, i, omega, Omega, nu] = rv2oe(r, v, mu, target)
    % Convert state vectors to orbital elements
    h = cross(r, v);
    n = cross([0;0;1], h);
    
    r_mag = norm(r);
    v_mag = norm(v);
    h_mag = norm(h);
    n_mag = norm(n);
    
    % Semi-major axis
    a = 1/(2/r_mag - v_mag^2/mu);
    
    % Eccentricity vector and magnitude
    e_vec = ((v_mag^2 - mu/r_mag)*r - dot(r,v)*v)/mu;
    e = norm(e_vec);
    
    % Inclination
    i = acosd(h(3)/h_mag);
    
    % Right ascension of ascending node
    Omega = acosd(n(1)/n_mag);
    if n(2) < 0
        Omega = 360 - Omega;
    end
    
    % Argument of periapsis
    omega = acosd(dot(n,e_vec)/(n_mag*e));
    if e_vec(3) < 0
        omega = 360 - omega;
    end
    
    % True anomaly of intercept
    nu = atan2(dot(r,target.v_peri),dot(r,target.r_peri));
    % if dot(r,v) < 0
    %     nu = 360 - nu;
    % end
end

function ECI_to_LVLH = get_ECI_to_LVLH(r, v)
    % Compute LVLH frame unit vectors
    % R_hat (z) points radially outward
    % H_hat (y) points along angular momentum
    % T_hat (x) completes right-handed system
    
    r_mag = norm(r);
    R_hat = r / r_mag;
    
    h = cross(r, v);
    h_mag = norm(h);
    H_hat = h / h_mag;
    
    T_hat = cross(H_hat, R_hat);
    T_hat = T_hat / norm(T_hat);
    
    % Assemble the DCM
    % Each column is a unit vector of LVLH frame expressed in ECI
    ECI_to_LVLH = [T_hat, H_hat, R_hat]';
end


function J = compute_jacobian(x0, lambda, base_error, target)
    % u_params contains [A; B; tf]
    % Extracts parameters
    B = lambda(1:3);
    lambdav = lambda(4:6);
    tf = lambda(7);
    
    % Baseline control calculations
    lambda_baseline = lambda;
    
    % Initialize Jacobian matrix
    n_params = length(lambda);  % 7 parameters: A(3), B(3), tf
    n_errors = length(base_error);
    J = zeros(n_errors, n_params);
    
    % Perturbation size
    epsilonA = 0.1;
    epsilonB = 0.01;
    epsilonT = 0.1;
    
    % Compute Jacobian columns for A parameters
    for i = 1:3
        B_perturbed = B;
        B_perturbed(i) = B(i) + epsilonB;
        
        % Compute perturbed control and lambda
        u_perturbed = [A; B_perturbed; tf];
        lambda_perturbed = compute_lambda(A, B_perturbed, tf);
        
        % Simulate with perturbed parameters
        %x_final = simulate_trajectory(x0, u_perturbed);
        [~, x_final] = ode45(@(t, x) dynamics(x, t, u_perturbed), [0 tf], x0);

        error_perturbed = compute_error_vector(x_final(end,:)', target, u_perturbed);
        
        % Compute Jacobian column
        J(:,i) = (error_perturbed - base_error) / lambda_perturbed(i);
    end

    % Compute Jacobian columns for A parameters
    for i = 1:3
        A_perturbed = A;
        A_perturbed(i) = A(i) + epsilonA;
        
        % Compute perturbed control and lambda
        u_perturbed = [A_perturbed; B; tf];
        lambda_perturbed = compute_lambda(A_perturbed, B, tf);
        
        % Simulate with perturbed parameters
        %x_final = simulate_trajectory(x0, u_perturbed);
        [~, x_final] = ode45(@(t, x) dynamics(x, t, u_perturbed), [0 tf], x0);

        error_perturbed = compute_error_vector(x_final(end,:)', target, u_perturbed);
        
        % Compute Jacobian column
        J(:,i+3) = (error_perturbed - base_error) / lambda_perturbed(i+3);
    end
    
    % Compute Jacobian column for tf
    tf_perturbed = tf + epsilonT;
    u_perturbed = [A; B; tf_perturbed];
    lambda_perturbed = compute_lambda(A, B, tf_perturbed);
    
   % x_final = simulate_trajectory(x0, u_perturbed);
    [~, x_final] = ode45(@(t, x) dynamics(x, t, u_perturbed), [0 tf], x0);

    error_perturbed = compute_error_vector(x_final(end,:)', target, u_perturbed);
    
    J(:,7) = (error_perturbed - base_error) / lambda_perturbed(7);
end

function u = compute_normalized_control(A, B, t)
    % Compute normalized thrust direction u = (A + Bt)/|A + Bt|
    unnormalized = A + B*t;
    u = unnormalized / norm(unnormalized);
end

function lambda = compute_lambda(A, B, t)
    % Compute primer vector lambda = [B; A + Bt]
    lambda = [B; A + B*t; t];
end

function plot_results(states, times, target)
    figure('Position', [100 100 1200 800]);
    
    % 3D trajectory plot
    subplot(2,2,1)
    plot3(states(:,1), states(:,2), states(:,3), 'b-', 'LineWidth', 2)
    hold on
    plot3(states(1,1), states(1,2), states(1,3), 'go', 'MarkerSize', 10)
    plot3(states(end,1), states(end,2), states(end,3), 'ro', 'MarkerSize', 10)
    grid on
    xlabel('X (km)')
    ylabel('Y (km)')
    zlabel('Z (km)')
    title('Transfer Trajectory')
    legend('Trajectory', 'Start', 'End')
    
    % Orbital elements evolution
    subplot(2,2,2)
    [a_hist, e_hist, i_hist, omega_hist, ~, ~] = deal(zeros(1,length(times)));
    global mu
    for j = 1:length(times)
        [a_hist(j), e_hist(j), i_hist(j), omega_hist(j), ~, ~] = ...
            rv2oe(states(j,1:3), states(j,4:6), mu, target);
    end
    
    plot(times, a_hist, 'b-', times, e_hist*1e4, 'r-', ...
         times, i_hist, 'g-', times, omega_hist, 'm-')
    grid on
    xlabel('Time (s)')
    ylabel('Value')
    title('Orbital Elements Evolution')
    legend('Semi-major axis (km)', 'Eccentricity (x10^{-4})', ...
           'Inclination (deg)', 'Arg. of Periapsis (deg)')
    
    % Velocity magnitude
    subplot(2,2,3)
    v_mag = sqrt(sum(states(:,4:6).^2));
    plot(times, v_mag)
    grid on
    xlabel('Time (s)')
    ylabel('Velocity (km/s)')
    title('Velocity Magnitude')
    
    % Position magnitude
    subplot(2,2,4)
    r_mag = sqrt(sum(states(:,1:3).^2));
    plot(times, r_mag)
    grid on
    xlabel('Time (s)')
    ylabel('Radius (km)')
    title('Orbital Radius')
end