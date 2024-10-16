function dragon_approach
    % Constants
    a = 6728; % Semi-major axis (km)
    mu = 398600.4418; % Gravitational parameter (km^3/s^2)
    n = sqrt(mu / a^3); % Mean motion (rad/s)
    
    % Initial conditions
    x0 = 1.0; % km
    y0 = 0.25; % km
    vx0 = 3 / 1000; % km/s (3 m/s)
    vy0 = 0.5 / 1000; % km/s (0.5 m/s)
    initial_state = [x0; y0; vx0; vy0];
    
    % Final conditions
    xf = 0.25; % km
    yf = 0; % km
    vxf = 0.2 / 1000; % km/s (0.2 m/s)
    vyf = 0; % km/s
    
    % Time of flight
    tf = 1200; % seconds
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
    
    % Feedback gain matrix K (example gains, adjust as needed)
    K = [-0.001, 0, -0.000005, 0;
          0, -0.001, 0, -0.000005] * .000001;
    
    % Define desired state
    desired_state = [xf; yf; vxf; vyf];
    
    % Solve ODE with feedback control
    [t, state] = ode45(@cw_dynamics_with_feedback, [0 tf], initial_state, options, n, K, desired_state);
    
    % Plot results
    figure;
    plot(state(:,1), state(:,2), '-b', 'DisplayName', 'Trajectory with Linear Feedback');
    hold on;
    plot(x0, y0, 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
    plot(xf, yf, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'End');
    xlabel('X (km)');
    ylabel('Y (km)');
    title('Trajectory with Linear Feedback Control');
    legend;
    grid on;
    hold off;

    % CW dynamics function with linear feedback control
    function dstate_dt = cw_dynamics_with_feedback(~, state, n, K, desired_state)
        x = state(1);
        y = state(2);
        vx = state(3);
        vy = state(4);
        
        % Calculate control input
        u = K * (desired_state - state);
        
        ax = 3 * n^2 * x + 2 * n * vy + u(1);
        ay = -2 * n * vx + u(2);
        
        dstate_dt = [vx; vy; ax; ay];
    end
end
