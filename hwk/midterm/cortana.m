clear
mu = 3.986004418e5;
a = 3728;           % semi-major axis of ISS
n = sqrt(mu/a^3);   %  mean motion of ISS
dt = 0.1;           % timestep 
tfinal = 60*20;
tsim = 0:dt:tfinal;   
tf = 1200;
dt = 0.1;
% function signatures

tsim = 0:0.1:1200;
%initial conditions
x0 = -1.0;
y0 = -0.25
xdot0 = 0.003;
ydot0 = -0.0005;
X0 = [x0; y0; 0;xdot0; ydot0;0];

% final conditions
xf = -0.25 ;
yf = 0. ;
xdotf = 0.0002 ;
ydotf = 0. ;
Xf = [xf; yf; 0;xdotf; ydotf;0];

r0 = [x0;y0;0];
v0 = [xdot0;ydot0;0];
rf = [xf;yf;0];
vf = [xdotf;ydotf;0];

dX = Xf - X0

    
    % Initial conditions
    state0 = [r0; v0];
    
    % Solve for the first Δv using STM
    A = cw_transfer_matrix(n, tf);
    delta_v1 = A \ ([rf; 0; 0; 0] - state0');
    v0_new = v0 + delta_v1(4:6)';
    state0_new = [r0; v0_new];
    
    % Propagate with first Δv
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
    [t1, state1] = ode45(@cw_dynamics, [0 tf], state0_new, options, n);
    
    % Calculate correction using STM
    pos_error = rf - state1(end, 1:3)';
    B = cw_transfer_matrix(n, tf - t1(end));
    delta_v2 = B \ [pos_error; 0; 0; 0];
    
    % Apply second impulse at t_correction
    t_correction = tf / 2;
    [t2, state2] = ode45(@cw_dynamics, [0 t_correction], state0_new, options, n);
    state2(end, 4:6) = state2(end, 4:6) + delta_v2(4:6)';
    
    % Propagate from correction to final time
    [t3, state3] = ode45(@cw_dynamics, [t_correction tf], state2(end, :)', options, n);
    
    % Combine results
    t = [t2; t3];
    state = [state2; state3];
    
    % Plot results
    figure;
    plot3(state(:,1), state(:,2), state(:,3), 'r', 'DisplayName', 'Adjusted Trajectory');
    hold on;
    plot3(r0(1), r0(2), r0(3), 'go', 'DisplayName', 'Start Point');
    plot3(rf(1), rf(2), rf(3), 'ro', 'DisplayName', 'End Point');
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Relative Motion with STM and Two Impulses');
    legend;
    grid on;
    hold off;

    % CW dynamics function
    function dstate_dt = cw_dynamics(~, state, n)
        x = state(1);
        y = state(2);
        z = state(3);
        vx = state(4);
        vy = state(5);
        vz = state(6);
        
        ax = 3 * n^2 * x + 2 * n * vy;
        ay = -2 * n * vx;
        az = -n^2 * z;
        
        dstate_dt = [vx; vy; vz; ax; ay; az];
    end

function A = cw_transfer_matrix(n, tf)
    A = [4 - 3 * cos(n * tf), 0, 0, sin(n * tf) / n, 2 * (1 - cos(n * tf)) / n, 0;
         6 * (sin(n * tf) - n * tf), 1, 0, -2 * (1 - cos(n * tf)) / n, (4 * sin(n * tf) - 3 * n * tf) / n, 0;
         0, 0, cos(n * tf), 0, 0, sin(n * tf) / n;
         3 * n * sin(n * tf), 0, 0, cos(n * tf), 2 * sin(n * tf), 0;
         6 * n * (cos(n * tf) - 1), 0, 0, -2 * sin(n * tf), 4 * cos(n * tf) - 3, 0;
         0, 0, -n * sin(n * tf), 0, 0, cos(n * tf)];
end
