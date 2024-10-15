close all
    clear
    mu = 3.986004418e5;
a = 6728;           % semi-major axis of ISS
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
X0 = [x0; y0; xdot0; ydot0];

% final conditions
xf = -0.25 ;
yf = 0. ;
xdotf = 0.0002 ;
ydotf = 0. ;
Xf = [xf; yf; xdotf; ydotf];

r0 = [x0;y0];
v0 = [xdot0;ydot0];
rf = [xf;yf];
vf = [xdotf;ydotf];

dX = Xf - X0
    A = cw_transfer_matrix(n, tf);
    B = cw_transfer_matrix(n,tf/2)
    % Solve for the first Δv
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
    tf = 1200; % Total time of flight (s)
    
     % Solve for the first Δv using STM

    [t1,X1] =  ode45(@cw_dynamics, 0:dt:tf, X0, options, n);
    vgo1 = Xf - X1(end,:)'
    deltav1= A\vgo1
    
    % simulate with dv1
    X1 = [r0; deltav1(3:4)+v0]
    [t2, X2] = ode45(@cw_dynamics, 0:dt:tf, X1, options, n);
   
    % now that velocity is on right path, look at position errors
    poserror = rf - X2(end,1:2)'

    deltav2 = B \ ([poserror*n; vf])
    
    % Propagate with first Δv
    [t3, X3] = ode45(@cw_dynamics, 0:dt:tf, X1, options, n);

    X3(end, 3:4) = X3(end, 3:4) + deltav2(3:4)';

    % Apply second impulse at t_correction
    t_correction = tf / 2;
    
    % Propagate from correction to final time
    [t4, X4] = ode45(@cw_dynamics, t_correction:dt:tf, X3(end, :)', options, n);
    
    % Combine results
    t = [t3; t4];
    state = [X3; X4];
    
    % Plot results
    figure(1);
    plot(X3(:,1), X3(:,2), 'r', 'DisplayName', 'Adjusted Trajectory');
    hold on;
    plot(r0(1), r0(2), 'go', 'DisplayName', 'Start Point');
    plot(rf(1), rf(2), 'ro', 'DisplayName', 'End Point');
    xlabel('X (km)');
    ylabel('Y (km)');
    title('Relative Motion with STM and Two Impulses');
    legend;
    grid on;
    hold off;

    % CW dynamics function
    function dstate_dt = cw_dynamics(~, state, n)
        x = state(1);
        y = state(2);
        vx = state(3);
        vy = state(4);
        
        ax = 3 * n^2 * x + 2 * n * vy;
        ay = -2 * n * vx;
        
        dstate_dt = [vx; vy; ax; ay];
    end


function A = cw_transfer_matrix(n, t)
    Phirr = @(t) [4-3*cos(n*t), 0;
                6*(sin(n*t)-n*t), 1];
    Phirv = @(t) [sin(n*t)/n, 2*(1-cos(n*t))/n;
                2*(cos(n*t) -1)/n, (4*sin(n*t) - 3*n*t)/n];
    Phivr = @(t) [3*n*sin(n*t), 0;
                6*n*(cos(n*t)-1), 0];
    Phivv = @(t) [cos(n*t), 2*sin(n*t);
                -2*sin(n*t), 4*cos(n*t)-3];

    A = [Phirr(t) Phirv(t);
            Phivr(t) Phivv(t)];
end
