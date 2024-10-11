% ASEN 6015
% MIDTERM EXAM
% Nicole James
% 10/10/24

clear
clc
close all
mu = 98600.4418; 
a = 3728;           % semi-major axis of ISS
n = sqrt(mu/a^3);   %  mean motion of ISS
dt = 0.001;           % timestep 

t = 0:dt:10;    % 50 seconds
%initial conditions
x0 = -1.0;
y0 = -0.25
xdot0 = 3;
ydot0 = -0.5;

% xddot - 2*n*ydot - 3*n^2*x = ux
% yddot + 2*n*xdot = uy
% X = [x; xdot; y; ydot]
% u = [ux; uy]

X0 = [x0; xdot0; y0; ydot0];
u0 = [0;0;]

A = [0 1 0 0;
        3*n^2 0 0 2*n;
        0 0 0 1;
        0 -2*n 0 0];

B = [0 0; 
    1 0;
    0 0;
    0 1];

F = eye(4) + A*dt;
G = B*dt;
%u = zeros(2, length(t));

dx = zeros(4,length(t));    % perturbation state vector


% theta = omega * t;
% xnom = r0*cos(n * t);
% ynom = r0*sin(n * t);
% xdotnom = -r0 * n * sin(n * t);
% ydotnom = r0 * n * cos(n * t);

U = zeros(2, length(t));
X(:,1) = X0;
U(:,1) = u0;
r(1) = sqrt(x0^2 + y0^2);
for i = 1:length(t)-1
    dx(:,i) = A * X(:,i) + B * U(:,i);
    X(:,i+1) = X(:,i) + dx(:,i)*dt;
    r(i+1) = sqrt(X(1,i+1)^2 + X(3,i+1)^2);
end

rmin= min(abs(r))   % r_min = 0.4114 km
tmin = find(r==rmin)*dt     % t_min = 0.312 sec

% Plot
figure;
subplot(2, 1, 1);
yyaxis left;
plot(t, X(1,:),t,X(3,:));
title('Total States vs. Time');
ylabel('x (km)');
grid on;
xtitle('Time (s)');
ytitle('States')
yyaxis right;
plot(t, X);
title('Total States vs. Time');
ylabel('x (km)');
legend('$x(t) (km)$','$\dot{x}(t) (km/s)$','$y(t) (km)$','$\dot{y}(t) (km/s)$','Location','best','Interpreter','latex');
grid on;
xtitle('Time (s)');
ytitle('States')

subplot(2, 1, 2);
plot(t, X(2,:)',t,X(4,:));
ylabel('\deltax (km/s)');
grid on;
% subplot(4, 1, 3);
% plot(t, X(3,:)');
% ylabel('y (km)');
% grid on;
% subplot(4, 1, 4);
% plot(t, X(4,:)');
% ylabel('\deltay (km/s)');
% xlabel('Time t (s)');
% grid on;
xlim([0 t(end)])
ylim([-5 30])


figure(2)
plot(X(1,:),X(3,:))
% Calculate perturbation vector
% Initialize perturbation state
% dx(:, 1) = [0, 0.075, 0, -0.021]';
% 
% dx(:,1) = [0; 0; 0; 0];
% % Compute perturbation states over time
% for i = 2:length(t)
%     dx(:, i) = F*dx(:,i-1) + G*u(:, i); 
% end
% 
% 
% 
% Atilde = [0 1 0 0 ; 2*mu/r0^3 0 0 0 ; 0 0 0 1; 0 0 -mu/r0^3 0];
% B = [0 0; 1 0; 0 0; 0 1];
% F = eye(4) + Atilde*dt;
% G = dt * B; 
% 
% % % Preallocate time and perturbation state vectors
% % %t = 0:dt:2*pi/omegaE/4; % Time vector (s)
% % t = 0:dt:14000;
% % dx = zeros(4, length(t)); % Perturbation state vector
% % % Define nominal state and control input over time
% % xnom = [r0*ones(1, length(t)); % x value (km)
% %     zeros(1, length(t)); % x velocity (km/s)
% %     zeros(1, length(t)); % y value (rad)
% %     y0dot*ones(1, length(t))]; % y velocity (rad/s)
% % u = zeros(2, length(t));
% % % Initialize perturbation state
% % dx(:, 1) = [0,0.75,0,-0.021];
% % % Compute perturbation states over time
% % for i = 2:length(ghjny8ttvy )
% %     dx(:, i) = F*dx(:, i - 1) + G*u(:, i);
% % end
% % % Compute total states over time
% 
% 
% x0 = 6678;
% y0 = 0;
% r0 = sqrt(x0^2+y0^2);
% x0dot = 0;
% y0dot = r0 * sqrt(mu/r0^3);
% omega = sqrt(mu/r0^3); % angular velocity of satellite
% % Atilde is time varying so need to create paramaterized variable.... 
% % -- pass in r and t (r0 for nominal)
% % Atilde = [0 1 0 0; 
% %     mu*(2*(cos(omega*t))^2-(sin(omega*t))^2)/r^3 0 3*mu*sin(omega*t)*cos(omega*t)/r^3 0;
% %     0 0 0 1;
% %     3*mu*sin(omega*t)*cos(omega*t)/r^3 0 mu*(2*(sin(omega*t)^2 -(cos(omega*t))^2))/r^3 0];
% A0 = [0 1 0 0 ; 2*mu/r0^3 0 0 0 ; 0 0 0 1; 0 0 -mu/r0^3 0];
% B = [0 0; 1 0; 0 0; 0 1];
% F = eye(4) + Atilde*dt; % pass in r and t to Atilde... 
% G = dt * B; 
% 
% % Define nonlinear model as anonymous function
% pert = [0, 0.075, 0, -0.021]'; % perturbation 
% xdot = @(x, u) [x(2);
%     -mu*x(1)/sqrt((x(1)^2+x(3)^2)^3) + u(1);
%     x(4);
%     -mu*x(3)/sqrt((x(1)^2+x(3)^2)^3) + u(2)];
% % Define initial total state
% X0 = [x0, x0dot, y0, y0dot]' +  pert; 
% % 
% % Define wrapper function for use in ODE45
% xdotwrap = @(t, x) xdot(x, [0; 0]);
% % Compute total state using ODE45
% options = odeset('RelTol', 1e-12); % Set ODE45 options
% [tout, xnl] = ode45(xdotwrap, t, X0, options);
% xnl = xnl'; % Take transpose so that dimensions match with xnom
% 
% % Plot
% figure;
% subplot(4, 1, 1);
% plot(t, xnl(1, :));
% title('Total States vs. Time, Nonlinear Dynanmics - ODE45');
% ylabel('x (km)');
% grid on;
% subplot(4, 1, 2);
% plot(t, xnl(2, :));
% ylabel('\deltax (km/s)');
% grid on;
% subplot(4, 1, 3);
% plot(t, xnl(3, :));
% ylabel('y (km)');
% grid on;
% subplot(4, 1, 4);
% plot(t, xnl(4, :));
% ylabel('\deltay (km/s)');
% xlabel('Time t (s)');
% grid on;
% % end of NL simulation
% 
% % % linearized nominal simulation
% % Define nominal state and control input over time
% dx = zeros(4, length(t)); % Perturbation state vector
% theta = omega * t;
% xnom = r0*cos(omega * t);
% ynom = r0*sin(omega * t);
% xdotnom = -r0 * omega * sin(omega * t);
% ydotnom = r0 * omega * cos(omega * t);
% 
% u = zeros(2, length(t));
% % Initialize perturbation state
% dx(:, 1) = [0, 0.075, 0, -0.021]';
% % Compute perturbation states over time
% for i = 2:length(t)
%     dx(:, i) = F*dx(:, i - 1) + G*u(:, i); % need to pass in r and t to Atilde...
% end
% 
% % Compute total states over time
% x = xnom + dx;
% 
% % Plot
% figure;
% subplot(4, 1, 1);
% plot(t, x(1, :));
% title('Total States vs. Time, Linearized DT SS Model');
% ylabel('x (km)'); ylim([-7000, 7000])
% grid on;
% subplot(4, 1, 2);
% plot(t, x(2, :));
% ylabel('\deltax (km/s)');
% grid on;
% subplot(4, 1, 3);
% plot(t, x(3, :));
% ylabel('y (km)');
% grid on;
% subplot(4, 1, 4);
% plot(t, x(4, :));
% ylabel('\deltay (km/s)');
% xlabel('Time t (s)');
% grid on;
% 
% figure;
% subplot(4, 1, 1);
% plot(t, x(1, :) - xnl(1, 1:length(t)));
% title('State Error vs. Time, Linearized DT SS Model vs. ODE45');
% ylabel('x (km)');
% grid on;
% subplot(4, 1, 2);
% plot(t, x(2, :) - xnl(2, 1:length(t)));
% ylabel('dx (km/s)');
% grid on;
% subplot(4, 1, 3);
% plot(t, x(3, :) - xnl(3, 1:length(t)));
% ylabel('y (km)');
% grid on;
% subplot(4, 1, 4);
% plot(t, x(4, :) - xnl(4, 1:length(t)));
% ylabel('dy (km/s)');
% xlabel('Time t (s)');
% grid on
