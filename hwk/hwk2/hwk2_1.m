clear
close all
init;
dt = 0.1;
k = 1;
N = [3,4,5];
y0 = 0;
% Definition of state space system [y ydot nt]
% yddot = nt-nc     ntdot = 0
A = [0 1 0;
    0 0 1;
    0 0 0];
B = [0 -1 0]';

vm = 3000;      % velocity of missile ft/s
vt = 1000;      % velocity of target ft/s
vc = vm + vt;   % closing speed ft/s
R0 = 40000;     % ft
y0 = 0 ;        % no vertical component
tgo = R0/vc;    % time to go
tspan = 0:dt:tgo;   

% part i: theta_he = 30 degrees, nt = 0;
% cos(theta_he) = dot(r,v)/R0

theta_he = deg2rad(30);      % heading angle in degrees
ydot = vm*sin(theta_he);
x0 = [0, ydot, 0]';
Xi = proNav(A,B,x0,R0,vc,tspan,dt,N,k,false);

% for N0 = N
% 
%     for i = 1:length(tspan)
%         nc(i,N-2) = -vm*deg2rad(theta_he)*N0/tgo*(1-i*dt/tgo)^(N0-2)
%     end
% end

% Plot results
figure;
subplot(3, 1, 1);
plot(tspan, Xi(:,1));
title(['Proportional Navigation with N0 = ', num2str(N)]);
ylabel('Vertical Displacement (ft)');

subplot(3, 1, 2);
plot(tspan, Xi(:,2));
ylabel('Vertical Velocity (ft/s)');

subplot(3, 1, 3);
plot(tspan, Xi(:,4));
ylabel('Normal Acceleration Command (ft/s^2)');
xlabel('Time (s)');

function states = proNav(A,B,x0,R0,vc,tspan,dt,Nvalues,k,augmented)
    
    for N = Nvalues
        % Initialize state vector
        x = x0;
        r = R0;
        
        for i = 1:length(tspan)
            % Calculate line-of-sight rate 
            lambda_dot = x(2) / r;
            lambda_ddot = x(3) / r;

            % Proportional Navigation command
            if(augmented)
                nc(i) = N * lambda_dot + k * lambda_ddot;
            else
                nc(i) = N * lambda_dot;
            end

            % State update
            r = r - vc*dt;
            xdot = A * x + B * nc(i);
            x = x + xdot * dt;
            % Store results
            y(i) = x(1);
            ydot(i) = x(2);
            nt(i) = x(3);

        end
        states = [y' ydot' nt' nc'];
    end
end 