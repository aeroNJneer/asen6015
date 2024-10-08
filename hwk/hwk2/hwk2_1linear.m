clear
close all
init;
dt = 0.1;
k = 1;
N = [3,4,5];
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

theta_he = deg2rad(30);      % heading angle in degrees (-rhat * vhat)
Vm = vm * [cos(theta_he) sin(theta_he)];
ydot = -vc*sin(theta_he);

x0 = [0, ydot, 0]';
[y, ydot, nt, nc, temp] = proNav(A,B,x0,R0,vc,tspan,dt,N,k,false);

% for N0 = N
% 
%     for i = 1:length(tspan)
%         nc(i,N-2) = -vm*deg2rad(theta_he)*N0/tgo*(1-i*dt/tgo)^(N0-2)
%     end
% end

% Plot results
figure(1)
plot(tspan,temp(1,:),tspan,temp(2,:),tspan,temp(3,:))
% figure;

% subplot(3, 1, 1);
% plot(tspan, Xi(:,1));
% title(['Proportional Navigation with N0 = ', num2str(N)]);
% ylabel('Vertical Displacement (ft)');
% 
% subplot(3, 1, 2);
% plot(tspan, Xi(:,2));
% ylabel('Vertical Velocity (ft/s)');
% 
% subplot(3, 1, 3);
% plot(tspan, Xi(:,4));
% ylabel('Normal Acceleration Command (ft/s^2)');
% xlabel('Time (s)');



function [y,ydot,nt,nc,temp] = proNav(A,B,x0,R0,vc,tspan,dt,Nvalues,k,augmented)
    theta_he = evalin('base','theta_he');
    for N = Nvalues
        % Initialize state vector
        x = x0;
        r = R0;
        theta = atan2(x(1),r);
        for i = 1:length(tspan)
            % Calculate line-of-sight rate 
            lambda_dot = x(2) / r(i);
       %     lambda_ddot = x(3) / r;     % 2 * lambda_dot *r_dot/r
                                        % nt_theta + (N-2) * r_dot * lambda_dot
            Vc = vc * [cos(theta) sin(theta)];
            % Proportional Navigation command
            if(augmented)
                nc = N * lambda_dot * Vc + 0.5 * N * x(3);
            else
                nc = N * lambda_dot * Vc;
            end

            % State update
            xdot = A * x + B * nc(i);
            x = x + xdot * dt;
            % Store results
            y(N-2,i) = x(1);
            ydot(N-2,i) = x(2);
            ZEM(N-2,i) = -x(1)/(tspan(end) - (i-1)*dt);
            nt(N-2,i) = x(3);
            theta(N-2,i) = atan2(y(i),r(i));
            temp(N-2,i) = nc(i)*tspan(end)/3000/theta(i);
            r(i+1) = r(i-1) - Vc*dt;

        end
        states = [y' ydot' nt' nc' temp'];
    end

    % -nc * tf / vm / theta_he vs t/tf

     
end 