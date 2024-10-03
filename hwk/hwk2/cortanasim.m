close all
clear

% Constants
vm = 3000;      % Missile velocity (ft/s)
vt = 1000;      % Target velocity (ft/s)
vc = vm + vt;   % Closing speed (ft/s)
R0 = 40000;     % Initial range (ft)
tgo = R0 / vc;  % Time to go (s)
dt = 0.1;      % Time step (s)
tfinal = tgo;  % Simulation time (s)

% Initial conditions
xm0 = 0;        % Missile initial x position (ft)
ym0 = 0;        % Missile initial y position (ft)
xt0 = R0;       % Target initial x position (ft)
yt0 = 0;        % Target initial y position (ft)
thetam = 0;     %deg2rad(30); % Missile initial heading (30 degrees in radians)
thetat = pi;     % Target initial heading (rad)
thetahe = deg2rad(30);  % Missile heading error
ydot0 = vm * sin(deg2rad(30)); % Initial y velocity (ft/s) for HE = 30 degrees
nt0 = 0;        % Initial target acceleration ft/s^2)

% Target velocity components
Vtx0 = vt*cos(thetat);
Vty0 = vt*sin(thetat);

% Relative positions
Rx = xt0 - xm0;
Ry = yt0 - ym0;

% LOS angle
lambda = atan2(Ry,Rx);

% lead angle
L = 0;  % asin(vt*sin(thetat + lambda)/vm);   % for linearized system this is 0

% mission velocity
Vmx0 = vm*cos(lambda + L + thetahe);
Vmy0 = vm*sin(lambda + L + thetahe);

X0 = [xm0 ym0 xt0 yt0 Vmx0 Vmy0 Vtx0 Vty0 thetam thetat thetahe nt0];

% Time vector
t = 0:dt:tfinal;
A = [0 1 0;
    0 0 1;
    0 0 0];
B = [0 -1 0]';

% Navigation constants
Nvalues = [3, 4, 5];

[Xm Ym Xt Xt Vmx Vmy Vtx Vty nc ] = proNavSimulation(X0,Nvalues,t,dt,false);

% Simulation loop for different guidance laws and initial conditions
function [xm,ym,xt,yt,nc,P] = proNavSimulation(X0,Nvalues,t,dt,augmented)
    vm = evalin('base','vm');
    vt = evalin('base','vt');

    % Preallocate arrays for results
    n = length(Nvalues);
    xm = zeros(length(t),n);
    ym = zeros(length(t),n);
    xt = zeros(length(t),n);
    yt = zeros(length(t),n);
    ydot = zeros(length(t),n);
    nc = zeros(length(t),n);
    
    for N = Nvalues
        % Initial conditions
        xm(1,N-2) = X0(1);
        ym(1,N-2) = X0(2);
        xt(1,N-2) = X0(3);
        yt(1,N-2) = X0(4);
        vmx(1,N-2) = X0(5);
        vmy(1,N-2) = X0(6);
        vtx(1,N-2) = X0(7);
        vty(1,N-2) = X0(8);
        thetam = X0(9);
        thetat = X0(10);
        thetahe = X0(11);
        nT = X0(12);
    
        for i = 2:length(t)
            Vt = norm([vtx(1,N-2) vty(1,N-2)]);
            Rx = xt(i-1,N-2) - xm(i-1,N-2);
            Ry = yt(i-1,N-2) - ym(i-1,N-2);
            R = norm([Rx Ry]);
            Vx = vtx(i-1,N-2) - vmx(i-1,N-2);
            Vy = vty(i-1,N-2) - vmy(i-1,N-2);
            
            % calculate lambda, lambda dot
            lambda = atan2(Ry,Rx);
            lambda_dot = (Rx*Vy + Ry*Vx)/R;

            % Calculate closing velocity
            Vc = (Rx*Vx + Ry*Vy)/R;

            % Calculate line-of-sight (LOS) angle
            LOSangle = atan2(yt(i-1) - ym(i-1), xt(i-1) - xm(i-1));
            
            % Calculate heading error
            thetahe = LOSangle - thetam;
            if augmented
                % Augmented ProNav
                nc(i,N-2) = N * thetahe * Vc + 0.5 * N * nT;
            else
                % Proportional Navigation
                nc(i,N-2) = N * Vc * thetahe;
            end
            
            % Update target position
            xt(i,N-2) = xt(i-1,N-2) + vt * cos(thetat) * dt;
            yt(i,N-2) = yt(i-1,N-2) + vt * sin(thetat) * dt;
               
            % Update missile heading
            thetam = thetam + nc(i,N-2) / vm * dt;
            
            % Update missile position
            xm(i,N-2) = xm(i-1,N-2) + vm * cos(thetam) * dt;
            ym(i,N-2) = ym(i-1,N-2) + vm * sin(thetam) * dt;

            % pos = [xt(i,N-2)-xm(i,N-2), yt(i,N-2)-ym(i,N-2)];
            % P(i,N-2) = (nc(i-1,N-2)*t(i)/(-vm*thetam*N));
            % if norm(pos)< 0.1
            %     break;
            % end
        end
    end
    % Plot results for PN
    figure(1);
    subplot(2, 1, 1);
    plot(t, ym);
    xlabel('Time (s)');
    ylabel('y (m)');
    title('Proportional Navigation - Vertical Displacement (ft)');
    legend('N=3', 'N=4', 'N=5','Location','west');
    grid on;

    subplot(2,1,2);
    plot(t,xm,t,xt)
    xlabel('Time (s)');
    ylabel('ProNav Control Acceleration nc (ft/s^2)');
    legend('N=3', 'N=4', 'N=5');
    subplot(2, 1, 2);
    plot(t, nc, 'DisplayName', ['N = ', num2str(N)]);
    xlabel('Time (s)');
    ylabel('ProNav Control Acceleration nc (ft/s^2)');
    legend('N=3', 'N=4', 'N=5');
    grid on;

    figure(2)
    subplot(2,1,1);
    plot(t,xm,t,xt)
    xlabel('Time (s)');
    ylabel('ProNav Control Acceleration - Horizontal Displacement (ft)');
    legend('N=3', 'N=4', 'N=5');
    
    subplot(2,1,2)
    plot(t/t(end),P)
    xlabel('Time (s)');
    ylabel('ProNav Control Acceleration nc (ft/s^2)');
    legend('N=3', 'N=4', 'N=5');

    figure(4)
    plot(xm,ym,xt,yt)

end 

