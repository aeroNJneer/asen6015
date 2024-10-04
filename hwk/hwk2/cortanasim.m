close all
clear
g = 32.174;         % gravitational accel
% Simulation parameters
nt = 0 * g;         % acceleration of target
Nvalues = [3 4 5];  % guidance gain
he = deg2rad(30);   % heading error
augmented = false;
dt = 0.01;
dtfine = 0.0001;

% Constants
vm = 3000;      % Missile velocity (ft/s)
vt = 1000;      % Target velocity (ft/s)
R0 = 40000;     % Initial range (ft)

% Initial conditions
xm = 0;        % Missile initial x position (ft)
ym = 0;        % Missile initial y position (ft)
xt = R0;       % Target initial x position (ft)
yt = 0;        % Target initial y position (ft)
beta = 0;       % Target FPA

for N = Nvalues
    % Target velocity components
    Vtx = -vt*cos(beta);
    Vty = vt*sin(beta);
    
    % Relative positions
    Rx = xt - xm;
    Ry = yt - ym;
    R = sqrt(Rx^2 + Ry^2);
    
    % calculate angles
    lambda = atan2(Ry,Rx);
    L = asin(vt*sin(beta+lambda)/vm);  
    theta = lambda + L;

    % mission velocity
    Vmx = vm*cos(theta + he);
    Vmy = vm*sin(theta + he);
    Vx = Vtx - Vmx;     % relative velocity, x_hat
    Vy = Vty - Vmy;     % relative velocity, y_hat
    Vc = -(Rx*Vx + Ry*Vy)/R;  % closing velocity
    lambdad = (Rx*Vy-Ry*Vx)/(R^2);

    i = 1;
    T = 0;
    while Vc>=0
        % Simulation loop for different guidance laws and initial conditions       
        i= i+1;
        if Vc <=1000
            dt = dtfine;
        end
        if augmented
            nc(i,N-2) = N*Vc*lambdad + 0.5*N*nt;
        else
            nc(i,N-2) = N*Vc*lambdad;
        end
        % acceleration components
        amx = -nc*sin(lambda);
        amy = nc*cos(lambda);
        
        betad = nt/vt;
        beta = beta + dt*betad;
        Vtx(i,N-2) = -vt*cos(beta);
        Vty(i,N-2) = vt*sin(beta);  
        Vmx(i,N-2) = Vmx(i-1,N-2) + dt*amx(i,N-2);
        Vmy(i,N-2) = Vmy(i-1,N-2) + dt*amy(i,N-2);
        
        xt(i,N-2) = xt(i-1) + dt*Vtx(i,N-2);
        yt(i,N-2) = yt(i-1) + dt*Vty(i,N-2);
        xm(i,N-2) = xm(i-1,N-2) + dt*Vmx(i,N-2);
        ym(i,N-2) = ym(i-1,N-2) + dt*Vmy(i,N-2);
        Rx = xt(i,N-2) - xm(i,N-2);
        Ry = yt(i,N-2) - ym(i,N-2);
        R = sqrt(Rx^2 + Ry^2);
        Vx = Vtx(i,N-2) - Vmx(i,N-2);
        Vy = Vty(i,N-2) - Vmy(i,N-2);
        lambda = atan2(Ry,Rx);
        lambdad = (Rx*Vy-Ry*Vx)/(R^2);
        T(i,N-2) = T(i-1,N-2) + dt;
        % Calculate closing velocity
        Vc = -(Rx*Vx + Ry*Vy)/R;     
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

