close all
clear
g = 32.174;         % gravitational accel
% Simulation parameters
nt = 3 * g;         % acceleration of target
Nvalues = [3 4 5];  % guidance gain
he = deg2rad(30);   % heading error
augmented = true;
dt = 0.01;
dtfine = 0.0001;

% Constants
vm = 3000;      % Missile velocity (ft/s)
vt = 1000;      % Target velocity (ft/s)
R0 = 40000;     % Initial range (ft)

for N = Nvalues
    % Initial conditions
    xm(1,N-2) = 0;        % Missile initial x position (ft)
    ym(1,N-2) = 0;        % Missile initial y position (ft)
    xt(1,N-2) = R0;       % Target initial x position (ft)
    yt(1,N-2) = 0;        % Target initial y position (ft)
    beta = 0;       % Target FPA
    i = 1;


    % Target velocity components
    Vtx(i,N-2) = -vt*cos(beta);
    Vty(i,N-2) = vt*sin(beta);
    
    % Relative positions
    Rx(i,N-2) = xt(i,N-2) - xm(i,N-2);
    Ry(i,N-2) = yt(i,N-2) - ym(i,N-2);
    R(i,N-2) = sqrt(Rx(i,N-2)^2 + Ry(i,N-2)^2);
    
    % calculate angles
    lambda = atan2(Ry(i,N-2),Rx(i,N-2));
    L = asin(vt*sin(beta+lambda)/vm);  
    theta = lambda + L;

    % mission velocity
    Vmx(i,N-2) = vm*cos(theta + he);
    Vmy(i,N-2) = vm*sin(theta + he);
    Vx = Vtx(i,N-2) - Vmx(i,N-2);     % relative velocity, x_hat
    Vy = Vty(i,N-2) - Vmy(i,N-2);     % relative velocity, y_hat
    Vc = -(Rx(i,N-2)*Vx + Ry(i,N-2)*Vy)/R(i,N-2);  % closing velocity
    lambdad = (Rx(i,N-2)*Vy-Ry(i,N-2)*Vx)/(R(i,N-2)^2);
    tgo = R(i,N-2)/Vc;
    ZEM = Vc * tgo;

    T(i,N-2) = 0;
    while Vc > 0
        % Simulation loop for different guidance laws and initial conditions       
        i= i+1;
        if Vc <=1000
            dt = dtfine;
        end
        if augmented
            nc(i,N-2) = N*Vc*lambdad + 0.5*N*nt;
        else
            nc(i,N-2) = N*Vc*lambdad;
            %nc(i,N-2) = N*ZEM/tgo^2;
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
        Rx(i,N-2) = xt(i,N-2) - xm(i,N-2);
        Ry(i,N-2) = yt(i,N-2) - ym(i,N-2);
        R(i,N-2) = sqrt(Rx(i,N-2)^2 + Ry(i,N-2)^2);
        Vx = Vtx(i,N-2) - Vmx(i,N-2);
        Vy = Vty(i,N-2) - Vmy(i,N-2);
        lambda = atan2(Ry(i,N-2),Rx(i,N-2));
        lambdad = (Rx(i,N-2)*Vy-Ry(i,N-2)*Vx)/(R(i,N-2)^2);
        T(i,N-2) = T(i-1,N-2) + dt;
        % Calculate closing velocity
        Vc = -(Rx(i,N-2)*Vx + Ry(i,N-2)*Vy)/R(i,N-2);     
%         tgo = R(i,N-2)/Vc;
%         ZEM = Vc * tgo;
        end 
end

output=[T,xm,ym,xt,yt,nc];
save 'c:\src\asen6015\datfil.txt' output -ascii -tabs

% Plot results for PN
figure(1);
subplot(3, 1, 1);
plot(T(:,1), nc(:,1)/g,T(:,1),nc(:,2)/g,T(:,1),nc(:,3)/g);
xlabel('Time (s)');
ylabel('Acceleration - n_c (g''s)');
title('ProNav Control Acceleration n_c (ft/s^2)');
legend('N=3', 'N=4', 'N=5','Location','southeast');
grid on;


subplot(3,1,2);
plot(T(:,1), ym(:,1),T(:,1),ym(:,2),T(:,1),ym(:,3));
xlabel('Time (s)');
ylabel('Vertical Position - y (m)');
title('ProNav Vertical Displacement (ft)');
legend('N=3', 'N=4', 'N=5','Location','northeast');
grid on;

subplot(3,1,3);
plot(T(:,1), Vmy(:,1),T(:,1),Vmy(:,2),T(:,1),Vmy(:,3));
title('ProNav Missile Vertical Velocity - \Deltay (ft/s)');
xlabel('Time (s)');
ylabel('Vertical Velocity (ft/s)')
legend('N=3', 'N=4', 'N=5');

% remove entries with zerosym3 = noZeros(ym(:,1));
xt3 = noZeros(xt(:,1));
xt4 = noZeros(xt(:,2));
xt5 = noZeros(xt(:,3));
xm3 = [0;noZeros(xm(:,1))];
xm4 = [0;noZeros(xm(:,2))];
xm5 = [0;noZeros(xm(:,3))];
yt3 = [0;noZeros(yt(:,1))];
yt4 = [0;noZeros(yt(:,2))];
yt5 = [0;noZeros(yt(:,3))];
ym3 = [0;noZeros(ym(:,1))];
ym4 = [0;noZeros(ym(:,2))];
ym5 = [0;noZeros(ym(:,3))];

figure(2)
subplot(3,1,1)
plot(xt3,yt3,xm3,ym3);
legend('Target','Missile');
title('ProNav Collision Path N=3');
xlabel('Horizontal Position (ft)');
ylabel('Vertical Position (ft)');
% ylim([-1000 4000]);

subplot(3,1,2)
plot(xt4,yt4,xm4,ym4);
legend('Target','Missile');
title('ProNav Collision Path N=4');
xlabel('Horizontal Position (ft)');
ylabel('Vertical Position (ft)')
% ylim([-1000 4000]);

subplot(3,1,3)
plot(xt5,yt5,xm5,ym5);
legend('Target','Missile');
title('ProNav Collision Path N=5');
xlabel('Horizontal Position (ft)');
ylabel('Vertical Position (ft)')
% ylim([-1000 4000]);

function x = noZeros(m)
    x = m(m~=0);
end