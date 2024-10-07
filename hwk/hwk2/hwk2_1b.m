% ASEN-6015 HWK #2 question 1 part b

close all
clear
g = 32.174;         % gravitational accel
% Simulation parameters
nt = 3 * g;         % acceleration of target
N = 4;  % guidance gain
he = deg2rad(30);   % heading error
augmented = true;
dt = 0.1;
dtfine = 0.001;
numsims = 100;   % number of siimulations

sigsq = 0.3;        % random error in tgo (sigma squared)
sigma = sqrt(sigsq);
% Constants
vm = 3000;      % Missile velocity (ft/s)
vt = 1000;      % Target velocity (ft/s)
R0 = 40000;     % Initial range (ft)


for j = 1:numsims
    % Initial conditions
    xm(1,j) = 0;        % Missile initial x position (ft)
    ym(1,j) = 0;        % Missile initial y position (ft)
    xt(1,j) = R0;       % Target initial x position (ft)
    yt(1,j) = 0;        % Target initial y position (ft)
    beta = 0;       % Target FPA
    i = 1;

    % Target velocity components
    Vtx(i,j) = -vt*cos(beta);
    Vty(i,j) = vt*sin(beta);
    
    % Relative positions
    Rx(i,j) = xt(i,j) - xm(i,j);
    Ry(i,j) = yt(i,j) - ym(i,j);
    R(i,j) = sqrt(Rx(i,j)^2 + Ry(i,j)^2);
    
    % calculate angles
    lambda = atan2(Ry(i,j),Rx(i,j));
    L = asin(vt*sin(beta+lambda)/vm);  
    theta = lambda + L;

    % mission velocity
    Vmx(i,j) = vm*cos(theta + he);
    Vmy(i,j) = vm*sin(theta + he);
    Vx = Vtx(i,j) - Vmx(i,j);     % relative velocity, x_hat
    Vy = Vty(i,j) - Vmy(i,j);     % relative velocity, y_hat
    Vc = -(Rx(i,j)*Vx + Ry(i,j)*Vy)/R(i,j);  % closing velocity
    lambdad = (Rx(i,j)*Vy-Ry(i,j)*Vx)/(R(i,j)^2);
    tgoerr = sigma*randn;
    tgo = R(i,j)/Vc + tgoerr;
    ZEM = Vc * tgo;

    T(i,j) = 0;
    while Vc > 0
        % Simulation loop for different guidance laws and initial conditions       
        i= i+1;
        if Vc <=100
            dt = dtfine;
        end
        if augmented
            nc(i) = N*ZEM/tgo^2 + 0.5*N*nt;
        else
            nc(i) = N*Vc*lambdad;
            %nc(i) = N*ZEM/tgo^2;
        end
        % acceleration components
        amx = -nc*sin(lambda);
        amy = nc*cos(lambda);
        
        betad = nt/vt;
        beta = beta + dt*betad;
        Vtx(i) = -vt*cos(beta);
        Vty(i) = vt*sin(beta);  
        Vmx(i) = Vmx(i-1) + dt*amx(i);
        Vmy(i) = Vmy(i-1) + dt*amy(i);
        
        xt(i) = xt(i-1) + dt*Vtx(i);
        yt(i) = yt(i-1) + dt*Vty(i);
        xm(i) = xm(i-1) + dt*Vmx(i);
        ym(i) = ym(i-1) + dt*Vmy(i);
        Rx(i) = xt(i) - xm(i);
        Ry(i) = yt(i) - ym(i);
        R(i) = sqrt(Rx(i)^2 + Ry(i)^2);
        Vx = Vtx(i) - Vmx(i);
        Vy = Vty(i) - Vmy(i);
        lambda = atan2(Ry(i),Rx(i));
        lambdad = (Rx(i)*Vy-Ry(i)*Vx)/(R(i)^2);
        T(i) = T(i-1) + dt;
        % Calculate closing velocity
        Vc = -(Rx(i)*Vx + Ry(i)*Vy)/R(i);     
        tgo = R(i)/Vc +tgoerr;
        ZEM = tgo * (Rx(i)*Vy - Ry(i)*Vx)/R(i);
    end
    ZEMfinal(j) = ZEM;
    Vcfinal(j) = Vc;
    Tcollision(j) = T(i);
end
   dev = std(ZEMfinal)
plot(1:length(Tcollision),Tcollision)


% % Plot results for PN
% figure(1);
% subplot(3, 1, 1);
% plot(T,nc/g);
% xlabel('Time (s)');
% ylabel('Acceleration - n_c (g''s)');
% title('ProNav Control Acceleration n_c (ft/s^2)');
% grid on;
% 
% 
% subplot(3,1,2);
% plot(T,ym);
% xlabel('Time (s)');
% ylabel('Vertical Position - y (m)');
% title('ProNav Vertical Displacement (ft)');
% grid on;
% 
% subplot(3,1,3);
% plot(T,Vmy);
% title('ProNav Missile Vertical Velocity - \Deltay (ft/s)');
% xlabel('Time (s)');
% ylabel('Vertical Velocity (ft/s)')
% 
% figure(2)
% plot(xt,yt,xm,ym);
% legend('Target','Missile');
% title('ProNav Collision Path N=4');
% xlabel('Horizontal Position (ft)');
% ylabel('Vertical Position (ft)');
% 
