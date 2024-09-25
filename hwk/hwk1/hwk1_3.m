clear
close all
clc
init;
err = 1; % degree in r_hat direction
% problem 3 - Hohmann transfer

r1 = 200 + Rearth;
r2 = 35786 + Rearth;
% circular orbit velocity
v1 = sqrt(Gearth/r1);
v2 = sqrt(Gearth/r2);

% orbital elements of transfer
e = EccfromRpRa(r1,r2);
a = (r1+r2)/2;
Tof = PeriodfromA(Gearth, a)/2;  % time of flight
Ti = PeriodfromA(Gearth, r1);   % period of lower orbit
Tf = PeriodfromA(Gearth, r2);   % period of higher orbit
rp = r1;
ra = r2;

vt1 = VfromRA(Gearth,r1,a);
vt2 = VfromRA(Gearth,r2,a);

hi = HfromRe(Gearth,r1,0)  ; % momentum initial orbit
hf = HfromRe(Gearth,r2,0)   ;% momentum final orbit
ht = HfromAe(Gearth,a,e)    ;% momentum transfer orbit
% delta v impulse burns
dv1 = vt1-v1;
dv2 = v2-vt2;

[Rx1,Ry1] = circularOrbit(Gearth,r1);

%% Plot the Hohmann Transfer Orbits

% find eccentric anomaly of transfer
t = linspace(0,Tof,1000);

EAt = []; TA = []; rt = [];
for i = 1:length(t)
    EAt(i) = FindEccAnomaly(t(i),a,e,Gearth);
    TAt(i) = deg2rad(TAfromEA(e,EAt(i)));
    rt(i) = RfromEA(a,e,EAt(i));
end

% data points for initial circular orbit
t = linspace(0,Ti,1000);
Tai = 2 * pi * t / Ti;

% data points for final circular orbit
t = linspace(0,Tf,1000);
Taf = 2 * pi * t / Tf;

% Plot the orbit
figure(1);
polarplot(Tai, r1 * ones(size(Tai)),Taf,r2*ones(size(Taf)),TAt,rt);
title('Ideal Initial, Transfer and Final Orbits');

%% Simulate transfer orbit
% No errors in thrust vector
k=Gearth;
w=sqrt(k/(r1^3));

[t,y] = hohmannTransfer(r1,k,w,Tof,Tf,dv1,dv2,0);

figure(2);
s=stackedplot(t(1:1000), y(1:1000,:));
xlabel('Time t (s)');
s.DisplayLabels={'Radial Position (km)', 'Radial Velocity (km/sec)', 'Angular Position (rad)', 'Angular Velocity (rad/sec)'};
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Ideal Nonlinear State Change After Burn 1 vs. Time (ode45)');

% Simulate ERROR CONDITIONS
[te,ye] = hohmannTransfer(r1,k,w,Tof,Tf,dv1,dv2,1);

figure(3);
plot(Rx1',Ry1',y(1:1000,1).*sin(y(1:1000,3)),y(1:1000,1).*cos(y(1:1000,3)),y(1001:end,1).*sin(y(1001:end,3)),y(1001:end,1).*cos(y(1001:end,3)),ye(1:1000,1).*sin(ye(1:1000,3)),ye(1:1000,1).*cos(ye(1:1000,3)),ye(1001:end,1).*sin(ye(1001:end,3)),ye(1001:end,1).*cos(ye(1001:end,3)))
title('Hohmann Transfer with and without thrust error')
legend("initial orbit", "transfer orbit", "final orbit", "transfer orbit w/error", "final orbit w/error")

dy = y-ye;
figure(4);
s=stackedplot(t, y-ye);
xlabel('Time t (s)');
s.DisplayLabels={'Radial Position (km)', 'Radial Velocity (km/sec)', 'Angular Position (rad)', 'Angular Velocity (rad/sec)'};
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Delta from Nominal States');

X = y(1,end) * [cos(y(3,end)) sin(y(3,end)) 0];
X_dot = y(2,end)* [cos(y(3,end)) sin(y(3,end)) 0];

%% Sensitivity study
% finite difference dy/dk (ypert-ynom)/(kpert-knom) 
% k = pointing error
% y = position on arrival of final orbit

dk = 1;
dy = [y(1001:end,1) y(1001:end,3)] - [ye(1001:end,1) ye(1001:end,3)];    % difference in position final states
dy = y(1001:end,1) - ye(1001:end,1);
sens_rf = mean(dy/dk)

% transfer orbit (0.1 degree error in radial direction)
dk1 = 0.1;
[t1,y1] = hohmannTransfer(r1,k,w,Tof,Tf,dv1,dv2,dk1);
dy1 = y(1001:end,1) - y1(1001:end,1);
exp_sens1 = sens_rf*dk1
sens_rf1 = mean(dy1/dk1)*dk1

% transfer orbit (-0.3 degree error in radial direction)
dk3 = -0.3;
[t3,y3] = hohmannTransfer(r1,k,w,Tof,Tf,dv1,dv2,dk3);
dy3 = y(1001:end,1) - y3(1001:end,1);
exp_sens2 = sens_rf*(-0.3)
sens_rf2 = mean(dy3/dk3)*dk3


%% finite burn

% [r,r_dot, theta, theta_dot, h, h_dot]
<<<<<<< Updated upstream
accel = 0.03;             % constant acceleration 30 m/s
dv = [0, 0, 0, dv1, 0, 0]';         % burn vector
btime = norm(dv)/norm(accel);
r0 = [r1, 0,0,0,0, 0]';
v0 = [0, v1, 0]';
=======
accel = 0.03;             % constant acceleration 30 m/s^
dv = [0, dv1, 0]';        % burn vector
btime = norm(dv)/norm(accel);
r0 = [r1, 0, 0]';
v0 = [0, v1, 0]';  % at periapsis
>>>>>>> Stashed changes

% solve ODEs
tspan = linspace(0,btime,1000);
x1 = [r0;v0/r1];
[tfin,yfin] = ode45(@(tfin,yfin) finiteburn(tfin,yfin,dv,accel), tspan, [r0 ;v0]);

<<<<<<< Updated upstream
r = solext.y()
=======
% coasting after finite burn using previous ode
tspan = linspace(0,Tof,1000);
r2 =[yfin(1000,1); yfin(1000,2); yfin(1000,3)]; % initial conditions
v2 = [yfin(1000,4); yfin(1000,5); yfin(1000,6)]; 
[t2fin,y2fin] = ode45(@(t2fin,y2fin) finiteburn(t2fin,y2fin,0,0), tspan, [r2 ;v2]);
>>>>>>> Stashed changes

yfinite = [yfin;y2fin];
X = [];
Q = [];
% Extract position components
for i = 1:length(yfinite)
    [raan inc theta] = RSW2Euler(yfinite(i,1:3),yfinite(i,4:6));
    C = Euler3132C([raan inc theta]);
    X = [X C*yfinite(i,1:3)'];
end

% Plot the 3D orbit
figure(5);
plot3(X(1,:),X(2,:),X(3,:));
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Orbit');

%% functions
function dy=myfunc(t,y,a)
    dy=zeros(4,1);
    dy(1)=y(2);
    dy(2)=y(1).*(y(4)^2)-(a./(y(1)^2));
    dy(3)=y(4);
    dy(4)=(-2.*y(2).*(y(4))./y(1));
end

function dy=dywithdv(t,y,ft,ut,dv1,dv2)
a = evalin("base",'k');
urdot = 0;
uvdot = 0;

   % urdot = interp1(ft,ut(:,1),t);
    %uvdot = interp1(ft,ut(:,2),t);
    if t==0
        urdot = 0;
        uvdot = dv1;
        y(4) = y(4) + uvdot/y(1);
    end

    dy=zeros(4,1);
    dy(1)=y(2);
    dy(2)=y(1).*(y(4)^2)-(a./(y(1)^2)) + urdot;
    dy(3)=y(4);
    dy(4)=(-2.*y(2).*(y(4))./y(1)) + uvdot;
end

function [x,y] = circularOrbit(mu,r)
init    
    T = PeriodfromA(mu,r);
    t = linspace(0,T,1000);
    n = MMfromP(T);
    x = r * sin(n.*t);
    y = r * cos(n.*t);
end

function [t,y] = hohmannTransfer(r1,a,w,tof,to2,dv1,dv2,err)
    if err ~= 0
        Dv1 = [dv1*sind(err) dv1*cosd(err)];
        Dv2 = [dv2*sind(err) dv2*cosd(err)];
    else
        Dv1 = [0 dv1];
        Dv2 = [0 dv2];
    end

    % transfer orbit (burn 1)
    tspan = linspace(0,tof,1000);
    x1 =[r1;Dv1(1);0;w+Dv1(2)/r1];  % initial conditions with error in radial direction
    [t1,y1] = ode45(@(t1,y1) myfunc(tspan,y1,a), tspan, x1);
   
    % final orbit after 2nd burn using output from previous ode
    tspan = linspace(0,to2,1000);
    x2=[y1(1000,1); y1(1000,2)+Dv2(1); y1(1000,3); y1(1000,4)+Dv2(2)/y1(1000,1)]; % initial conditions
    [t2,y2] = ode45(@(t2,y2) myfunc(tspan,y2,a), tspan, x2);

    t = [t1; t1(end)+t2];
    y = [y1; y2];
end

function dydt = finiteburn(t,y,deltaV,accel)
    % Constants
    init;
    mu = Gearth; % Gravitational parameter for Earth (km^3/s^2)
    r = y(1:3);
    v = y(4:6);
    DV = deltaV/norm(deltaV);   % direction of deltav
    Accel = accel*v/norm(v);           % acceleration along deltaV vector
    r_hat = r/norm(r);          % radial direction
    % EOM
    r_dot = v;
    v_dot = (-mu/(norm(r)^3)).*r_hat  + Accel;
    dydt = [r_dot; v_dot];
end

%%