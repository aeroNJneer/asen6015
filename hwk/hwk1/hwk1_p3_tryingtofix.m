clear
close all
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

% hi = HfromRe(Gearth,r1,0)  ; % momentum initial orbit
% hf = HfromRe(Gearth,r2,0)   ;% momentum final orbit
% ht = HfromAe(Gearth,a,e)    ;% momentum transfer orbit
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
mu=Gearth;
w=sqrt(mu/(r1^3));
x0=[r1;0;0;w+(dv1/r1)]; % initial conditions

% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
opts = odeset('RelTol',1e-32,'AbsTol',1e-32);
[t,y] = ode45(@(t,y) myfunc(t,y,mu), tspan, x0, opts);

% final orbit (burn 2) using output from previous ode
tspan = linspace(0,Tf,1000);
xf = [y(1000,1), y(1000,2),y(1000,3),y(1000,4)+(dv2/y(1000,1))];
[tf,yf] = ode45(@(t,y) myfunc(t,y,mu), tspan, xf, opts);

figure(2);
s=stackedplot(t, y);
xlabel('Time t (s)');
s.DisplayLabels={'Radial Position (km)', 'Radial Velocity (km/sec)', 'Angular Position (rad)', 'Angular Velocity (rad/sec)'};
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Ideal Nonlinear State Change After Burn 1 vs. Time (ode45)');

[th,yh] = hohmannTransfer(r1,mu,w,Tof,Tf,dv1,dv2,0);
% Simulate ERROR CONDITIONS

% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
x1e=[r1;dv1*sind(err);0;w+(dv1*cosd(err)/r1)]; % initial conditions
[t1e,y1e] = ode45(@(t,y) myfunc(t,y,mu), tspan, x1e, opts);

% final orbit (burn2) using output from previous ode
tspan = linspace(0,Tf,1000);
x2e=[y1e(1000,1);y1e(1000,2) + (dv2*sind(err));y1e(1000,3);y(1000,4)+(dv2*cosd(err)/y1e(1000,1))]; % initial conditions
[t2e,y2e] = ode45(@(t,y) myfunc(t,y,mu), tspan, x2e, opts);

figure(3);
plot(Rx1',Ry1',y(:,1).*sin(y(:,3)),y(:,1).*cos(y(:,3)),yf(:,1).*sin(yf(:,3)),yf(:,1).*cos(yf(:,3)),y1e(:,1).*sin(y1e(:,3)),y1e(:,1).*cos(y1e(:,3)),y2e(:,1).*sin(y2e(:,3)),y2e(:,1).*cos(y2e(:,3)))
title('Hohmann Transfer with and without thrust err')
legend("initial orbit", "transfer orbit", "final orbit", "transfer orbit w/error", "final orbit w/error")

figure(4);
s=stackedplot(t, y-y1e);
xlabel('Time t (s)');
s.DisplayLabels={'Radial Position (km)', 'Radial Velocity (km/sec)', 'Angular Position (rad)', 'Angular Velocity (rad/sec)'};
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Delta from Nominal States Transfer Orbit');

figure(5);
s=stackedplot(tf, yf-y2e);
xlabel('Time t (s)');
s.DisplayLabels={'Radial Position (km)', 'Radial Velocity (km/sec)', 'Angular Position (rad)', 'Angular Velocity (rad/sec)'};
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Delta from Nominal States - Final Orbit');


%% Sensitivity study
% finite difference dy/dk (ypert-ynom)/(kpert-knom) 

% k = pointing error
% y = position on arrival of final orbit
dk = 1;
dy = [yf(:,1) yf(:,3)] - [y2e(:,1) y2e(:,3)];    % difference in position final states
dy = yf-y2e;
sens = dy/dk;
E = mean(sens)
ddv = deg2rad(0.1);  % -0.3

% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
x1ee=[r1;dv1*sind(ddv);0;w+(dv1*cosd(ddv)/r1)]; % initial conditions
[t,y1ee] = ode45(@(t,y) myfunc(t,y,mu), tspan, x1ee, opts);

% final orbit (burn2)
tspan = linspace(0,Tf,1000);
x2ee=[y1e(1000,1);y1e(1000,2) + (dv2*sind(ddv));y1e(1000,3);y(1000,4)+(dv2*cosd(ddv)/y1e(1000,1))]; % initial conditions
[t,y2ee] = ode45(@(t,y) myfunc(t,y,mu), tspan, x2ee, opts);

YY = mean(y2ee-yf)

Vi = sqrt(var(y2e))
Ve3 = sqrt(var(y2ee))

%% finite burn study




%% functions
function dy=myfunc(t,y,a)
    dy=zeros(4,1);
    dy(1)=y(2);
    dy(2)=y(1).*(y(4)^2)-(a./(y(1)^2));
    dy(3)=y(4);
    dy(4)=(-2.*y(2).*(y(4))./y(1));
end

function [x,y] = circularOrbit(mu,r)
init    
    T = PeriodfromA(mu,r);
    t = linspace(0,T,1000);
    n = MMfromP(T);
    x = r * sin(n.*t);
    y = r * cos(n.*t);
end

function [t,y] = hohmannTransfer(r1,mu,w,tof,to2,dv1,dv2,err)
    if err ~= 0
        Dv1 = [dv1*sind(err) dv1*cosd(err)];
        Dv2 = [dv2*sind(err) dv2*cosd(err)];
    else
        Dv1 = [0 dv1];
        Dv2 = [0 dv2];
    end

    x0=[r1;0;0;w+(dv1/r1)]; % initial conditions
    
    % transfer orbit (burn 1)
    tspan = linspace(0,tof,1000);
    opts = odeset('RelTol',1e-32,'AbsTol',1e-32);
    [t,y] = ode45(@(t,y) myfunc(t,y,mu), tspan, x0, opts);
    
    % final orbit (burn 2) using output from previous ode
    tspan = linspace(0,to2,1000);
    xf = [y(1000,1), y(1000,2),y(1000,3),y(1000,4)+(dv2/y(1000,1))];
    [tf,yf] = ode45(@(t,y) myfunc(t,y,mu), tspan, xf, opts);
    t = [t; t(end)+tf];
    y = [y; yf];
end

function compare_burns(v0,r0,deltaV,duration,acceleration)
    % Constants
    mu = 398600; % Gravitational parameter for Earth (km^3/s^2)
    r0 = 7000; % Initial orbit radius (km)
    v0 = sqrt(mu / r0); % Initial orbit velocity (km/s)
    % burn_duration  Duration of the finite burn (s)
    % acceleration  Constant acceleration (km/s^2)

    % Normalize the burn vector
    burn_vector = deltaV / norm(deltaV);

    % Time settings
    t_impulsive = 0; % Impulsive burn time
    t_finite = linspace(0, burn_duration, 1000); % Finite burn time array

    % Impulsive burn
    v_impulsive = v0 + deltaV * burn_vector(2); % Only tangential component affects velocity
    r_impulsive = r0 + deltaV * burn_vector(1); % Only radial component affects position

    % Finite burn
    v_finite = v0 + acceleration * t_finite * burn_vector(2); % Tangential velocity change
    r_finite = r0 + 0.5 * acceleration * t_finite.^2 * burn_vector(1); % Radial position change

    % Plotting
    figure;
    hold on;
    plot_orbit(r_impulsive, v_impulsive, mu, 'Impulsive Burn');
    plot_orbit(r_finite, v_finite, mu, 'Finite Burn');
    legend('Impulsive Burn', 'Finite Burn');
    title('Comparison of Impulsive and Finite Burns');
    xlabel('X (km)');
    ylabel('Y (km)');
    axis equal;
    grid on;
    hold off;
end

function plot_orbit(r, v, mu, label)
    % Calculate orbital elements
    a = r / (2 - (r * v^2 / mu));
    e = sqrt(1 - (r * v^2 / mu));
    theta = linspace(0, 2 * pi, 1000);
    r_orbit = a * (1 - e^2) ./ (1 + e * cos(theta));

    % Convert to Cartesian coordinates
    x = r_orbit .* cos(theta);
    y = r_orbit .* sin(theta);

    % Plot orbit
    plot(x, y, 'DisplayName', label);
end
