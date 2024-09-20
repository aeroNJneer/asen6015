
clear
close all
init;
error = 1; % degree in r_hat direction
% problem 3 - Hohmann transfer


r1 = 200 + Rearth;
r2 = 35786;
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
x0=[r1;0;0;w+(dv1/r1)]; % initial conditions

% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
opts = odeset('RelTol',1e-15,'AbsTol',1e-15);
[t,y] = ode45(@(t,y) myfunc(t,y,k), tspan, x0, opts);

% final orbit (burn 2) using output from previous ode
tspan = linspace(0,Tf,1000);
xf = [y(1000,1), y(1000,2),y(1000,3),y(1000,4)+(dv2/y(1000,1))];
[tf,yf] = ode45(@(t,y) myfunc(t,y,k), tspan, xf, opts);

figure(3);
s=stackedplot(t, y);
xlabel('Time t (s)');
s.DisplayLabels={'Radial Position (km)', 'Radial Velocity (km/sec)', 'Angular Position (rad)', 'Angular Velocity (rad/sec)'};
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Ideal Nonlinear State Change After Burn 1 vs. Time (ode45)');

% Simulate ERROR CONDITIONS

% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
x1e=[r1;dv1*sind(error);0;w+(dv1*cosd(error)/r1)]; % initial conditions
[t,y1e] = ode45(@(t,y) myfunc(t,y,k), tspan, x1e, opts);

% final orbit (burn2) using output from previous ode
tspan = linspace(0,Tf,1000);
x2e=[y1e(1000,1);y1e(1000,2) + (dv2*sind(error));y1e(1000,3);y(1000,4)+(dv2*cosd(error)/y1e(1000,1))]; % initial conditions
[t,y2e] = ode45(@(t,y) myfunc(t,y,k), tspan, x2e, opts);

figure(4);

plot(Rx1',Ry1',y(:,1).*sin(y(:,3)),y(:,1).*cos(y(:,3)),yf(:,1).*sin(yf(:,3)),yf(:,1).*cos(yf(:,3)),y1e(:,1).*sin(y1e(:,3)),y1e(:,1).*cos(y1e(:,3)),y2e(:,1).*sin(y2e(:,3)),y2e(:,1).*cos(y2e(:,3)))
title('Hohmann Transfer with and without thrust error')
legend("initial orbit", "transfer orbit", "final orbit", "transfer orbit w/error", "final orbit w/error")

figure(5);
s=stackedplot(t, y-y1e);
xlabel('Time t (s)');
s.DisplayLabels={'Radial Position (km)', 'Radial Velocity (km/sec)', 'Angular Position (rad)', 'Angular Velocity (rad/sec)'};
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Delta from Nominal States Transfer Orbit');


figure(8);
s=stackedplot(t, yf-y2e);
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
[t,y1ee] = ode45(@(t,y) myfunc(t,y,k), tspan, x1ee, opts);

% final orbit (burn2)
tspan = linspace(0,Tf,1000);
x2ee=[y1e(1000,1);y1e(1000,2) + (dv2*sind(ddv));y1e(1000,3);y(1000,4)+(dv2*cosd(ddv)/y1e(1000,1))]; % initial conditions
[t,y2ee] = ode45(@(t,y) myfunc(t,y,k), tspan, x2ee, opts);

YY = mean(y2ee-yf);

Vi = sqrt(var(y2e));
Ve3 = sqrt(var(y2ee));

%% input control added to ode

totalt = int32(Tf);
dt = 10;
timespan1 =linspace(0,Tof,1000);
timespan2 = linspace(Tof,Tf,1000);
timespan = [timespan1 timespan2]; 
ut = zeros(length(timespan),2);
ut(1,2) = dv1;
ut(1000,2) = dv2;

xx=[r1;0;0;w]; % initial conditions
%for i = 1:length(t)
    [t,y1ee] = ode45(@(t,y) dywithdv(t,y,timespan,ut,dv1,dv2), t, xx, opts);
 %   xx = y1ee(end,:);
%end


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