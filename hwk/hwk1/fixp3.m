clear
close all
init;


error = 1; % degree in r_hat direction
% problem 3 - Hohmann transfer
r1 = 200 + Rearth;


k=Gearth;
w=sqrt(k/(r1^3));
x0=[r1;0;0;w+(dv1/r1)]; % initial conditions



% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
opts = odeset('RelTol',1e-32,'AbsTol',1e-32);
[t,y] = ode45(@(t,y) myfunc(t,y,k), tspan, x0, opts);

% final orbit (burn 2)
tspan = linspace(0,Tf,1000);
xf = [y(1000,1), y(1000,2),y(1000,3),y(1000,4)+(dv2/y(1000,1))];
[tf,yf] = ode45(@(t,y) myfunc(t,y,k), tspan, xf, opts);



figure(3);
s=stackedplot(t, y);

figure(4);
plot(y(:,1).*sin(y(:,3)),y(:,1).*cos(y(:,3)),yf(:,1).*sin(yf(:,3)),yf(:,1).*cos(yf(:,3)))
title('Transfer and Final Orbit')



% Simulate ERROR CONDITIONS

% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
x1e=[r1;dv1*sind(error);0;w+(dv1*cosd(error)/r1)]; % initial conditions
[t,y1e] = ode45(@(t,y) myfunc(t,y,k), tspan, x1e, opts);

% final orbit (burn2)
tspan = linspace(0,Tf,1000);
x2e=[y1e(1000,1);y1e(1000,2) + (dv2*sind(error));y1e(1000,3);y(1000,4)+(dv2*cosd(error)/y1e(1000,1))]; % initial conditions
[t,y2e] = ode45(@(t,y) myfunc(t,y,k), tspan, x2e, opts);

figure(6);


@@ -175,6 +187,18 @@ s.DisplayLabels={'Radial Position (km)', 'Radial Velocity (km/sec)', 'Angular Po
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Delta from Nominal States - Final Orbit');














%% Sensitivity study

@@ -192,12 +216,12 @@ error = deg2rad(0.1);  % -0.3
% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
x1ee=[r1;dv1*sind(error);0;w+(dv1*cosd(error)/r1)]; % initial conditions
[t,y1ee] = ode45(@(t,y) myfunc(t,y,k), tspan, x1ee, opts);

% final orbit (burn2)
tspan = linspace(0,Tf,1000);
x2ee=[y1e(1000,1);y1e(1000,2) + (dv2*sind(error));y1e(1000,3);y(1000,4)+(dv2*cosd(error)/y1e(1000,1))]; % initial conditions
[t,y2ee] = ode45(@(t,y) myfunc(t,y,k), tspan, x2ee, opts);

YY = mean(y2ee-yf)


@@ -206,12 +230,14 @@ Ve3 = sqrt(var(y2ee))


%% functions
function dy=myfunc(t,y,a)


    dy=zeros(4,1);
    dy(1)=y(2);
    dy(2)=y(1).*(y(4)^2)-(a./(y(1)^2));
    dy(3)=y(4);
    dy(4)=(-2.*y(2).*(y(4))./y(1));
end

function [x,y] = circularOrbit(mu,r)

@@ -221,4 +247,42 @@ init
    n = MMfromP(T);
    x = r * sin(n.*t);
    y = r * cos(n.*t);
end