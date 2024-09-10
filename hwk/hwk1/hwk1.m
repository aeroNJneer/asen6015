clear
close all
init;

u = [0;0];
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

% simulate hohmann transfer trajectory

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
% 
% i = linspace(0,Tof,1000);
% EAs = []; TAs = []; rts = []; xs = []; ys = [];
% 
% % initial conditions
% dv1err = 0;
% 
% h0 = HfromP(Gearth,r1); % h of first orbit
% TA0 = 0;
% rx0 = r1;  
% ry0 = 0;
% vx0 = (v1 + dv1) * sin(TA0 + dv1err);
% vy0 = (v1 + dv1) * cos(TA0 + dv1err);
% 
% R = [r1 0 0];
% V = [vx0 vy0 0];
% 
% H = HfromRV(R,V)
% V0 = sqrt(vx0^2 + vy0^2);
% a = AfromVR(Gearth,V0, r1);
% e = EccfromArp(a,r1);
% t = linspace(0,Tof,1000);
% for i = 1:length(t)
%     EAt(i) = FindEccAnomaly(t(i),a,e,Gearth);
%     TAt(i) = deg2rad(TAfromEA(e,EAt(i)));
%     rt(i) = RfromEA(a,e,EAt(i));
%     vs(i) = VfromRA(Gearth,rt(i),a);
% end
% vxf = (vs(i) + dv2) * sin(TAt(i));
% vyf = (vs(i) + dv2) * cos(TAt(i));
% Rx = rt.*cos(TAt);
% Ry = rt.*sin(TAt);
% 
  [Rx1,Ry1] = circularOrbit(Gearth,r1);
% [Rx2,Ry2] = circularOrbit(Gearth,rt(i));
% 
% 
% % initial conditions
% dv1err = deg2rad(1);
% dv1err = 0;
% vx0 = (v1 + dv1) * sin(TA0 + dv1err);
% vy0 = (v1 + dv1) * cos(TA0 + dv1err);
% V0 = sqrt(vx0^2 + vy0^2);
% a = AfromVR(Gearth,V0, r1);
% e = EccfromArp(a,r1);
% t = linspace(0,Tof,1000);
% for i = 1:length(t)
%     EAe(i) = FindEccAnomaly(t(i),a,e,Gearth);
%     TAe(i) = deg2rad(TAfromEA(e,EAe(i)));
%     rt(i) = RfromEA(a,e,EAe(i));
%     vs(i) = VfromRA(Gearth,rt(i),a);
% end
% vxe = (vs(i) + dv2) * sin(TAe(i)+dv1err);
% vye = (vs(i) + dv2) * cos(TAe(i)+dv1err);
% Rexf = rt.*cos(TAe);
% Reyf = rt.*sin(TAe);
% 
% [finalx,finaly] = circularOrbit(Gearth, sqrt(Rexf(i)^2+Reyf(i)^2));
% figure(2);
% plot(Rx1,Ry1,Rx2,Ry2,Rx,Ry,Rexf,Reyf,finalx,finaly)

% Simulate transfer orbit
% No errors in thrust vector
k=Gearth;
w=sqrt(k/(r1^3));
x0=[r1;0;0;w+(dv1/r1)]; % initial conditions



% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
opts = odeset('RelTol',1e-32,'AbsTol',1e-32);
[t,y] = ode45(@(t,y) myfunc(t,y), tspan, x0, opts);

% final orbit (burn 2)
tspan = linspace(0,Tf,1000);
xf = [y(1000,1), y(1000,2),y(1000,3),y(1000,4)+(dv2/y(1000,1))];
[tf,yf] = ode45(@(t,y) myfunc(t,y), tspan, xf, opts);




% %%% TRY TO START W/ MULTIPLE XO'S
% tspan = linspace(0,Tf,1000);
% X0 = [x0 ; xf];
% [tf,y999] = ode45(@(t,y) myfunc(t,y,k), tspan, xf, opts);
% 


figure(3);
s=stackedplot(t, y);
xlabel('Time t (s)');
s.DisplayLabels={'Radial Position (km)', 'Radial Velocity (km/sec)', 'Angular Position (rad)', 'Angular Velocity (rad/sec)'};
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Ideal Nonlinear State Change After Burn 1 vs. Time (ode45)');

figure(4);
plot(y(:,1).*sin(y(:,3)),y(:,1).*cos(y(:,3)),yf(:,1).*sin(yf(:,3)),yf(:,1).*cos(yf(:,3)))
title('Transfer and Final Orbit')

X0 = [x0; xf]
tspan = linspace(0,Tf,1000);
xf = [y(1000,1), y(1000,2),y(1000,3),y(1000,4)+(dv2/y(1000,1))];
[tf,y9999] = ode45(@(t,y) myfunc(y,u), tspan, xf, opts);


% Simulate ERROR CONDITIONS

% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
x1e=[r1;dv1*sind(error);0;w+(dv1*cosd(error)/r1)]; % initial conditions
[t,y1e] = ode45(@(t,y) myfunc(y,u), tspan, x1e, opts);

% final orbit (burn2)
tspan = linspace(0,Tf,1000);
x2e=[y1e(1000,1);y1e(1000,2) + (dv2*sind(error));y1e(1000,3);y(1000,4)+(dv2*cosd(error)/y1e(1000,1))]; % initial conditions
[t,y2e] = ode45(@(t,y) myfunc(y,u), tspan, x2e, opts);

figure(6);

plot(Rx1',Ry1',y(:,1).*sin(y(:,3)),y(:,1).*cos(y(:,3)),yf(:,1).*sin(yf(:,3)),yf(:,1).*cos(yf(:,3)),y1e(:,1).*sin(y1e(:,3)),y1e(:,1).*cos(y1e(:,3)),y2e(:,1).*sin(y2e(:,3)),y2e(:,1).*cos(y2e(:,3)))
title('Hohmann Transfer with and without thrust error')
legend("initial orbit", "transfer orbit", "final orbit", "transfer orbit w/error", "final orbit w/error")

figure(7);
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
% 
% 
% u = [0;0]
% 
% % Define nonlinear model as anonymous function
% % 
% % Define wrapper function for use in ODE45
% myfunc = @(t, x) xdot(x,u);
% % Compute total state using ODE45
% options = odeset('RelTol', 1e-12); % Set ODE45 options
% [tout, xnl] = ode45(myfunc, t, X0, options);
% xnl = xnl'; % Take transpose so that dimensions match with xnom


%% Sensitivity study
% finite difference dy/dk (ypert-ynom)/(kpert-knom) 

% k = pointing error
% y = position on arrival of final orbit
dk = error;
dy = yf-y2e;

sens = dy/dk;
E = mean(sens);
error = deg2rad(0.1);  % -0.3

% transfer orbit (burn 1)
tspan = linspace(0,Tof,1000);
x1ee=[r1;dv1*sind(error);0;w+(dv1*cosd(error)/r1)]; % initial conditions
[t,y1ee] = ode45(@(t,y) myfunc(y,u), tspan, x1ee, opts);

% final orbit (burn2)
tspan = linspace(0,Tf,1000);
x2ee=[y1e(1000,1);y1e(1000,2) + (dv2*sind(error));y1e(1000,3);y(1000,4)+(dv2*cosd(error)/y1e(1000,1))]; % initial conditions
[t,y2ee] = ode45(@(t,y) myfunc(y,u), tspan, x2ee, opts);

YY = mean(y2ee-yf)

Vi = sqrt(var(y2e))
Ve3 = sqrt(var(y2ee))


%% functions
function dy=myfunc(t,y,u)
init;
a = Gearth;
    dy=zeros(4,1);
    dy(1)=y(2);
    dy(2)=y(1).*(y(4)^2)-(a./(y(1)^2))+u(:,1);
    dy(3)=y(4);
    dy(4)=(-2.*y(2).*(y(4))./y(1))+u(:,2);
end

function [x,y] = circularOrbit(mu,r)
init    
    T = PeriodfromA(mu,r);
    t = linspace(0,T,1000);
    n = MMfromP(T);
    x = r * sin(n.*t);
    y = r * cos(n.*t);
end

% 
% 
% % Define nonlinear model as anonymous function
% pert = [0, 0.075, 0, -0.021]'; % perturbation 
% xdot = @(x, u) [x(2);
%     -mu*x(1)/sqrt((x(1)^2+x(3)^2)^3) + u(1);
%     x(4);
%     -mu*x(3)/sqrt((x(1)^2+x(3)^2)^3) + u(2)];
% % Define initial total state
% X0 = [x0, x0dot, y0, y0dot]' +  pert; 
% % 
% % Define wrapper function for use in ODE45
% myfunc = @(t, x) xdot(x, [0; 0]);
% % Compute total state using ODE45
% options = odeset('RelTol', 1e-12); % Set ODE45 options
% [tout, xnl] = ode45(myfunc, t, X0, options);
% xnl = xnl'; % Take transpose so that dimensions match with xnom
% 
% % % linearized nominal simulation
% % Define nominal state and control input over time
% dx = zeros(4, length(t)); % Perturbation state vector
% theta = omega * t;
% xnom = r0*cos(omega * t);
% ynom = r0*sin(omega * t);
% xdotnom = -r0 * omega * sin(omega * t);
% ydotnom = r0 * omega * cos(omega * t);
% 
% u = zeros(2, length(t));
% % Initialize perturbation state
% dx(:, 1) = [0, 0.075, 0, -0.021]';
% % Compute perturbation states over time
% for i = 2:length(t)
%     dx(:, i) = F*dx(:, i - 1) + G*u(:, i); % need to pass in r and t to Atilde...
% end
% 
% % Compute total states over time
% x = xnom + dx;
