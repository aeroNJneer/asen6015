
clear
close all
init;
error = 1; % degree in r_hat direction
% problem 3 - Hohmann transfer
r1 = 200 + Rearth;
r2 = 35786;
mu = Gearth;
dt = 10;
t = 0:dt:80000;

% circular orbit velocity
v1 = sqrt(Gearth/r1);
v2 = sqrt(Gearth/r2);

x0 = r1;
y0 = 0;
r0 = sqrt(x0^2+y0^2);
x0dot = 0;
y0dot = r0 * sqrt(mu/r0^3);
omega = sqrt(mu/r0^3); % 

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



% Define nonlinear model as anonymous function
pert = [0, 0.075, 0, -0.021]'; % perturbation 
xdot = @(y, ur,ut) [y(2);
    y(1).*(y(4)^2)-(mu./(y(1)^2))+ur;
    y(4);
    (-2.*y(2).*(y(4))./y(1))+ut];

% INITIAL CONDITIONS CAN BE AN ARRAY!!



% simulate hohmann transfer trajectory

u = [0;0]

% Define nonlinear model as anonymous function
pert = [0, 0, 0, 0]'; % perturbation 
% xdot = @(x, u) [x(2);
%     -mu*x(1)/sqrt((x(1)^2+x(3)^2)^3) + u(1);
%     x(4);
%     -mu*x(3)/sqrt((x(1)^2+x(3)^2)^3) + u(2)];
% Define initial total state

k=Gearth;
w=sqrt(k/(r1^3));
X0=[r1;0;0;w+(dv1/r1)]; % initial conditions

% 
% Define wrapper function for use in ODE45
xdotwrap = @(t, x) xdot(x, u(1),u(2));
% Compute total state using ODE45
options = odeset('RelTol', 1e-12); % Set ODE45 options
[tout, xnl] = ode45(xdotwrap, t, X0, options);
xnl = xnl'; % Take transpose so that dimensions match with xnom

% % linearized nominal simulation
% Define nominal state and control input over time
dx = zeros(4, length(t)); % Perturbation state vector
theta = omega * t;
xnom = r0*cos(omega * t);
ynom = r0*sin(omega * t);
xdotnom = -r0 * omega * sin(omega * t);
ydotnom = r0 * omega * cos(omega * t);

u = zeros(2, length(t));
% Initialize perturbation state
dx(:, 1) = [0, 0.075, 0, -0.021]';
% Compute perturbation states over time
% for i = 2:length(t)
%     dx(:, i) = F*dx(:, i - 1) + G*u(:, i); % need to pass in r and t to Atilde...
% end

% Compute total states over time
x = xnom + dx;


% Plot
figure;
subplot(4, 1, 1);
plot(t, xnl(1, :));
title('Total States vs. Time, Linearized DT SS Model');
ylabel('x (km)'); %ylim([-7000, 7000])
grid on;
subplot(4, 1, 2);
plot(t, xnl(2, :));
ylabel('\deltax (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, xnl(3, :));
ylabel('y (km)');
grid on;
subplot(4, 1, 4);
plot(t, xnl(4, :));
ylabel('\deltay (km/s)');
xlabel('Time t (s)');
grid on;




%% functions
function dy=myfunc(t,y,ur,ut)
    dy=zeros(4,1);
    dy(1)=y(2);
    dy(2)=y(1).*(y(4)^2)-(a./(y(1)^2))+ur;
    dy(3)=y(4);
    dy(4)=(-2.*y(2).*(y(4))./y(1))+ut;
end

function [x,y] = circularOrbit(mu,r)
init    
    T = PeriodfromA(mu,r);
    t = linspace(0,T,1000);
    n = MMfromP(T);
    x = r * sin(n.*t);
    y = r * cos(n.*t);
end
