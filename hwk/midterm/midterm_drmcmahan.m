% midterm - Dr McMahon

clear
close all
init;
global mu
mu = Gearth;
a = 6728;
n = sqrt(mu/a^3);
x0 = [-1.5 -0.5 .003 0];
xf = [-.05 0 0 0];
global A
A = [0 0 1 0;
    0 0 0 1;
    3*n^2 0 0 2*n;
    0 0 -2*n 0];
global B
B = [0 0;
    0 0;
    1 0;
    0 1];

% natural dynamics
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[T,X] = ode45(@derivs,[0 3*60*60],x0,options);

for i = 1:length(T)
    r(i) = sqrt(X(i,1)^2 + X(i,2)^2);
end

[minr, minrind] = min(r);
global Atilde;
Atilde = [A -B*B'; zeros(4) -A'];
Phi0 = eye(8);
Phi0 = reshape(Phi0,1,64);

% STM
[Tp, Phi] = ode45(@derivsPhi, [0 20*60], Phi0,options);
Phi_tf = reshape(Phi(end,:),8,8);
lambda0 = inv(Phi_tf(1:4,5:8)) * (xf' - Phi_tf(1:4,1:4)*x0');

% optimal trajectory
[Taug,Xaug] = ode45(@derivsAug, [0 20*60], [x0 lambda0'], options);

for i = 1:length(Taug)
    uopt(i,:) = -B * Xaug(i,5:8);
end

[Tcont, Xcont] = ode45(@derivs, [0 3*60*60], Xaug(end,1:4),options);

% % neighboring optimal case
% x0pert = [-0.98 -0.4 3.7e-3 0.3e-3];
% [Tpert,Xpert] = ode45(@derivspert,[0 20*60],x0pert,options);
% 
% delxf = zeros(4,1);
% delx0 = (x0pert - x0);
% dellambda0 = inv(Phi_tf(1:4,5:8))* (delxf - Phi_tf(1:4,1:4)*delx0);
% [Tneigh,Xneigh] = ode45(@derivsAug,Taug,[x0pert (lambda0+delLambda0)'],options);
% 
% for i = length(Tneigh)
%     uneighbor(i,:) = -B'*Xneigh(i,5:8);
% end

delu = uneighbor - uopt;

function dx = derivs(t,x)
    global A;global B;
    u = [0;0];
    dx = A*x +B*u;
end

function dx = derivsPhi(t,x)
    global Atilde;
    phi = reshape(x,8,8);
    dphi = Atilde*phi;
    dx = reshape(dphi,64,1);
end

function dx = derivsAug(t,x)
    global Atilde;
    dx = Atilde*x;
end

% function dx = derivspert(t,x)
%     u1 = interpl(Taug,uopt(:,1),t);
%     u2 = interpl(Taug,uopt(:,2),t);
