init;

g = @(r) Gearth/r^2;
n_c = @(Np, Vc,l_dot) Np*Vc*lambda_dot;
R = @(Rt,Rm) Rt-Rm;
V = @(Vt,Vm) Vt-Vm;

Vm = @(vm,theta_he,lambda) [vm*cos(L+theta_he + lambda);
    vm*sin(L+theta+lambda)];

Vt = @(vt,beta) [-vt*cos(beta);
    vt*sin(beta)];

a = @(at,am) at - am;        % no differential gravity
at = @(nt,beta) [nt*sin(beta)
    nt*cos(beta)];

beta_dot = @(nt,vt) nt/vt;

lambda_dot = @(R,V) (R(1)*V(2) - R(2)*V(1))/norm(R^2);     % LOS rate
theta_he = @(R,V) -acos(dot(R,V));
Herror = cos(theta_he(R,V));
L = @(vt,vm,beta,lambda) asin((vt*sin(beta+lambda))/vm);    % lead angle
A = [0 1 0;
    0 0 1;
    0 0 0];
B = [0 -1 0]';

vm = 3000;      % simulation constants ft/s
vt = 1000;      % simulation constants ft/s

Vc = Vm + Vt;   % closing speed ft/s
R0 = 40000;     % ft
y0 = 0;

tgo = R/Ve      % 


ydot; 
ydot0 ;         % determed by initial heading earror

x = [y ydot nt]';

theta_i = deg2rad(30);
nu_i = 0;

theta_ii = deg2rad(0);
nu_ii = 3*g;

theta_iii = deg2rad(30);
nu_iii = 3*g;

Nprime = [3 4 5];

Theta_he = [theta_i; theta_ii; theta_iii];

linearHeadingError = @(vm,theta_he,Nprime,t,tf);

HE0 = linearHeadingError(vm,theta_he,Nprime,0,100)
%% part a
%% proportional navigation



%% augmented proportional navigation

