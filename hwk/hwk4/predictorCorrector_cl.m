clear all
clc

    global mu acc
    mu = 3.986004415e5; % Earth's gravitational parameter [km^3/s^2]
    acc = 0.03; % Thrust acceleration [km/s^2]
    
    % Initial state
    x0 = [8276; 5612; 5; -3.142; 4.672; 0];
    tf_guess = 23;
    % Target orbit parameters
    target.a = 10000;  % km
    target.i = 0;      % degrees
    target.e = 0.0001;
    target.w = 35 * pi/180;  % radians
    [tf, A, B] = predictorCorrector(x0,tf_guess,acc,target)

function [tf, A_final, B_final] = predictorCorrector(x0, tf_guess, acc,target)
% Initial conditions
global mu
r0 = x0(1:3);
v0 = x0(4:6);

% Constants
mu = 398600.4418; % Earth's gravitational parameter (km^3/s^2)
tolerance = 1e-6;  % Convergence tolerance
max_iterations = 50;
dt = 0.1;  % Integration time step

% Initial guess for control parameters
A = [1; 0; 0];  % Initial A guess in LVLH
B = [0; 0; 0];  % Initial B guess in LVLH
tf = tf_guess;
u = [A;B;tf];



% Main predictor-corrector loop
for iter = 1:max_iterations
    % Predict step - Integrate nominal trajectory
    [t, x, lambda] = integrateTrajectory(r0, v0, u, acc, dt);
    
    % Calculate error vector
    E = calculateError(x(:,end), u, lambda(:,end), tf, target);
    E0_norm = norm(E);
    
    if E0_norm < tolerance
        break;
    end
    
    % Calculate Jacobian through finite differencing
    J = calculateJacobian(r0, v0, A, B, tf, acc, dt, E, target);
    
    % Update control parameters using line search
    alpha = lineSearch(r0, v0, A, B, tf, acc, dt, E, J, target);
    
    % Apply updates with line search parameter
    delta = -alpha * pinv(J) * E;
    
    % Extract updates for A, B, and tf
    dA = delta(1:3);
    dB = delta(4:6);
    dtf = delta(7);
    
    % Update parameters
    A = A + dA;
    B = B + dB;
    tf = tf + dtf;
end

A_final = A;
B_final = B;
end

function [t, x, lambda] = integrateTrajectory(r0, v0, u, acc, dt)
% Initialize state and time vectors
global mu
A = u(1:3);
B = u(4:6);
tf = u(7);
x0 = [r0; v0];
t = 0:dt:tf;
N = length(t);
x = zeros(6, N);
lambda = zeros(3, N);
x(:,1) = x0;

% RK4 integration
for i = 1:N-1
    k1 = dynamics(t(i), x(:,i), A, B, tf, acc);
    k2 = dynamics(t(i) + dt/2, x(:,i) + dt*k1/2, A, B, tf, acc);
    k3 = dynamics(t(i) + dt/2, x(:,i) + dt*k2/2, A, B, tf, acc);
    k4 = dynamics(t(i) + dt, x(:,i) + dt*k3, A, B, tf, acc);
    
    x(:,i+1) = x(:,i) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    
    % Calculate lambda (costate) at each step
    lambda(:,i) = B + (A + B*t(i))/norm(A + B*t(i));
end
lambda(:,N) = B + (A + B*tf)/norm(A + B*tf);
end

function dxdt = dynamics(t, state, A, B, tf, acc)
% Unpack state
global mu
r = state(1:3);
v = state(4:6);

% Calculate LVLH frame
h = cross(r, v);
rhat = r/norm(r);
hhat = h/norm(h);
that = cross(hhat, rhat);

% Current control direction in LVLH
u_lvlh = (A + B*t)/norm(A + B*t);

% Convert control to ECI
R_lvlh2eci = [rhat that hhat];
u_eci = R_lvlh2eci * u_lvlh;

% Gravitational acceleration
g = -mu*r/norm(r)^3;

% State derivatives
dr = v;
dv = g + acc*u_eci;

dxdt = [dr; dv];
end

function E = calculateError(xf, u, lambda_f, tf, target)
% Calculate orbital elements from final state
global mu
B = u(4:6);
[a, e, i, ~, omega, nu] = rv2oe(xf(1:3), xf(4:6));

% Calculate target orbit vectors in ECI
[rp, vp] = getTargetPeriapseVectors(target);

% Calculate error vector components using equations from homework
E = zeros(7,1);
E(1:3) = vp - xf(4:6);  % Velocity error ECI
E(4) = norm(xf(1:3)) - target.a*(1-target.e^2)/(1+target.e*cos(nu));  % Radius error
E(5) = dot(cross(rp, vp)/norm(cross(rp, vp)), xf(1:3));  % Orbit plane error
E(6) = dot(lambda_f, -mu*xf(1:3)/norm(xf(1:3))^3) - dot(B, xf(4:6));  % Transversality
E(7) = norm(lambda_f) - 1;  % Normalization
end

function J = calculateJacobian(r0, v0, A, B, tf, acc, dt, E0, target)
% Calculate Jacobian through forward finite differencing
global mu
epsilonA = 0.1;
epsilonB = 0.01;
epsilont = 0.1;

n_params = 7;  % 3 for A, 3 for B, 1 for tf
J = zeros(length(E0), n_params);

% Perturb A components
for i = 1:3
    A_perturbed = A;
    A_perturbed(i) = A(i) + epsilonA;
    u_perturbed = [A_perturbed; B; tf];
    [~, xf, lambda_f] = integrateTrajectory(r0, v0, u_perturbed, acc, dt);
    E_perturbed = calculateError(xf(:,end), u_perturbed, lambda_f(:,end), tf, target);
    J(:,i) = (E_perturbed - E0)/epsilonA;
end

% Perturb B components
for i = 1:3
    B_perturbed = B;
    B_perturbed(i) = B(i) + epsilonB;
    u_perturbed = [A; B_perturbed; tf];
    [~, xf, lambda_f] = integrateTrajectory(r0, v0, u_perturbed, acc, dt);
    E_perturbed = calculateError(xf(:,end), u_perturbed, lambda_f(:,end), tf, target);
    J(:,i+3) = (E_perturbed - E0)/epsilonB;
end

% Perturb tf
u_perturbed = [A; B; tf+epsilont];
[~, xf, lambda_f] = integrateTrajectory(r0, v0, u_perturbed, acc, dt);
E_perturbed = calculateError(xf(:,end), u_perturbed, lambda_f(:,end), tf + epsilont, target);
J(:,7) = (E_perturbed - E0)/epsilont;
end

function alpha = lineSearch(r0, v0, A, B, tf, acc, dt, E0, J, target)
% Line search to find step size that reduces error
global mu
alpha = 1.0;
alpha_min = 1e-6;
rho = 0.5;  % Reduction factor
c = 0.1;    % Sufficient decrease parameter

E0_norm = norm(E0);
delta = -(pinv(J) * E0);

while alpha > alpha_min
    % Try update with current alpha
    A_new = A - alpha*delta(1:3);
    B_new = B - alpha*delta(4:6);
    tf_new = tf - alpha*delta(7);
    u_new = [A_new; B_new; tf_new];
    % Check new error
    [~, xf, lambda_f] = integrateTrajectory(r0, v0, u_new, acc, dt);
    E_new = calculateError(xf(:,end), u_new, lambda_f(:,end), tf_new, target);
    
    if norm(E_new) < (1 - c*alpha)*E0_norm
        break;
    end
    
    alpha = rho*alpha;
end

if alpha <= alpha_min
    alpha = alpha_min;
end
end

function [a, e, i, Omega, omega, nu] = rv2oe(r, v)
% Convert state vectors to orbital elements
global mu
h = cross(r, v);
n = cross([0; 0; 1], h);
e_vec = cross(v, h)/mu - r/norm(r);
e = norm(e_vec);
a = 1/(2/norm(r) - dot(v,v)/mu);

i = acos(h(3)/norm(h));
Omega = acos(n(1)/norm(n));
if n(2) < 0
    Omega = 2*pi - Omega;
end

omega = acos(dot(n,e_vec)/(norm(n)*e));
if e_vec(3) < 0
    omega = 2*pi - omega;
end

nu = acos(dot(e_vec,r)/(e*norm(r)));
if dot(r,v) < 0
    nu = 2*pi - nu;
end
end

function [rp, vp] = getTargetPeriapseVectors(target)
% Calculate periapse position and velocity vectors for target orbit
global mu
p = target.a*(1-target.e^2);
rp_mag = p/(1+target.e);
vp_mag = sqrt(mu*(2/rp_mag - 1/target.a));

% Calculate in perifocal frame
rp_pqw = rp_mag*[1; 0; 0];
vp_pqw = vp_mag*[0; 1; 0];

% Rotation matrices
R3_omega = [cos(target.w) -sin(target.w) 0;
            sin(target.w) cos(target.w)  0;
            0         0          1];
R1_i = [1  0        0;
        0  cos(target.i)  -sin(target.i);
        0  sin(target.i)   cos(target.i)];

% Transform to ECI
R_pqw2eci = R3_omega*R1_i;
rp = R_pqw2eci*rp_pqw;
vp = R_pqw2eci*vp_pqw;
end