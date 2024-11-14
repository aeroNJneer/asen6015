% hwk 4

clear
clc
close all
init;const;
set(0,'DefaultAxesXGrid','on', ...
    'DefaultAxesYGrid','on', ...
    'DefaultLineLineWidth',2, ...
    'DefaultAxesLineWidth',1.5);
set(groot, 'defaultAxesFontSize',12, ...
    'defaultTextFontSize',12, ...
    'defaultLegendFontSize',14, ...
    'defaultAxesTitleFontSizeMultiplier', 1.5, ...
    'defaultAxesTitleFontSizeMultiplier', 1.2);

DCM313 = [ 0.8192 -0.5736 0;
         0.5736 0.8192  0 ;
         0 0 1 ];
mu = Gearth;
athrust = 0.03; %   km/s^2
tfguess =  23;

% function 
theta = @(Ax,Az,Bx,Bz,t) atan2(Az+Bz*t,Ax+Bx*t);
alpha = @(Ax,Ay,Bx,By,t) atan2(Ay+By*t,Ax+Bx*t);
uhat = @(A,B,t) (A + B*t)/norm(A+B*t);

% LVLH is aligned with desired orbital plane 
rhat = @(r) -r/norm(r);                          % ECI to LVLH
hhat = @(r,v) -cross(r,v)/norm(cross(r,v));      % ECI to LVLH @ tf
vhat = @(y,z) cross(y,z)/norm(cross(y,z));       % ECI to LVLH

lvlh2eci = @(r,h,v) [r,h,v];

TA = @(r) atan2(dot(r,Vperf_lvlh),dot(r,Rperf_lvlh));
Rreq = @(r) pf/(1+eccf*cos(TA(r)));
Vreq = @(r) -1/rperf * sqrt(mu/pf) * sin(TA(r)*Rperf_lvlh) + (1 - rperf/pf*(1-cos(TA(r)))*Vperf_lvlh);

% [ai,ecci,inci,Raani,wi,true_anomaly] = ECI2classical(r0,v0,Gearth);
% rperi = Rp(ai,ecci);                % initial radius of periapsis
% vperi = VfromRA(mu,rperi,ai);       % initial velocity at periapsis
% a_delta = af - ai;                  % delta SMA
% ecc_delta = eccf - ecci;            % delta eccentricity
% p_i = PfromAe(ai,ecci);             % initial semi-latus rectum
% % w_delta = wf - wi;
% h_delta = hf - norm(hi);            % delta momentum
% vp_delta = vperf-vperi;             % delta velocity
% rp_delta = rperf - rperi;           % delta radius of periapsis

r = r0;
v = v0;
delu = 0.1;
ts = linspace(0,tfguess,100);
% for i = 1:length(ts)
%     A_eci = lvlh2eci(zhat(r), yhat(r,v), xhat(yhat(r,v),zhat(r)))  * A0;
%     lamba_v = A_eci + B_eci*t;
%     
% end

% Solve the ODE
[t, state] = ode45(@(t, y) satellite_dynamics(t, y, A0, B0, athrust, mu), ts, X0);

figure(1)  % no correction ... 
plot3(state(:,1),state(:,2),state(:,3))

% now calculate differences 
% calculate difference from predicted to desired
rpred = state(end,1:3)
vpred = state(end,4:6)
pred_nu = TA(rpred)
pred_rreq_lvlh = Rreq(rpred)     
pred_vreq_lvlh = Vreq(vpred)

pred_rreq_eci = DCM313 * pred_rreq_lvlh

% Function to compute satellite dynamics
function dydt = satellite_dynamics(t, y, A, B, burn_acceleration, mu)
init;
    r = y(1:3);
    v = y(4:6);
    z = zhat(r);
    y = yhat(r,v);
    x = xhat(y,z);

    C = lvlh2eci(z,y,x);
    u_hat = (A+B)/(norm(A)+norm(B));
   
    r_norm = norm(r);
    a_gravity = -mu / r_norm^3 * r;

    % burn will be in lvlh frame: u_hat
    a_burn_lvlh = burn_acceleration * u_hat;
    
    a_burn = C * a_burn_lvlh;
    dydt = [v; a_gravity + a_burn];
end

function [E, J] = compute_error_and_jacobian(u)
    % Compute error vector E and Jacobian matrix J
    E = compute_error(u);
    J = compute_jacobian(A,B,t);
end

function J = compute_jacobian(lambda)
    % Compute the Jacobian matrix J using finite differencing
    delta = [1e-2;1e-2;1e-2;1e-3;1e-3;1e-3];
    
    lambda = [A+B*t; B];
    lamdapert = [A + delA +B*t; B+0.01];
    for i = 1:length(lambda)
        dellambda(i) = lambda(i)/(lambdapert(i)-lambda(i));
    end

    for i = 1:n
        lambda_perturb = lambda;
        lambda_perturb(i) = lambda_perturb(i) + delta(i);
        E_perturb = compute_error(lambda_perturb);
        dellambda = lambda(i)/(lambda_perturb(i)-lambda(i));

        J(:, i) = (E_perturb - compute_error(lambda_perturb(i))) / dellambda;
    end
end

function predictor_corrector_guidance()
    % Initial conditions
    u0 = [initial_guess]; % Initial guess for control vector
    tol = 1e-6; % Tolerance for convergence
    max_iter = 100; % Maximum number of iterations

    % Predictor-Corrector loop
    for k = 1:max_iter
        % Prediction step
        [E, J] = compute_error_and_jacobian(u0);

        % Check for convergence
        if norm(E) < tol
            disp('Converged');
            break;
        end

        % Correction step
        delta_u = -J \ E;
        u0 = u0 + delta_u;
    end

    % Final control vector
    u_final = u0;
    disp('Final control vector:');
    disp(u_final);
end



function E = compute_error(lambda)
    % final conditions
    const;
    TA = @(r) atan2(dot(r,Vperf_lvlh),dot(r,Rperf_lvlh));
    Rreq = @(r) pf/(1+eccf*cos(TA(r)));
    Vreq = -1/rperf * sqrt(mu/pf) * sin(TA(r)*Rperf_lvlh) + (1 - rperf/pf*(1-cos(TA(r)))*Vperf_lvlh);
    G = -mu * r / norm(r)^3;

    % Simulate the orbit with the current control vector u
    [t, state] = ode45(@(t, y) satellite_dynamics(t, y, A0, B0, athrust, mu), ts, X0);

    % Compute the error vector E
    E = [vreq(1) - v(1);
        vreq(2) - v(2);
        vreq(3) - v(3);
        rreq - norm(r);
        dot(Hf,r);
        lambda * G - dot(B,v);
        sqrt(lambaf*lambdaf) - 1;
    ];
end

