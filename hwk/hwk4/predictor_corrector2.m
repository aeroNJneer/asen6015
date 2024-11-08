function [A, B, t_burn, nu, state_history, lambda_r, lambda_v] = predictor_corrector2(x0, lambda_r0, lambda_v0, A, B, thrust_acc, t_burn_guess)
% Constants
mu = 3.986e14;  % Earth's gravitational parameter (m^3/s^2)

% Target orbit parameters
a_target = 10000e3;  % semi-major axis (m)
i_target = 0 * pi/180;  % inclination (rad)
e_target = 0.0001;  % eccentricity
omega_target = 35 * pi/180;  % argument of periapsis (rad)

% Integrate equations of motion with initial guess
t_span = linspace(0, t_burn_guess, 100);
[t, state] = ode45(@(t, x) two_body_ode(t, x, A, B, lambda_r0, lambda_v0, thrust_acc), t_span, [x0; lambda_r0; lambda_v0]);
r_final = state(end, 1:3);
v_final = state(end, 4:6);
lambda_r_final = state(end, 7:9);
lambda_v_final = state(end, 10:12);

% Compute target orbit parameters
h_target = cross(r_final, v_final);
p_target = norm(h_target)^2 / mu;
e_target_vec = (norm(v_final)^2 / mu - 1 / norm(r_final)) * r_final - v_final / norm(v_final) * sqrt(norm(v_final)^2 * p_target / mu);
e_target = norm(e_target_vec);
i_target = acos(h_target(3) / norm(h_target));
omega_target = atan2(e_target_vec(2), e_target_vec(1));

% Compute error quantities
v_r1 = v_final(1);
v_r2 = v_final(2);
v_r3 = v_final(3);
r_r = norm(r_final);
h_r_hat = h_target / norm(h_target);
lambda_t = A + B * t_burn_guess;
g_t = -mu * r_final / r_r^3;

error1 = v_r1 - v_final(1);
error2 = v_r2 - v_final(2);
error3 = v_r3 - v_final(3);
error4 = r_r - a_target / (1 - e_target^2);
error5 = h_r_hat.' * r_final;
error6 = lambda_t.' * g_t - B.' * v_final;
error7 = sqrt(lambda_t.' * lambda_t) - 1;

% Compute true anomaly
nu = atan2(r_final.' * v_final / sqrt(mu * p_target), r_r - p_target / r_r);

% Iteratively update guidance parameters
tol = 1e-6;
max_iter = 20;
for i = 1:max_iter
    x0 = [A; B; t_burn_guess];
    [x, fval, exitflag] = fsolve(@(x) error_func(x, v_r1, v_r2, v_r3, r_r, h_r_hat, g_t, v_final, lambda_r_final, lambda_v_final), x0, optimoptions('fsolve', 'TolFun', tol));
    
    A = x(1:3);
    B = x(4:6);
    t_burn = x(7);
    
    % Integrate equations of motion with updated guidance parameters
    t_span = linspace(0, t_burn, 100);
    [t, state] = ode45(@(t, x) two_body_ode(t, x, A, B, lambda_r_final, lambda_v_final, thrust_acc), t_span, [x0(1:6); lambda_r_final; lambda_v_final]);
    r_final = state(end, 1:3);
    v_final = state(end, 4:6);
    lambda_r_final = state(end, 7:9);
    lambda_v_final = state(end, 10:12);
    
    % Recompute error quantities
    v_r1 = v_final(1);
    v_r2 = v_final(2);
    v_r3 = v_final(3);
    r_r = norm(r_final);
    h_r_hat = h_target / norm(h_target);
    lambda_t = A + B * t_burn;
    g_t = -mu * r_final / r_r^3;
    
    error1 = v_r1 - v_final(1);
    error2 = v_r2 - v_final(2);
    error3 = v_r3 - v_final(3);
    error4 = r_r - a_target / (1 - e_target^2);
    error5 = h_r_hat.' * r_final;
    error6 = lambda_t.' * g_t - B.' * v_final;
    error7 = sqrt(lambda_t.' * lambda_t) - 1;
    
    nu = atan2(r_final.' * v_final / sqrt(mu * p_target), r_r - p_target / r_r);
    
    if max(abs([error1, error2, error3, error4, error5, error6, error7])) < tol
        break;
    end
end

state_history = [t, state];
end

function dx = two_body_ode(t, x, A, B, lambda_r, lambda_v, thrust_acc)
r = x(1:3);
v = x(4:6);
lambda_r = x(7:9);
lambda_v = x(10:12);

tan_theta = (A(3) + B(3)*t) / (A(1) + B(1)*t);
tan_alpha = (A(2) + B(2)*t) / (A(1) + B(1)*t);
u = (A + B*t) / norm(A + B*t);  % lvlh frame
a_thrust = thrust_acc * u;

% need DCM here ... LVLH to ECI
a_grav = -mu * r / norm(r)^3;   % ECI frame
a_thrust_eci = C * a_thrust_lvlh

dr_dt = v;
dv_dt = a_grav + a_thrust;
dlambda_r_dt = A + B*t;
dlambda_v_dt = B;

dx = [dr_dt; dv_dt; dlambda_r_dt; dlambda_v_dt];
end

function f = error_func(x, v_r1, v_r2, v_r3, r_r, h_r_hat, g_t, v_final, lambda_r, lambda_v)
A = x(1:3);
B = x(4:6);
t_burn = x(7);

lambda_t = A + B * t_burn;

f = [A(1) + B(1)*t_burn - v_r1;
     A(2) + B(2)*t_burn - v_r2;
     A(3) + B(3)*t_burn - v_r3;
     sqrt(A(1)^2 + A(2)^2 + A(3)^2) - r_r;
     A(1)*h_r_hat(1) + A(2)*h_r_hat(2) + A(3)*h_r_hat(3);
     lambda_t.' * g_t - B.' * v_final;
     sqrt(lambda_t.' * lambda_t) - 1];
end