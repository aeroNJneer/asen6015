
diags = 1.2e-10;
offdiags = -0.1*1e-10;
Qguess = [diags, offdiags; offdiags, diags];

x_EKF = zeros(n, T, n_trials);
P_EKF = zeros(n, n, T, n_trials);
y_EKF = zeros(p, nS, T, n_trials);
S_EKF = zeros(p*nS, p*nS, T, n_trials);

Fx = @(x) eye(n)+Dt*A(atan2(x(3), x(1)), sqrt(x(1)^2+x(3)^2));
Hx = @(k, x, i) C(k*Dt, x, i);

fprintf('Running EKF over MC Simulations: ')
tic
for i = 1:n_trials
    [x_t, P_t, y_m, S_t] = EKF(...
        n, y_TMT(:, :, :, i), Qguess, 0.9*Rtrue, mu0, P0, Fx, Omegak,...
        Hx, Dt, rho, drho, phi, mu);
    x_EKF(:, :, i) = x_t;
    P_EKF(:, :, :, i) = P_t;
    y_EKF(:, :, :, i) = y_m;
    S_EKF(:, :, :, i) = S_t;
end
toc

e_x_EKF = x_TMT - x_EKF;
e_y_EKF = y_TMT - y_EKF;

NEES_EKF = zeros(T, n_trials);
NIS_EKF = NaN * zeros(T, n_trials);

fprintf('Running NEES/NIS over EKF Data: ')
tic
for i = 1:n_trials
    e_x_dat = e_x_EKF(:, :, i);
    P_dat = P_EKF(:, :, :, i);
    e_y_dat = e_y_EKF(:, :, :, i);
    S_dat = S_EKF(:, :, :, i);
    for j = 1:T
        
        % Run NEES Calculation
        e_x = e_x_dat(:, j);
        eps_x = e_x' * (P_dat(:, :, j) \ e_x);
        NEES_EKF(j, i) = eps_x;
        
        % Run NIS Calculation
        S_i = S_dat(:, :, j);
        S_inds = find(~isnan(S_i(1, :)), 1, 'last');
        S_i = S_i(1:S_inds, 1:S_inds);
        e_y = e_y_dat(:, :, j);
        e_y = e_y(:);
        e_y = e_y(~isnan(e_y));
        [r, ~] = size(S_i);
        eps_y = NaN;
        if r == length(e_y)
            e_y_norm_one_meas = S_i \ e_y;
            if r > 3
                eps_y = e_y(1:3)'*e_y_norm_one_meas(1:3);
            else
                eps_y = e_y' * e_y_norm_one_meas; 
            end 
        end
        NIS_EKF(j, i) = eps_y;
        
    end
end