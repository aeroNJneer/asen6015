%initial conditions
x0 = -1.0;
y0 = -0.25
xdot0 = 0.003;
ydot0 = -0.0005;
mu = 3.986004418e5;

% final conditions
xf = -0.25 ;
yf = 0. ;
xdotf = 0.0002 ;
ydotf = 0. ;
Xf = [xf; yf; xdotf; ydotf];

r0 = [x0;y0;0];
v0 = [xdot0;ydot0];
rf = [xf;yf;0];
vf = [xdotf;ydotf];


[v1,v2] = lambert_solver(r0,rf,20*60,mu)

function [v1, v2] = lambert_solver(r1, r2, tof, mu)
    % r1 and r2 are the position vectors at the initial and final points
    % tof is the time of flight
    % mu is the gravitational parameter
    
    % Calculate the norms
    r1_norm = norm(r1);
    r2_norm = norm(r2);
    
    % Angle between the two vectors
    cos_theta = dot(r1, r2) / (r1_norm * r2_norm);
    
    % Calculate the chord length
    c = sqrt(r1_norm^2 + r2_norm^2 - 2 * r1_norm * r2_norm * cos_theta);
    
    % Calculate the semi-perimeter
    s = (r1_norm + r2_norm + c) / 2;
    
    % Calculate the universal anomaly
    % This involves solving Kepler's equation iteratively
    % Let's assume an initial guess for the universal anomaly
    chi = sqrt(mu) * tof / s;
    tol = 1e-6; % Tolerance for the iteration
    
    % Start iteration
    while true
        z = chi^2;
        C = stumpff_C(z);
        S = stumpff_S(z);
        
        % Update the guess
        chi_new = sqrt(s^3 / (8*mu)) * (1 - (r1_norm + r2_norm) / (2*s)) + ...
                  sqrt(mu / s^3) * chi^2 * (S * c + 0.5 * chi^2 * C) / c;
        
        if abs(chi_new - chi) < tol
            chi = chi_new;
            break;
        end
        chi = chi_new;
    end
    
    % Compute velocity vectors
    f = 1 - (chi^2 / r1_norm) * C;
    g_dot = 1 - (chi^2 / r2_norm) * C;
    g = tof - chi^3 / sqrt(mu) * S;
    
    v1 = (r2 - f * r1) / g;
    v2 = (g_dot * r2 - r1) / g;
end

function C = stumpff_C(z)
    if z > 0
        C = (1 - cos(sqrt(z))) / z;
    elseif z < 0
        C = (1 - cosh(sqrt(-z))) / z;
    else
        C = 1 / 2;
    end
end

function S = stumpff_S(z)
    if z > 0
        S = (sqrt(z) - sin(sqrt(z))) / (sqrt(z))^3;
    elseif z < 0
        S = (sinh(sqrt(-z)) - sqrt(-z)) / (sqrt(-z))^3;
    else
        S = 1 / 6;
    end
end


