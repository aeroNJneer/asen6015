% 
% 
% compute_B_matrix()
% 
% function B = compute_B_matrix(x_d, u)
%     % compute_B_matrix - Computes the B matrix by linearizing around the desired state x_d
%     %
%     % Inputs:
%     %   x_d - Desired state vector [a, e, i, RAAN, Argument of periapsis, MeanAnomaly]
%     %   u - Control input vector [a_r, a_theta, a_h]
%     %
%     % Output:
%     %   B - Jacobian matrix of the system with respect to the control inputs
% 
%     % Define symbolic variables for the state and control inputs
    syms a e i RAAN omega M a_r a_theta a_h  p TA r h mu

    p = a*(1-e^2);
    h = sqrt(mu*p);
    % State vector
    x = [a; e; i; RAAN; omega; M];

    % Control input vector
    u = [a_r; a_theta; a_h];

    % Define the equations of motion (Gauss' Variational equations)
    % Note: These equations need to be defined based on the specific system dynamics
    % For demonstration, we use placeholder equations
    f = [a+(2*a^2*e*sin(TA)/h)*a_r + 2*a^2*p/(h*r)*a_theta;
        e + (p*sin(TA)/h)*a_r + (((p+r)*cos(TA))/h)*a_theta + r*e*a_theta; 
        i + (r*cos(TA+omega)/h)*a_h; 
        RAAN + (r*sin(TA+omega)/(h*sin(i)))*a_h ; 
        omega - (p*cos(TA)/(h*e))*a_r + ((p+r)*sin(TA)/(h*e))*a_theta - (r*sin(TA+w)*cos(i)/(h*sin(i)))*a_h; 
        (M - (3/2)*sqrt(mu/(r^3))) + (sqrt(1-e^2)/(h*e)*(p*cos(TA+omega)-(2*r*e)))*a_r + (sqrt(1-e^2)/(h*e))*(p+r)*sin(TA)*a_theta];

    % Compute the Jacobian matrix of the system with respect to the control inputs
    B = jacobian(f, u)

    % Substitute the desired state into the Jacobian matrix
    B = subs(B, x, x_d)

    % Convert symbolic matrix to numeric matrix
%    B = double(B);
%end
