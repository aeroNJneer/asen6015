init
mu = Gearth;
a = 7000;
e = 0.1;
i = deg2rad(1);
w = 0;
raan = deg2rad(1);
TA = 0;
p = PfromAe(a,e);
h = HfromP(mu,p);
r = RfromP(p,e,TA);
k = 0.1;
M = 0;
<<<<<<< Updated upstream
u_rv = [0 0; 1 0; 0 0 ; 0 1];
u_coe = 

X = [a; e; i; raan; w; n];
=======
X = [a; e; i; raan; w; M];
>>>>>>> Stashed changes
Xd = [a+100; e+0.01; i+deg2rad(5); raan + deg2rad(10); w+deg2rad(15); M+deg2rad(20)]
A = zeros(6,6);
A(6,1) = -(3/2)*sqrt(mu/(a^5));

B = [2*a^2*e*sin(TA)/h 2*a^2*p/h/r 0;
    p*sin(TA)/h (p+r)*cos(TA)+(r*e) 0;
    0 0 r*cos(w+TA)/h;
    0 0 r*sin(w+TA)/(h*sin(i));
    -p*cos(TA)/(h*e) (p+r)*sin(TA)/(h*e) -r*sin(w+TA)*cos(i)/(h*sin(i));
    sqrt(1-e^2)*p*cos(TA-2*r*e)/(h*e) -sqrt(1-e^2)*(p+r)*sin(TA)/(h*e) 0];

K = k*eye(3);

% closed loop
X_dot = (A-B'*K)*X + B*K*Xd



function B = compute_B_matrix(x_d, u)
    % compute_B_matrix - Computes the B matrix by linearizing around the desired state x_d
    %
    % Inputs:
    %   x_d - Desired state vector [a, e, i, RAAN, Argument of periapsis, MeanAnomaly]
    %   u - Control input vector [a_r, a_theta, a_h]
    %
    % Output:
    %   B - Jacobian matrix of the system with respect to the control inputs

    % Define symbolic variables for the state and control inputs
    syms a e i RAAN omega M a_r a_theta a_h real

    % State vector
    x = [a; e; i; RAAN; omega; M];

    % Control input vector
    u = [a_r; a_theta; a_h];

    % Define the equations of motion (Gauss' Variational equations)
    % Note: These equations need to be defined based on the specific system dynamics
    % For demonstration, we use placeholder equations
    f = [a + a_r; e + a_theta; i + a_h; RAAN; omega; M];

    % Compute the Jacobian matrix of the system with respect to the control inputs
    B = jacobian(f, u);

    % Substitute the desired state into the Jacobian matrix
    B = subs(B, x, x_d);

    % Convert symbolic matrix to numeric matrix
    B = double(B);
end


