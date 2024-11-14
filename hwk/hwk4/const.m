init;
mu = Gearth;

DCM313 = [ 0.8192 -0.5736 0;
         0.5736 0.8192  0 ;
         0 0 1 ];
athrust = 0.03; %   km/s^2
% initial guess
A0 = [1;0;0];
B0 = [0;0;0];
tfguess =  23;
lambdaguess = A0 + B0*tfguess;

% initial conditions
r0 = [8276; 5612; 5];
v0 = [-3.142; 4.672; 0];
X0 = [ r0; v0];
t0 = 0;
h_initial = cross(r0,v0);

% final conditions

a_final = 10000;
inc_final = 0;
ecc_final = 0.0001;
w_final = deg2rad(35);
h_final = HfromAe(mu,a_final,ecc_final);
H_final = [0; h_finalf; 0];                    % angular momentum vector in lvlh
p_final = PfromAe(a_final,ecc_final);
rad_peri_final = Rp(a_final,ecc_final);                 % radius of periapsis
Rad_peri_final_lvlh = [rad_peri_final; 0; 0];         % Rp vector in lvlh
vel_peri_final = VfromRA(mu,rad_peri_final,a_final);       % velocity at periapsis
Vel_peri_final_lvlh = [0; 0; vel_peri_final];         % Vp in lvlh
