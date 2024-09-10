init
mu = Gearth;
a = 7000;
e = 0.01;
i = 0;
w = 0;
raan = 0;
TA = 0;
p = PfromAe(a,e);
h = HfromP(mu,p);
r = RfromP(p,e,TA);
k = 0.1;
M = 0;
u_rv = [0 0; 1 0; 0 0 ; 0 1];
u_coe = 

X = [a; e; i; raan; w; n];
Xd = [a+100; e+0.01; i+deg2rad(5); raan + deg2rad(10); w+deg2rad(15); M+deg2rad(20)]
A = zeros(6,6);
A(6,1) = -(3/2)*sqrt(mu/(a^5));

B = [2*a^2*e*sin(TA)/h 2*a^2*p/h/r 0;
    p*sin(TA)/h (p+r)*cos(TA)+(r*e) 0;
    0 0 r*cos(w+TA)/h;
    0 0 r*sin(w+TA)/(h*sin(i));
    -p*cos(TA)/(h*e) (p+r)*sin(TA)/(h*e) -r*sin(w+TA)*cos(i)/(h*sin(i));
    sqrt(1-e^2)*p*cos(TA-2*r*e)/(h*e) -sqrt(1-e^2)*(p+r)*sin(TA)/(h*e) 0];

B = [B; zeros(3,2)]
K = k*eye(6);
%K = [k;k;k;k;k;k]

X_dot = (A-B'*K)*X + B*K*Xd


