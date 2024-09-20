init;

mu = Gearth;

a = Rearth + 35786;    % GEO altitude
e = 0.01;                % near circular orbit
i = 0;                  % inclination
Raan = 0;               % RAAN
w = 0;                  % argument of periapsis

n =  MMfromMua(mu,a);   % mean motion
M = @(t) n*t;           % mean anomaly
tof = PeriodfromA(mu,a); % orbital period
t = 0:60:tof;
X = @(t) [a e i Raan w M(t)];

X0 = [a e i Raan w M(0)]';
Xd0 = [a e i Raan w M(0)]';
Mt = M(t)';
Xd = @(t) [a e i Raan w M(t)];
A = @(t) X(t) .* eye(6);
A0 = X0 .* eye(6);
B = [0 0; 1 0; 0 0 ; 1 0];
k = eye(6);
Xnom = 0;
t1 = 0:0.1:1;
for i = 2:length(t)
    A = [a e i Raan w Mt(i)]'.*eye(6);
    %[t,x= = @(t,x) ode45(odefunc(t,x,X0, dx = (A-B*k)*X + k*Xd;
   [t,y] = ode45(@(t,y) odefunc(t,y,A,B,k,Xd(i-1)), t1, X0);
    
end 
figure(1)
plot(x(1,:),x(3,:));
function dx = odefunc(t,x,A,B,k,Xd)
    dx = (A-B*k)*X + k*Xd;
end