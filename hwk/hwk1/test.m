clear
init
k = Gearth;
% problem 3 - Hohmann transfer
r1 = 200 + Rearth;
r2 = 35786;
% circular orbit velocity
v1 = sqrt(Gearth/r1);
v2 = sqrt(Gearth/r2);

% orbital elements of transfer

e = EccfromRpRa(r1,r2)
a = (r1+r2)/2
Tof = PeriodfromA(Gearth, a)/2  % time of flight
Ti = PeriodfromA(Gearth, r1);   % period of lower orbit
Tf = PeriodfromA(Gearth, r2);   % period of higher orbig
rp = r1;
ra = r2;

vt1 = VfromRA(Gearth,r1,a)
vt2 = VfromRA(Gearth,r2,a)

hi = HfromRe(Gearth,r1,0)   % momentum initial orbit
hf = HfromRe(Gearth,r2,0)   % momentum final orbit
ht = HfromAe(Gearth,a,e)    % momentum transfer orbit
% delta v impulse burns
dv1 = vt1-v1
dv2 = v2-vt2
r = r1;
w=sqrt(k/(r^3))

A=[0,1,0,0;3*w^2,0,0,2*r*w;0,0,0,1;0,-2*w/r,0,0]
B=[0,0;1,0;0,0;0,1/r]
C=zeros(2,4)
D=zeros(2)
%C = [0 1 0 0; 0 0 0 1];
% AB=[A,B;C,D]
% FG=expm(AB*10)
% C1=[1,0,0,0;0,0,1,0]

sys=ss(A,B,C,D)
t = linspace(1,Tof,1000)';
u = zeros(length(t),2);
u(1,2) = dv1;

lsim(ss,u,t)