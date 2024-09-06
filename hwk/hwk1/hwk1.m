clear
init;


% problem 3

r1 = 200 + Rearth;
r2 = 35786;

e = EccfromRpRa(r1,r2)
a = (r1+r2)/2
P = PeriodfromA(Gearth, a)/2
dv1 = sqrt((Gearth/a)*(1+e)/(1-e)) -sqrt(Gearth/r1)
dv2 = sqrt(Gearth/r2) - sqrt((Gearth/a)*(1-e)/(1+e))

% simulate hohmann transfer trajectory

E0 = 0;
time = [0:60:P];
for i = 1:length(time)
      EA(i) = FindEccAnomaly(i, a, e, Gearth);
      R(i) = RfromEA(a,e,EA(i));
end