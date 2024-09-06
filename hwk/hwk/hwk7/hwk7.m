% asen 5050 hwk #7
init

% part a - elliptical hohmann transfer from Earth to Jupier
rperi = AofEarthHelio * AU
rapo = AofJupiterHelio * AU
rclosest = 1e6      % r periapsis of flyby to Jupiter
a = 0.5 * (rapo + rperi)
e = (AofJupiterHelio - AofEarthHelio)/(AofJupiterHelio + AofEarthHelio) 
TOFa = 0.5 * PeriodfromA(Gsun,a)
TOFayears = TOFa / Year
Vt_peri = VfromRA(Gsun,rperi,a) % after first maneuver
Vt_apo = VfromRA(Gsun,rapo,a) % before second maneuver

vearth = sqrt(Gsun/rperi)
Ein = EfromVr(Gsun,Vt_apo,rapo)
hin = sqrt(Gsun*rapo*(1-e))
deltavE = Vt_peri - vearth

% part b hyperbolic orbit of flyby
Vin = Vt_apo
vjupiter = sqrt(Gsun/rapo)
vinfj = Vin - vjupiter
ahypj = -Gjupiter/(vinfj^2)
ecchj = ecchfromArp(ahypj,rclosest)
turnangle = -delta(ecchj)

cosbeta = (vinfj^2 + vjupiter^2 - Vin^2)/(2*vinfj*vjupiter)
betain = acosd(cosbeta)
betaout = betain + turnangle

% v_out calculations
voutj = sqrt( vinfj^2 + vjupiter^2 - (2*vinfj*vjupiter*cosd(betaout)))
[a1 a2 a3] = triangle(-vinfj,vjupiter,voutj)

fpaj = a1
vthetaj = voutj * cosd(fpaj)
vrj = voutj * sind(fpaj)
Voutj = [vrj; vthetaj; 0]

% part d 
R = [rapo; 0; 0];
H = cross(R,Voutj)
Es = EfromVr(Gsun,voutj,rapo)  % was energy increased?


%% prob 2 flyby titan

clear;
init
mtitan = 1.3455e23;
rotitan = 1.221830e6;
asaturn = AofSaturnHelio * AU
rp_sc = 6e5;
ra_sc = 1.8e6;
a_sc = (rp_sc + ra_sc)/2
e_sc = EccfromRpRa(rp_sc,ra_sc)

% find true anomaly
rtside = a_sc*(1-e_sc^2)
TAt = -acosd(((rtside/rotitan) -1)*2)

h = HfromTa(Gsaturn,rotitan,e,TAt)
vr = Vradial(Gsaturn,e,h,TAt)
vt = Vtrans(Gsaturn,e,h,TAt)

% v-inf titan
Vin = [vr;vt;0]
vin = norm(Vin)
Esin = EfromVr(Gsaturn,vin,rotitan)
vsaturn = [0;sqrt(Gsun/asaturn); 0]
Vtitan = [0; sqrt(Gsaturn/rotitan); 0]
vtitan = norm(Vtitan)
Vinftitan = Vin - Vtitan
vinftitan = norm(Vinftitan)

% hyperbolic orbit wrt titan
Gtitan = G*mtitan
ahtitan = -Gtitan/(vinftitan^2)
ehtitan = ecchfromArp(ahtitan,3000)
deltatitan = delta(ehtitan)

% find v out
bin = acosd( (vinftitan^2 + vtitan^2 - vin^2)/(2*vinftitan*vtitan))
bout = bin - deltatitan
vouttitan = vout(vinftitan,vtitan,bout)

[a1 a2 a3] = triangle(vinftitan,vtitan,vouttitan)
fpa = -a1
vtheta = vouttitan * cosd(fpa)
vr = vouttitan * sind(fpa)
Vout = [vr; vtheta; 0]
vout = norm(Vout)

% part d resultant elliptical orbit wrt saturn after flyby
R = [rotitan; 0; 0];
H = cross(R,Vout)
Es = EfromVr(Gsaturn,vout,rotitan)
aout = AfromE(Gsaturn,Es)
eout = EccfromHE(Gsaturn,norm(H),Es)
TAout = TAfromHr(Gsaturn,norm(H),rotitan,eout)
rpout = Rp(aout,eout)
raout = Ra(aout,eout)

% part e delta v equivalent
deltaVeq = Vout - Vin
deltav = norm(deltaVeq)
deltaveq = deltavequiv(vinftitan,deltatitan)

%% prob 3 GMAT problem
clear;
init

perialtJ = 9.829589e6;
iJ = 3.378;
vinfJ = 10.272;
rpJ = perialtJ + Rjupiter
aJ = -Gjupiter/(vinfJ^2)
eJ = 1-(rpJ/aJ)
deltaJ = delta(eJ)

perialtS = 4.036874e6;
rpS = perialtS + Rsaturn
vinfS = 5.3572;
aS = -Gsaturn/(vinfS^2)

eS = 1-(rpS/aS)
deltaS = delta(eS)

rpSs = 20732.114 + Rsaturn



Vx = 5.289894007442628;
Vy = 2.209289843371483E+01;
Vz = -2.627701845286099E-01;

   CassiniVX = 5.28832476237
   CassiniVY = 22.0972626796
   CassiniVZ = -0.262665585359


   deltavx = CassiniVX - Vx
   deltavy = CassiniVY - Vy
   deltavz = CassiniVZ - Vz

   vpa(deltavx)
   vpa(deltavy)
   vpa(deltavz)