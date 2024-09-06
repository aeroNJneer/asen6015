init

%% prob 1 lunar orbit

a1 = 7045;
e1 = 0.23;
TA1 = -142;
% find semi latus rectum to find r

p1 = PfromAe(a1,e1)
r1 = RfromP(p1,e1,TA1)
R1 = [r1;0;0];
% find specific ang momentum
h1 = HfromTa(Gmoon,r1,e1,TA1)
v1 = VfromRA(Gmoon, r1, TA1)
v1rad = Vradial(Gmoon,e1,h1,360+TA1)
v1tran = Vtrans(Gmoon,e1,h1,360+ TA1)
V1 = [v1rad;v1tran;0]
DeltaV1 = [0.3; -0.1;0];
deltav1 = norm(DeltaV1)
V2 = DeltaV1 + V1

v2 =  norm(V2)
a2 = AfromVR(Gmoon,v2,r1)
H2 = cross(R1,V2)
h2 = norm(H2)
E2 = EccfromHRV(Gmoon,R1,V2,H2)
e2 = norm(E2)

TA2 = TAfromRe(E2,R1) * ta_sign(R1,V2)
deltaW = 360+TA1 - TA2

% part f
Msc = 1224;
Mfuel = 156;
Minit = Msc;
Is = 212;

Mf = Mfinal(Minit,deltav1,Is)
deltaM = Minit - Mf

%% problem 2

init

r1 = 6500;
E1 = pi/2;
rp1 = 5915;
rp2 = 5712;
ra2 = 7888;
a1 = r1;
e = EccfromArp(a1,rp1)
v1 = VfromRA(Gmars,r1,a1)
p = PfromAe(a1,e)
TA1 = TAfromEA(e,E1)

h1 = HfromP(Gmars,p)
vrad1 = Vradial(Gmars,e,h1,TA1)
Vtran1 = Vtrans(Gmars,e,h1,TA1)
V1 = [vrad1;Vtran1;0]
v1check = norm(V1)

a2 = (rp2+ra2)/2

e2 = EccfromRpRa(rp2,ra2)
p2 = PfromAe(a2,e2)
h2 = HfromP(Gmars,p2)
TA2 = TAfromHr(Gmars,h2,r1,e2)

vr2 = Vradial(Gmars,e2,h2,TA2)
vt2 = Vtrans(Gmars,e2,h2,TA2)
v2 = VfromRA(Gmars,r1,a2)
V2 = [vr2;vt2;0];
v2check = norm(V2)
DeltaV2 = V2-V1
deltav2 = norm(DeltaV2)


%% problem 3

clear
init;
format long;
% Hohmann transfers

rpt = AofEarthHelio * AU;
rat = AofSaturnHelio * AU;

a_t = 0.5 * (AofSaturnHelio + AofEarthHelio)*AU
v1i = v1_init(Gsun,rpt)

v2f = v2_final(Gsun,rat)
e = EccfromRpRa(rpt,rat)

Et = EfromA(Gsun,a_t)
vatransfer = VfromRA(Gsun,rat,a_t)
vptransfer = VfromRA(Gsun,rpt,a_t)

deltav1 = vptransfer - v1i
deltav2 = v2f - vatransfer
deltavtot = deltav1+deltav2
TOF = 0.5 * PeriodfromA(Gsun,a_t)
TOFday = TOF/Day
TOFyear = TOF/Year

mmsc = MMfromMua(Gsun,a_t)
mmsaturn = MMfromMua(Gsun,rat)
synod1 = synodFromMm(mmsc,mmsaturn)
synodcheck = synodFromMua(Gsun,a_t,rat)

synodyear = synod1/Year

phase = pi -( mmsaturn * TOF)
phasedeg = rad2deg(phase)

% bi elliptic transfer
rb = 11*AU;
r1i = rpt;
v1i = v1_init(Gsun,r1i)
r3f = rat;
v3f = v2_final(Gsun,rat)
 
at1 = 0.5*(rpt + rb)
et1 = (rb - r1i)/(rb + r1i)


at2 = 0.5*(rb+rat)
et2 = EccfromRpRa(rat,rb)

et2check = EccfromRpRa(9.555,11)

vpt1 = VfromRA(Gsun,rpt,at1)
vat1 = VfromRA(Gsun,rb,at1)

vpt2 = VfromRA(Gsun,rat,at2)
vat2 = VfromRA(Gsun,rb,at2)
deltav1t1 = vpt1 - v1i
deltav2t1 = vat2 - vat1
deltav3t2 = v3f - vpt2

deltavb = abs(deltav1t1) + abs(deltav2t1) + abs(deltav3t2)

TOFbi = pi * (sqrt(at1^3/Gsun) + sqrt(at2^3)/Gsun)
TOFbiYear = TOFbi/Year
deltaV = deltavtot - deltavb
deltaTOF = TOF - TOFbi
deltaTOF / Year
deltaTOF / Day

% part 
deltavtot_ = deltavtot*1000
mratio = Mfinal(1,deltavtot,212)

%% problem 4

clear 
init
a = 6600;
e = 0.46;
i = 74.2;
raan = 207;
w = 10;
TA = 180;
mu = 42838.314;

p = PfromAe(a,e)
h = HfromP(mu,p)
E = EfromA(mu,a )
EA = EAfromTatan(e,TA)
TOF = 0.08639;
TOFsec = TOF*Day
TOFmin = TOFsec / 60
TOFhour = TOFsec / Hour

t2 = 27113.58638948632;1
t3 = 27113.60193770004 ;
dt = t3-t2
dtsec = dt*Day
dtmin = dtsec / Minute
dthour = dtsec / Hour
period = PeriodfromA(mu,a)
halfperiodhour = period / Hour /2
mmcheck = MMfromP(period)
mm = MMfromMua(mu,a )

EA2 = 333.1636;
EA3 = 22.55;
ea2 = EA2-360
EA2rad = deg2rad(EA2)-(2*pi)
EA3rad = deg2rad(EA3)
T500 = T_bfroma(mm,e,EA2rad,EA3rad)
T500i = T_pastPeri(mm,e,EA2rad)
T500f = T_pastPeri(mm,e,EA3rad)
Tcheck = abs(T500f)+abs(T500i)
Tcheck / Minute
T500min = T500/Minute