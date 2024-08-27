% HWK #2
 % Problem 1

r = 6678; % km
k = 398600; % km^3/s^2
w = sqrt(k/r^3);
deltat = 10;
x0=[10;-0.5;0;2.5e-5];

A1=[0,1,0,0; 3*(w^2), 0,0,2*r*w; 0,0,0,1; 0,-2*w/r, 0,0];

B1=[0,0;1,0;0,0;0,1/r];

C1=[1,0,0,0; 0,0,1,0];

D=zeros();

sys1=ss(A1,B1,C1,D); 

Ahat=[A1,B1;zeros(2,4),zeros(2,2)];

eAhat=expm(Ahat*deltat);



F=eAhat(1:4,1:4);
G=eAhat(1:4,5:end);

sysDis=ss(F,G,C1,D);

