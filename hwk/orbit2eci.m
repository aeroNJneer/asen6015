function C = orbit2eci(raan,i,theta)
% DCM for orbit (r,theta,h) frame to ECI (x,y,z) frame
%   raan: Right ascension of the ascending node (rad)
%   i: inclination (rad)
%   theta: argument of periapsis (w) + true anomaly (TA)

st1 = sin(raan);
ct1 = cos(raan);
st2 = sin(i);
ct2 = cos(i);
st3 = sin(theta);
ct3 = cos(theta);

C(1,1) = ct3*ct1-st3*ct2*st1;
C(1,2) = ct3*st1+st3*ct2*ct1;
C(1,3) = st3*st2;
C(2,1) = -st3*ct1-ct3*ct2*st1;
C(2,2) = -st3*st1+ct3*ct2*ct1;
C(2,3) = ct3*st2;
C(3,1) = st2*st1;
C(3,2) = -st2*ct1;
C(3,3) = ct2;

end