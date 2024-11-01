function [a,e,i,Raan,argument_of_periapsis,true_anomaly] = ECI2classical(position,velocity,mu)
%this function converts position and velocity in ECI coordinate to for a
%given eliptical orbit. Position and velocity must be given as row vector. 
%classical orbital elements [a,e,i,omega,argument_of_perigee,true_anomaly]

r = norm(position);   % scalar norm of positon 
v = norm(velocity);  % scalar norm of velocity 

H = cross(position,velocity);    % vector H
h = norm(H);           %scalar norm of H

a = 1/(-v^2/mu + 2/r);
p = h^2/mu;

e_vec  = 1/mu*(cross(velocity,H)-mu*position/r);  % eccentricity vector 
e= norm(e_vec);  %eccentricity 

u_H  = H/h;  %unit vector in H direction 
u_e = e_vec/e; %unit vector in e direction 
u_m = cross(u_H,u_e);

i = acos(u_H(3));
% argument_of_perigee= atan(u_H(1)/u_m(3));

% line of nodes
n = cross([0 0 1],H);
n = n/norm(n);

if dot(position,velocity) > 0
    true_anomaly = acos(dot(u_e,position/norm(position)));
else 
   true_anomaly = 2*pi-acos(dot(u_e,position/norm(position)));
end

argument_of_periapsis = acos(dot(n,u_e));
if e_vec(3)>0
    argument_of_periapsis = acos(dot(n,u_e));
else 
    if acos(dot(n,u_e))>0   
         argument_of_periapsis = 2*pi()-acos(dot(n,u_e));
    end
end

if  n(1)<0 & n(2)<0
    Raan = atan(n(2)/n(1))+pi;
elseif n(1)<0 & n(2)>0;
    Raan = atan(n(2)/n(1))+pi;
elseif n(2)<0 & n(1)>0;
    Raan = atan(n(2)/n(1))+2*pi;
else 
    Raan = atan(n(2)/n(1));
end
