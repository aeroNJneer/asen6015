function [positions,velocities]  =  simulate_orbit(mu,r_eci,v_eci,n)
%SIMULATE_ORBIT simulates a satellite trajectory in cartesian coordinate by solving the kepler's
%equation. 
% INPUT:
%   r_eci = > initial position vector in ECI cooridnate system 
%   v_eci = > initial velocity vecotr in ECI coordinate system 
%   n     = > n is the number of periods over which the simulation is carried out 
%
% OUTPUT: 
%   Positon     = > Position[t] is the ECI Position of the satellitecomputed at time t  
%   Velocities  = > velocities[t] is the ECI velocity of the satellitecomputed at time t  

% Initial values, parameters and variable definitions
init 
H = cross(r_eci,v_eci);
p = norm(H)^2/mu;
[a,e,i,Raan,argument_of_periapsis,true_anomaly] = ECI2classical(r_eci,v_eci);    %ECI to classical orbital elements convertion Problem I-1
T = 2*pi()*a^1.5/mu^0.5;   %Orbital period  Problem I-2
runtime =  n*floor(T);       %for how long to simulate the orbit propagation. Five orbital period length.
t_a = [];                  % t_a[t+1] is the true anomaly at time t
pi_dec = 0;

c_w = [cos(argument_of_periapsis) sin(argument_of_periapsis) 0
    -sin(argument_of_periapsis) cos(argument_of_periapsis) 0
    0 0 1];
c_i = [1 0 0
    0 cos(i) sin(i)
    0 -sin(i) cos(i)];
c_omega = [cos(Raan) sin(Raan) 0
    -sin(Raan) cos(Raan) 0
    0 0 1];
C = c_w*c_i*c_omega;  


% simulates a satellite trajectory in cartesian coordinate by solving the kepler's
% equation. 

% Handles possible issues with quadrant ambiguity
if true_anomaly >= pi
        true_anomaly =  true_anomaly-pi;
        pi_dec = 1;
end

t = 0;
t_a(t+1) = true_anomaly;

positions = [];
velocities = [];
while t<runtime
    if pi_dec
        true_anomaly = true_anomaly+pi;
        E = 2*atan(tan(true_anomaly/2)*sqrt((1-e)/(1+e)));
        r = a*(1-e^2)/(1+e*cos(t_a(t+1)));
        x = a*(cos(E)-e);
        y = a*sqrt(1-e^2)*sin(E);
        xdot = -sqrt(mu*a)*sin(E)/r;
        ydot = sqrt(mu*a*(1-e^2))*cos(E)/r;
        position = C'*[x;y;0];
        velocity = C'*[xdot;ydot;0];
        positions(t+1,1:3) = position';
        velocities(t+1,1:3) = velocity';
        true_anomaly = true_anomaly-1;   
    else 
        E = 2*atan(tan(true_anomaly/2)*sqrt((1-e)/(1+e)));
        r = a*(1-e^2)/(1+e*cos(t_a(t+1)));
        x = a*(cos(E)-e);
        y = a*sqrt(1-e^2)*sin(E);
        xdot = -sqrt(mu*a)*sin(E)/r;
        ydot = sqrt(mu*a*(1-e^2))*cos(E)/r;
        position = C'*[x;y;0];
        velocity = C'*[xdot;ydot;0];
        positions(t+1,1:3) = position';
        velocities(t+1,1:3) = velocity';

    end

    M =  E-e*sin(E);
    delta_t = M*a^1.5/sqrt(mu);
    Tp = t-delta_t;
    t = t+1;
    delta_t = t-Tp;
    E = FindEccAnomaly(delta_t,a,e,mu)
    true_anomaly = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
    t_a(t+1) = true_anomaly;
end

% applies correction if prior corrections were employed to avoid quadrant
% ambiguity 
if pi_dec
         t_a = t_a+ones(1,length(t_a))*pi();
end


end