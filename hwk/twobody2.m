function [rf, vf] = twobody2(tau, ri, vi)
  
% solution of the two body initial value problem
  
% Shepperd's method
% input
%  tau = propagation time interval (seconds)
%  ri  = initial eci position vector (kilometers)
%  vi  = initial eci velocity vector (km/sec)
% output
%  rf = final eci position vector (kilometers)
%  vf = final eci velocity vector (km/sec)
% Orbital Mechanics with MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global mu
% convergence criterion
tol = 1.0e-10;
rdotv = dot(ri, vi);
   
rm1 = sqrt(dot(ri, ri));
   
beta = (2 * mu / rm1) - dot(vi, vi);
u = 0;
    
if (beta ~= 0)
    
   umax = 1 / sqrt(abs(beta));
   
else
    
   umax = 1d+24;
   
end
umin = -umax;
if (beta <= 0)
    
   delu = 0;
   
else
    
   p = 2 * pi * mu * beta ^ (-1.5);
   n = fix((1 / p) * (tau + 0.5 * p - 2 * rdotv / beta));
   delu = 2 * n * pi * beta ^ (-2.5);
   
end
tsav = 1.0e99;
niter = 0;
% kepler iteration loop
while (1)
    
  niter = niter + 1;
      
  q = beta * u * u;
  q = q / (1 + q);
  u0 = 1 - 2 * q;
  u1 = 2 * u * (1 - q);
      
  % continued fraction iteration
  n = 0;
  l = 3;
  d = 15;
  k = -9;
  a = 1;
  b = 1;
  g = 1;
  while (1)
    gsav = g;
    k = -k;
    l = l + 2;
    d = d + 4 * l;
    n = n + (1 + k) * l;
    a = d / (d - n * a * q);
    b = (a - 1) * b;
    g = g + b;
    
    if (abs(g - gsav) < tol)
       break;
    end   
  end  
  
  uu = (16 / 15) * u1 * u1 * u1 * u1 * u1 * g + delu;
       
  u2 = 2 * u1 * u1;
  u1 = 2 * u0 * u1;
  u0 = 2 * u0 * u0 - 1;
  u3 = beta * uu + u1 * u2 / 3;
      
  r = rm1 * u0 + rdotv * u1 + mu * u2;
  t = rm1 * u1 + rdotv * u2 + mu * u3;
       
  dtdu = 4 * r * (1 - q);
     
  % check for time convergence
  if (abs(t - tsav) < tol)
     break;
  end   
  usav = u;
  tsav = t;
  terr = tau - t;
  if (abs(terr) < abs(tau) * tol)
     break;
  end   
  du = terr / dtdu;
      
  if (du < 0)
     umax = u;
     u = u + du;
     if (u < umin)
        u = 0.5 * (umin + umax);
     end   
  else
     umin = u;
     u = u + du;
     if (u > umax)
        u = 0.5 * (umin + umax);
     end   
  end
  % check for independent variable convergence
  if (abs(u - usav) < tol)
     break;
  end   
  
  % check for more than 20 iterations
  
  if (niter > 20)
     return;
  end   
end
uc = u;
fm = -mu * u2 / rm1;
ggm = -mu * u2 / r;
f = 1 + fm;
g = rm1 * u1 + rdotv * u2;
ff = -mu * u1 / (rm1 * r);
gg = 1 + ggm;
% compute final state vector
for i = 1:1:3
    rf(i) = f * ri(i) + g * vi(i);
    vf(i) = ff * ri(i) + gg * vi(i);
end
