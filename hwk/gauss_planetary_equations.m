function [a, e, i, RAAN, omega, M] = gauss_planetary_equations(r, v, mu, perturbations, dt)
    % r: Position vector [x; y; z] in meters
    % v: Velocity vector [vx; vy; vz] in meters per second
    % mu: Gravitational parameter (GM) in m^3/s^2
    % perturbations: Perturbative accelerations [fR; fT; fN] in m/s^2
    % dt: Time step in seconds
    init    % load function signatures

    % Calculate theta angle between 

    % Calculate the specific angular momentum
    H = HfromRV(r,v);
    h_norm = norm(H);

    % Calculate the eccentricity vector
    E_vec = EccfromHRV(mu,r,v,H);
    e = norm(E_vec);

    % Calculate the semi-major axis
    a = afromVR(mu,norm(v),norm(r));

    % Calculate semi-latus rectum
    p = pfromAe(a,e);
    % Calculate the inclination
    i = IfromH(H);

    % Calculate the right ascension of the ascending node (RAAN)
    LON = LONfromH(H);
    Raan = RaanfromLON(LON);
    if raan_sign(LON) < 0
        Raan = 2*pi - Raan;
    end

    % Calculate the argument of periapsis
    w = WfromLONe(LON,E);
    if w_sign(E) <0
        w = 2 * pi - w;
    end

    % Calculate the mean and true anomaly
    E = EfromEr(e_vec,r);
    if E_sign(r,v) < 0
        E = 2 * pi - E;
    end

    M = E - e * sin(E);
    TA = TAfromEA(e, E);

    % Perturbative accelerations
    fR = perturbations(1);
    fT = perturbations(2);
    fN = perturbations(3);

    % Update the orbital elements using Gauss planetary equations
    da_dt = 2 * sqrt(a^3 / mu) * (e * sin(E) * fR + (1 + e * cos(E)) * fT);
    de_dt = sqrt(a / mu) * (sin(E) * fR + (cos(E) + cos(E) * e) * fT - sin(E) * e * fN);
    di_dt = sqrt(a / mu) * cos(E + omega) * fN;
    dRAAN_dt = sqrt(a / mu) * sin(E + omega) / sin(i) * fN;
    domega_dt = sqrt(a / mu) * (cos(E) * fR - sin(E) * fT) / e - dRAAN_dt * cos(i);
    dM_dt = sqrt(mu / a^3) + sqrt(a / mu) * (-(2 + e * cos(E)) * fR + (1 - e^2) * fT / e);

    da_dt = (2*a^2/h) * (e*sin(TA)*fR + p*fT/(norm(r)));
    de_dt = (p*sin(TA)*fR + ((P+r)*cos(TA) + (norm(r)*e)*fT))/h;
    di_dt = norm(r) * cos(TA+omega) * fN/h;
    dRAAN_dt = norm(r)* sin(TA+ omega) *fN/(h* sin(i));
    domega_dt = (-p*cos(TA)*fR + (p+r)*sin(TA)*fT)/(h*e) - norm(r)*sin(TA+2)*cos(i)*fH/(h*sin(i));
    dM_dt = sqrt(mu / a^3) + a*sqrt(1-e^2) * ((p*cos(TA) - 2*norm(r)*e)*fR) - (p+r)*sin(TA)*fN;

    % Integrate over time step dt
    a = a + da_dt * dt;
    e = e + de_dt * dt;
    i = i + di_dt * dt;
    RAAN = RAAN + dRAAN_dt * dt;
    omega = omega + domega_dt * dt;
    M = M + dM_dt * dt;
end
