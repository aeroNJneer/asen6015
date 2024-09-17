function deltaCOE = propogate_orbital_elements(COE, mu, perturbations, dt)
    % COE: Array of iniital classical orbital elements
    % mu: Gravitational parameter (GM) in m^3/s^2
    % perturbations: Perturbative accelerations [fR; fT; fN] in m/s^2
    % dt: Time step in seconds
    init    % load function signatures

    a = COE(1);
    e = COE(2);
    i = COE(3);
    Raan = COE(4);
    omega = COE(5);
    M = COE(6);
    
    % Calculate semi latus rectum
    p = PfromAe(a,e);

    % Calculate h
    h = HfromP(mu,p);

    % Calculate mean motion
    n = MMfromMua(mu,a);

    % Calculate the eccentric and true anomaly
    E = M;
    dE = 0;
    while (dE-E) > 1e-6
        E = M - E + e*sin(E);
        dFe = -1 +e*cos(E);
        dE = E - E/dFe;
    end

    TA = TAfromEA(e, E);

    % Calculate r
    r = RfromP(p,e,TA);

    % Perturbative accelerations
    fR = perturbations(1);
    fT = perturbations(2);
    fH = perturbations(3);

    da_dt = (2*a^2/h) * (e*sin(TA)*fR + p*fT/r);
    de_dt = (p*sin(TA)*fR + (p+r)*cos(TA) + (r*e*fT))/h;
    di_dt = r * cos(TA+omega) * fH/h;
    dRaan_dt = r* sin(TA+ omega) *fH/(h* sin(i));
    domega_dt = (-p*cos(TA)*fR + (p+r)*sin(TA)*fT)/(h*e) - r*sin(TA+2)*cos(i)*fH/(h*sin(i));
    dM_dt = n + sqrt(1-e^2) * ((p*cos(TA) - 2*r*e)*fR - (p+r)*sin(TA)*fT)/(h*e);

    % Integrate over time step dt
    a = a + da_dt * dt;
    e = e + de_dt * dt;
    i = i + di_dt * dt;
    Raan = Raan + dRaan_dt * dt;
    omega = omega + domega_dt * dt;
    M = M + dM_dt * dt;
    deltaCOE = [a, e, i, Raan, omega, M]';
end
