function [theta, raan, inc] = RSW2Euler(r,v)
init;    
    mu = Gearth;
    H = cross(r,v);
    inc = IfromH(H);
    lon = LONfromH(H);
    raan = RaanfromLON(lon) * raan_sign(lon);
    if isnan(raan)
        raan = 0;
    end
    e = EccfromHRV(mu,r,v,H);
    w = WfromLONe(e,lon) * w_sign(e);
    if isnan(w)
        w = 0;
    end
    TA = TAfromRe(r,e);
    if isnan(TA)
        TA = 0;
    end
    theta = w+TA;
    
end