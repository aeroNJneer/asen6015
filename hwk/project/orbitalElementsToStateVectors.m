function [r, v] = orbitalElementsToStateVectors(mu, E)

    a = E(0);
    e = E(1);
    i = E(2);
    Raan = E(3);
    omega = E(4);
    EA = E(5);
    
    % semi-major axis is assumed to be positive here we apply the convention of having it negative as for
    % comegautations to result in higher elegance
    if (e > 1) 
        a = -a;
    end
    
    % 1 - We start by evaluating position and velocity in the perifocal reference system
    if (e < 1.0) % EA is the eccentric anomaly
       
        b = a * sqrt(1 - e * e);
        n = sqrt(mu / (a * a * a));
        
        xper = a * (cos(EA) - e);
        yper = b * sin(EA);
        xdotper = -(a * n * sin(EA)) / (1 - e * cos(EA));
        ydotper = (b * n * cos(EA)) / (1 - e * cos(EA));
   else % EA is the Gudermannian
        b = -a * sqrt(e * e - 1);
        n = sqrt(-mu / (a * a * a));
        
        dNdZeta = e * (1 + tan(EA) * tan(EA)) - (0.5 + 0.5 * pow(tan(0.5 * EA + M_PI_4), 2)) / tan(0.5 * EA + M_PI_4);
        
        xper = a / cos(EA) - a * e;
        yper = b * tan(EA);
        
        xdotper = a * tan(EA) / cos(EA) * n / dNdZeta;
        ydotper = b / pow(cos(EA), 2) * n / dNdZeta;
    end
    
    % 2 - We then built the rotation matrix from perifocal reference frame to inertial
    
    cosRaan = cos(Raan);
    cosomega = cos(omega);
    sinRaan = sin(Raan);
    sinomega = sin(omega);
    cosi = cos(i);
    sini = sin(i);
    
    R1(1) = cosRaan * cosomega - sinRaan * sinomega * cosi;
    R1(2) = -cosRaan * sinomega - sinRaan * cosomega * cosi;
    R1(3) = sinRaan * sini;
    R2(1) = sinRaan * cosomega + cosRaan * sinomega * cosi;
    R2(2) = -sinRaan * sinomega + cosRaan * cosomega * cosi;
    R2(3) = -cosRaan * sini;
    R3(1) = sinomega * sini;
    R3(3) = cosomega * sini;
    R3(2) = cosi;
    
    % 3 - We end by transforming according to this rotation matrix
    
    temp(3) = [xper, yper, 0.0];
    temp2(3) = [xdotper, ydotper, 0];
    
    for (j = 1:3) 
        r0(j) = 0.0;
        v0(j) = 0.0;
        for (k = 1:3)
            r0(j) += R(j)(k) * temp(k);
            v0(j) += R(j)(k) * temp2(k);
        end
    end
r = r0;
v = v0;
end
  