function EccAnom = FindEccAnomFromMeanAnom(M,ecc)
    if (M > pi)
        E(1) = M-ecc;
    else
        E(1) = M+ecc;
    end
    
    % E = M as initial guess to see if it converges faster
    E(1) = M;
    
    for i = 1:100
        % Kepler's equation
        Fe = M - E(i) + ecc*sin(E(i));
        % derivative of KE
        dFe = -1 + ecc*cos(E(i));
        % use Newton's method En+1 = En - Fe/dFe
        E(i+1) = E(i) - Fe/dFe;
        % check if value is converging. if so break out
        if abs(E(i+1) - E(i)) < 1e-6
           break;
        end
        if i==100
            disp("could not find value for eccentric anomaly")
        end
    end
%disp("Number of iterations: " + i)
    EccAnom = E(i);
end




