function [E] = FindEccAnomFromMeanAnom(M,e)
    E = M;
    k = 1;
    err = 1e-10;
    % stabilito intervallo di ricerca (sx,dx) e max valore derivata dymax;
    while (k>err)
        y = e*sin(E)+M;
        k = abs(abs(E)-abs(y));
        E = y;
    end
end