function E = stateVector2orbitalElements(mu, r0, v0)

% From cartesian to osculating Keplerian
%%
% Transforms cartesian coordinates (r,v) to Keplerian elements (a,e,i,W,w,E).
% Note that we use the eccentric anomaly (or Gudermannian if e > 1)
%
% @param[in] r0 cartesian position vector.
% @param[in] v0 cartesian velocity vector.
% @param[in] mu gravitational parameter.
%
% @param[out] E the osculating Keplerian parameters (a,e,i,W,w,E).
%
% Note: this method will be singular if the elements are not defined or around parabolic motion
%

    % 1 - compute h: the orbital angular momentum vector
    h= cross(r0, v0);

    % 2 - compute p: the orbital parameter
    p = dot(h, h) / mu; 

    % 3 - compute n: the vector of the node line
    % This operation is singular when inclination is zero, in which case the orbital parameters
    % are not defined
    n = cross([0;0;1], h);

    % 4 - compute evett: the eccentricity vector
    R0 = norm(r0);
    Ecc = cross(v0, h)/mu - r0/R0;

    % The eccentricity is calculated and stored as the second orbital element
    E(2)= norm(Ecc);

    % The semi-major axis (positive quantity) is calculated and stored as the first orbital element
    E(1) = abs(p / (1 - E(1) * E(1)));

    % Inclination is calculated and stored as the third orbital element
    E(3) = acos((h(3)) / norm(h));

    % Argument of periapsis is calculated and stored as the fifth orbital element
    E(5) = acos(dot(n,Ecc)/E(1));
    if (Ecc(3) < 0) 
        E(5) = 2 * pi - E(5);
    end

    % Argument of longitude is calculated and stored as the fourth orbital element
    E(4) = acos(n(1)/norm(n));
    if (n(2) < 0) 
        E(4) = 2 * pi - E(4);
    end

    % 4 - We compute ni: the true anomaly (in 0, 2*PI)
    ni = acos(dot(Ecc, r0) / E(2) / R0);    % Eccentric Anomaly

    if (dot(r0,v0) < 0.0) 
        ni = 2 * pi - ni;
    end

    % Eccentric anomaly or the gudermannian is calculated and stored as the sixth orbital element
    if (E(2) < 1.0)
        E(5) = 2.0 * atan(sqrt((1 - E(2)) / (1 + E(2))) * tan(ni / 2.0));  % algebraic kepler's equation
    else
        E(5) = 2.0 * atan(sqrt((E(2) - 1) / (E(2) + 1)) * tan(ni / 2.0));    % algebraic equivalent of kepler's equation in terms of the Gudermannian
    end
end

