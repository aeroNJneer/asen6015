function Atilde = getAtilde(xNom)
    % Constant
    mu = 398600;
    
    Atilde = [0 1 0 0;...
        3*mu*xNom(1)^2/(xNom(1)^2+xNom(3)^2)^2.5 - mu/(xNom(1)^2+xNom(3)^2)^1.5 0 3*mu*xNom(1)*xNom(3)/(xNom(1)^2+xNom(3)^2)^2.5 0;...
        0 0 0 1;
        3*mu*xNom(1)*xNom(3)/(xNom(1)^2+xNom(3)^2)^2.5 0 3*mu*xNom(3)^2/(xNom(1)^2+xNom(3)^2)^2.5 - mu/(xNom(1)^2+xNom(3)^2)^1.5 0];
end