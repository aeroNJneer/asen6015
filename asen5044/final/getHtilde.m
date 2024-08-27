function Htilde = getHtilde(xNom, stationState)
    % Compute relative range
    rho = sqrt((xNom(1)-stationState(1))^2 + (xNom(3)-stationState(3))^2);

    % Compute elements of Htilde
    H11 = (xNom(1)-stationState(1))/((xNom(1)-stationState(1))^2 + (xNom(3)-stationState(3))^2)^0.5;
    H12 = 0;
    H13 = (xNom(3)-stationState(3))/((xNom(1)-stationState(1))^2 + (xNom(3)-stationState(3))^2)^0.5;
    H14 = 0;
    H21 = (xNom(2)-stationState(2))/rho - (xNom(1)-stationState(1))*((xNom(1)-stationState(1))*(xNom(2)-stationState(2))+(xNom(3)-stationState(3))*xNom(4)-stationState(4))/rho^3;
    H22 = (xNom(1)-stationState(1))/rho;
    H23 = (xNom(4)-stationState(4))/rho - (xNom(3)-stationState(3))*((xNom(1)-stationState(1))*(xNom(2)-stationState(2))+(xNom(3)-stationState(3))*xNom(4)-stationState(4))/rho^3;
    H24 = (xNom(3)-stationState(3))/rho;
    H31 = (stationState(3)-xNom(3))/rho^2;
    H32 = 0;
    H33 = (xNom(1)-stationState(1))/rho^2;
    H34 = 0;

    % Output
    Htilde = [H11 H12 H13 H14;...
        H21 H22 H23 H24;...
        H31 H32 H33 H34];
end