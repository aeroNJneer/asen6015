function [meas] = getNLMeasurement(t, satState, stationID)
    % Compute station state
    stationState = getStationState(t, stationID);

    % Compute nominal NL measurement
    rho = sqrt((satState(1)-stationState(1))^2 + (satState(3)-stationState(3))^2);
    meas = [rho;...
        ((satState(1)-stationState(1))*(satState(2)-stationState(2)) + (satState(3)-stationState(3))*(satState(4)-stationState(4)))/rho;...
        atan2(satState(3)-stationState(3), satState(1)-stationState(1))];
end