function [Y] = getMeasurements(Y, t, X, stationNums, SvR)

% Compute measurements
for j = 1:length(stationNums)
    for i = 1:length(t)
            % Identify the current station's state vector
            stationID = stationNums(j);
            stationState = getStationState(t(i), stationID);
    
            % Compute satellite elevation angle and station angle to check
            % when satellite becomes visible
            satElevation = atan2(X(3, i)-stationState(3), X(1, i)-stationState(1));
            stationTheta = atan2(stationState(3), stationState(1));
    
            % Compute smallest angle between satellite and tracking station
            angleDiff = satElevation - stationTheta;
            angleDiff = abs(mod((angleDiff + pi), 2*pi) - pi);
    
            % Check to see if satellite is visible to station
            if angleDiff <= pi/2
                % Compute measurement with process noise
                outY = getNLMeasurement(t(i), X(:,i), stationID) + SvR*randn(size(SvR, 1), 1);
                
                % Query station map, add tracking data, and return to map array
                stationData = Y{stationNums==stationID};
                stationData(t(i)) = outY;
                Y{stationNums==stationID} = stationData;
            end
        end
    end
end