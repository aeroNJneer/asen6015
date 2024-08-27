function [measY] = parseGivenData(tvec, data)
    % This function takes the given class data and turns it into a map that
    % can be read by our filters
    % Initialize
    stationNums = 1:12;
    measY = {};
    for i = 1:length(stationNums)
        measY = [measY; {containers.Map('KeyType', 'double','ValueType', 'any')}];
    end

    for i = 1:length(data)
        % Extract cell at time
        trackAtTime = data{i};

        % Check to see if data exists for that time step
        if ~isnan(trackAtTime)
            % Iterate through all tracking data at timestep
            for j = 1:size(trackAtTime,2)
                trackPoint = trackAtTime(:,j);

                % Query station map, add tracking data, and return to map array
                stationData = measY{stationNums==trackPoint(4)};
                stationData(tvec(i)) = trackPoint(1:3);
                measY{stationNums==trackPoint(4)} = stationData;
            end
        end
    end
end