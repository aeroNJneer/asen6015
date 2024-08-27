function [y, stationIDs] = getStackedMeasAtTime(t, meas)
    % Initialize empty measurements and station IDs
    y = [];
    stationIDs = [];

    for i=1:length(meas)
        % For each station, check to see if a measurement exists at t
        station = meas{i};
        if isKey(station, t)
            % Measurement was taken, so stack it and note the station ID
            y = [y; station(t)];
            stationIDs = [stationIDs; i];
        end
    end
end