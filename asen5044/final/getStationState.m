% This function returns the station position [Xi_s(t); Xdoti_s(t); Yi_s(t); Ydoti_s(t)]
% given time t and station number i
function state = getStationState(t, stationNum)
    % Constants
    rE = 6378;
    omegaE = 2*pi/86400;
    
    theta = (stationNum-1) * pi/6;
    state = rE* [cos(omegaE*t+theta);...
        -omegaE*sin(omegaE*t+theta);...
        sin(omegaE*t+theta);...
        omegaE*cos(omegaE*t+theta)];
end