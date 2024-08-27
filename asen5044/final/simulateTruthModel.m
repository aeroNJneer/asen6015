function [X, Y] = simulateTruthModel(t, x0, Q, R, stationNums)
    % This function performs a complete nonlinear simulation of satellite
    % dynamics and tracking station measurements with process and
    % measurement noise. There is no forcing input.
    %
    % Inputs:
    % t = Matrix of time steps that you want to simulate at
    % x0 = 4x1 initial state matrix
    % Q = process noise covariance matrix
    % R = measurement noise covariance matrix
    % stationNums = matrix containing the IDs of the stations that are
    % generating measurements
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
    %
    % Outputs:
    % X = time-tagged truth satellite state vector corresponding to input t
    % Y = cell array of containing time-tagged maps of truth nonlinear
    % measurements from tracking stations corresponding to input t
    %
    % truthY is a cell array containing maps for each tracking station. For
    % example, truthY[1] is the tracking data for station 1, and truthY[10]
    % is the tracking data for station 10. Each element is a map whose key
    % is the timestep and value is the 3x1 measurement matrix.

    % Constant
    mu = 398600;
    
    % Define nonlinear model as anonymous function and set options
    xdot = @(x, u, w) [x(2);
        -mu*x(1)/sqrt((x(1)^2+x(3)^2)^3) + u(1) + w(1);
        x(4);
        -mu*x(3)/sqrt((x(1)^2+x(3)^2)^3) + u(2) + w(2)];
    options = odeset('RelTol', 1e-12);

    % Initialize truth state vector and Cholesky decomp of Q matrix
    if Q == zeros(2)
        % No noise case
        SvQ = Q;
    else
        % Noisy case
        SvQ = chol(Q, 'lower');
    end

    % Initialize truth state vector and tracking data maps
    X = x0;
    Y = {};
    for i = 1:length(stationNums)
        Y = [Y; {containers.Map('KeyType', 'double','ValueType', 'any')}];
    end

    % Perform Cholesky decomp on R matrix
    if R == zeros(3)
        % No noise case
        SvR = R;
    else
        % Noisy case
        SvR = chol(R, 'lower');
    end
    
    % Integrate orbit trajectory stepwise with process noise
    for i = 2:length(t)
        % Draw random noise
        w = SvQ * randn(size(SvQ, 1), 1);
        
        % Define wrapper function for use in ODE45
        xdotwrap = @(t, x) xdot(x, [0; 0], w);

        % Compute next step using ODE45
        [~, xout] = ode45(xdotwrap, [t(i-1) t(i)], X(:, end), options);

        % Store state output
        X = [X xout(end, :)'];
        
        % Compute measurements
        for j = 1:length(stationNums)
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