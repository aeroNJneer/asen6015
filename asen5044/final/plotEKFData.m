function plotEKFData(t, x, estErr, estErrCov, trackingData) 
    %% Plot nonlinear truth dynamics model simulation

    % Total states
    figure;
    subplot(4, 1, 1);
    plot(t, x(1, :));
    title('Noisy Truth States vs. Time, Typical EKF Run');
    ylabel('X (km)'); ylim([-10000, 10000])
    grid on;
    subplot(4, 1, 2);
    plot(t, x(2, :));
    ylabel('Xdot (km/s)'); ylim([-10, 10])
    grid on;
    subplot(4, 1, 3);
    plot(t, x(3, :));
    ylabel('Y (km)'); ylim([-10000, 10000])
    grid on;
    subplot(4, 1, 4);
    plot(t, x(4, :));
    ylabel('Ydot (km/s)'); ylim([-10, 10])
    xlabel('Time t (s)');
    grid on;
    
    %% Plot nonlinear model data simulation
    figure;
    for j = 1:length(trackingData)
        % Instantiate arrays for scatter plot for time and each measurement
        time = [];
        range = [];
        rangeRate = [];
        elevation = [];

        % Get tracking data for current station
        data = trackingData{j};
        for i=1:length(t)
            % Iterate through map to get all the times where a measurement
            % is available
            if isKey(data, t(i))
                trackPoint = data(t(i));
                time = [time; t(i)];
                range = [range; trackPoint(1,1)];
                rangeRate = [rangeRate; trackPoint(2,1)];
                elevation = [elevation; trackPoint(3,1)];
            end
        end

        % Plot with single color for this tracking station
        subplot(3,1,1);
        hold on
        scatter(time, range);
        subplot(3,1,2);
        hold on
        scatter(time, rangeRate);
        subplot(3,1,3);
        hold on
        scatter(time, elevation);
    end

    % Add labels
    subplot(3,1,1);
    sgtitle('Noisy Simulated Measurements, Typical EKF Run');
    ylabel('Range (km)');
    subplot(3,1,2);
    ylabel('Range Rate (km/s)');
    subplot(3,1,3);
    ylabel('Elevation (rad)');
    xlabel('Time t (s)');
    hold off;

    %% Plot EKF estimation errors
    figure;
    threeSigma1 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma1 = [threeSigma1 3*sqrt(estErrCov(1,i*size(estErrCov, 1)-3))];
    end
    subplot(4,1,1);
    hold on;
    plot(t(2:end), estErr(1,:), 'b');
    plot(t(2:end), threeSigma1, '--b');
    plot(t(2:end), -threeSigma1, '--b');
    sgtitle('State Estimation Errors vs. Time, Typical EKF Run');
    ylabel('X Est. Error [km]');xlim([1000,14000])
    legend('Est. Error', '+/- 3\sigma bounds')
    hold off;

    threeSigma2 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma2 = [threeSigma2 3*sqrt(estErrCov(2,i*size(estErrCov, 1)-2))];
    end
    subplot(4,1,2);
    hold on;
    plot(t(2:end), estErr(2,:), 'b');
    plot(t(2:end), threeSigma2, '--b');
    plot(t(2:end), -threeSigma2, '--b');
    ylabel('Xdot Est. Error [km/s]');
    hold off;

    threeSigma3 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma3 = [threeSigma3 3*sqrt(estErrCov(3,i*size(estErrCov, 1)-1))];
    end
    subplot(4,1,3);
    hold on;
    plot(t(2:end), estErr(3,:), 'b');
    plot(t(2:end), threeSigma3, '--b');
    plot(t(2:end), -threeSigma3, '--b');
    ylabel('Y Est. Error [km]');
    hold off;

    threeSigma4 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma4 = [threeSigma4 3*sqrt(estErrCov(4,i*size(estErrCov, 1)))];
    end
    subplot(4,1,4);
    hold on;
    plot(t(2:end), estErr(4,:), 'b');
    plot(t(2:end), threeSigma4, '--b');
    plot(t(2:end), -threeSigma4, '--b');
    ylabel('Ydot Est. Error [km/s]');
    xlabel('Time t (s)');
    hold off;
end