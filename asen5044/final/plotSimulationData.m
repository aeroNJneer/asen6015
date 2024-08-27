function plotSimulationData(t, x, xNoise, y, yNoise, trackingData)
    % Total states
    figure;
    subplot(4, 1, 1);
    plot(t, x(1, :));
    sgtitle(['Monte Carlo Truth Sta' ...
        'tes vs. Time']);
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
    
        figure;
    subplot(4, 1, 1);
    plot(t, xNoise(1, :));
    sgtitle(['Monte Carlo Noisy Truth Sta' ...
        'tes vs. Time']);
    ylabel('X (km)'); ylim([-10000, 10000])
    grid on;
    subplot(4, 1, 2);
    plot(t, xNoise(2, :));
    ylabel('Xdot (km/s)'); ylim([-10, 10])
    grid on;
    subplot(4, 1, 3);
    plot(t, xNoise(3, :));
    ylabel('Y (km)'); ylim([-10000, 10000])
    grid on;
    subplot(4, 1, 4);
    plot(t, xNoise(4, :));
    ylabel('Ydot (km/s)'); ylim([-10, 10])
    xlabel('Time t (s)');
    grid on;
    
    %% Plot nonlinear model data simulation
    figure;
    for j = 1:length(yNoise)
        % Instantiate arrays for scatter plot for time and each measurement
        time = [];
        range = [];
        rangeRate = [];
        elevation = [];

        % Get tracking data for current station
        data = yNoise{j};
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
    sgtitle('Noisy Simulated Measurements, Monte Carlo Run');
    ylabel('Range (km)');
    subplot(3,1,2);
    ylabel('Range Rate (km/s)');
    subplot(3,1,3);
    ylabel('Elevation (rad)');
    xlabel('Time t (s)');
    hold off;

  
end

