function generatePlots(type, t, title, data, varargin)
    if (type == "msmt")
        plotMeasurements(t, data, title)
    elseif (type == "states")
        plotStates(t, data, title)
    elseif (type == "errors")
        data1 = varargin{1,:};
        plotEstimationErrors(t, title, data, data1 )
    elseif (type == "2sigma")
        plot2Sigma(t, title, data)
    end
end

function plotStates(t, x, title)
    figure;
    sgtitle(title);
    subplot(4, 1, 1);
        plot(t, x(1, :));
        ylabel('X (km)');
        grid on;
    subplot(4, 1, 2);
        plot(t, x(2, :));
        ylabel('Xdot (km/s)');
        grid on;
    subplot(4, 1, 3);
        plot(t, x(3, :));
        ylabel('Y (km)');
        grid on;
    subplot(4, 1, 4);
        plot(t, x(4, :));
        ylabel('Ydot (km/s)');
        xlabel('Time t (s)');
    grid on;
end
 
function plotMeasurements(t, trackingData, title)
% Plot measurements
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

        % Plot with with single color for this tracking station
        subplot(3,1,1);
            hold on
            scatter(time, range,5);
        subplot(3,1,2);
            hold on
            scatter(time, rangeRate,5);
        subplot(3,1,3);
            hold on
            scatter(time, elevation,5);
    end

    % Add labels
    sgtitle(title);
    subplot(3,1,1);
        ylabel('\rho Range (km)');
    subplot(3,1,2);
        ylabel('\rhodot Range Rate (km/s)');
    subplot(3,1,3);
        ylabel('\phi Elevation (rad)');
        xlabel('Time t (s)');
    hold off;
end

function plotEstimationErrors(t, title, estErr, estErrCov)
    figure;
    threeSigma1 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma1 = [threeSigma1 2*sqrt(estErrCov(1,i*size(estErrCov, 1)-3))];
    end
    subplot(4,1,1);
    hold on;
    plot(t(2:end), estErr(1,:), 'b');
    plot(t(2:end), threeSigma1, '--b');
    plot(t(2:end), -threeSigma1, '--b');
    sgtitle(title);
    ylabel('X Est. Error [km]');
    legend('Est. Error', '+/- 2\sigma bounds')
    hold off;

    threeSigma2 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma2 = [threeSigma2 2*sqrt(estErrCov(2,i*size(estErrCov, 1)-2))];
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
        threeSigma3 = [threeSigma3 2*sqrt(estErrCov(3,i*size(estErrCov, 1)-1))];
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
        threeSigma4 = [threeSigma4 2*sqrt(estErrCov(4,i*size(estErrCov, 1)))];
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


function plot2Sigma(t, title, estErrCov)
    figure;
    threeSigma1 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma1 = [threeSigma1 2*sqrt(estErrCov(1,i*size(estErrCov, 1)-3))];
    end
    subplot(4,1,1);
    plot(t(2:end), threeSigma1);
    sgtitle(title);
    ylabel('X Est. Error [km]');
    ylim([0 0.5])

    threeSigma2 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma2 = [threeSigma2 2*sqrt(estErrCov(2,i*size(estErrCov, 1)-2))];
    end
    subplot(4,1,2);
    hold on;
    plot(t(2:end), threeSigma2);
    ylabel('Xdot Est. Error [km/s]');
    ylim([0.0005 0.002])

    threeSigma3 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma3 = [threeSigma3 2*sqrt(estErrCov(3,i*size(estErrCov, 1)-1))];
    end
    subplot(4,1,3);
    plot(t(2:end), threeSigma3);
    ylabel('Y Est. Error [km]');
    ylim([0 0.5])

    threeSigma4 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma4 = [threeSigma4 2*sqrt(estErrCov(4,i*size(estErrCov, 1)))];
    end
    subplot(4,1,4);
    plot(t(2:end), threeSigma4);
    ylabel('Ydot Est. Error [km/s]');
    xlabel('Time t (s)');
    ylim([0.0005 0.002])

end

