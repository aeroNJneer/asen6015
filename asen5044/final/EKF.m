function [trialTotalStates, trialEstimationErrs, trialEstimationErrsCov, trialMeasInnovations, trialMeasInnovationsCov, xnoisy, measY, trialNEES, trialNIS] = ...
EKF (t, nomX0, pPlus0, Qkf, R, omegaTilde, givenData, dt, stationNums, sim)

    pPlus = pPlus0;    % Selected such that the satellite generally stays outside the Earth
    Sv = chol(pPlus, 'lower');
    
    % Draw intial state with perturbation and conduct truth model sim
    X0 = nomX0 + Sv * randn(size(Sv, 1), 1);

    % Monte Carlo simulation for measurements
    [xnoisy, measY] = simulateTruthModel(t, X0, givenData.Qtrue, R, stationNums);
    if (~sim)
       % Parse data from file
        measY = parseGivenData(givenData.tvec, givenData.ydata);
    end
    xhatPlus = X0;

    % Initialize this run's outputs
    trialTotalStates = {};
    trialEstimationErrs = {};
    trialEstimationErrsCov = {};
    trialMeasInnovations = {};
    trialMeasInnovationsCov = {};
    trialNEES = {};
    trialNIS = {};
    scale=1;
    for j = 2:length(t)           
        Ftilde = eye(4) + getAtilde(xhatPlus)*dt;
        pMinus = Ftilde*pPlus*Ftilde' + omegaTilde*Qkf*omegaTilde'*scale;
        [xhatMinus, ~] = simulateTruthModel([t(j-1) t(j)], xhatPlus, zeros(2), zeros(3), stationNums);
        xhatMinus = xhatMinus(:,end);
        
        % Get stacked station measurements
        [y, stationIDs] = getStackedMeasAtTime(t(j), measY);
        yhat = [];
        HtildeStacked = [];
        RStacked = [];
        if ~isempty(y)
            % If there are measurements, then perform correction step by
            % stacking nominal measurements y* per station and getting corresponding
            % stacked Htilde and R matrices
            for k=1:length(stationIDs)
                % Compute a measurement for station j at time t(i)
                yhat = [yhat; getNLMeasurement(t(j), xhatMinus, stationIDs(k))];
                HtildeStacked = [HtildeStacked; getHtilde(xhatMinus, getStationState(t(j), stationIDs(k)))];
                RSize = size(R, 1);   % R is square
                RStacked = [RStacked [zeros((k-1)*RSize, RSize); R; zeros((length(stationIDs)-k)*RSize, RSize)]];
            end

            % Compute Kalman gain
            K = pMinus * HtildeStacked' / (HtildeStacked*pMinus*HtildeStacked' + RStacked);
    
            % Update estimation error covariance and compute measurement innovation 
            pPlus = (eye(size(K, 1))-K*HtildeStacked)*pMinus;
            dy = y-yhat;
            xhatPlus = xhatMinus + (K * dy);

            % Compute measurement innovation covariance matrix
            S = HtildeStacked*pMinus*HtildeStacked' + RStacked;
        else
            % No measurements available, so carry forward pMinus and state
            pPlus = pMinus;
            xhatPlus = xhatMinus;
            % dy and S are null
            dy = [nan; nan; nan];
            S = [nan nan nan; nan nan nan; nan nan nan];
        end
        % Store outputs
        trialTotalStates = [trialTotalStates xhatPlus];
        trialEstimationErrs = [trialEstimationErrs xnoisy(:, j)-xhatPlus];
        trialEstimationErrsCov = [trialEstimationErrsCov pPlus];
        trialMeasInnovations = [trialMeasInnovations dy];
        trialMeasInnovationsCov = [trialMeasInnovationsCov S];
        trialNEES = [trialNEES (xnoisy(:, j)-xhatPlus)' * (pPlus)^(-1) * (xnoisy(:, j)-xhatPlus)];
        trialNIS = [trialNIS dy'*(S)^(-1)*dy];
     end
end