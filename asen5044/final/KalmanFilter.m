classdef KalmanFilter
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        InitialState
        NominalState
        TruthState
        PerturbedState
        TotalState
        ModelState
        EstimationError
        EstimationErrorCovariance
        MeasurementInnovation
        MeasurementInnovationCovariance
        NonlinearModel
    end

    methods
        function obj = set.ModelState(obj, modelState)
            obj.ModelState = modelState;
        end
        
        function [tout, modelState] = modelState(t, nonlinearModel, initialState, mu, w) 
            % Define nonlinear model as anonymous function and set options
            xdot = nonlinearModel;
            options = odeset('RelTol', 1e-12);
        
            % Define wrapper function for use in ODE45
            xdotwrap = @(t, x) xdot(x, [0; 0], w);
            % Compute total state using ODE45
            [tout, modelState] = ode45(xdotwrap,t, initialState, options);
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end