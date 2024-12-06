% Low-Thrust LEO to GEO Transfer Guidance Scheme
% Based on Kluever (2012) paper in Journal of Guidance

classdef LowThrustTransfer4
    properties (Constant)
        % Mission parameters
        MU = 3.986e14;           % Earth's gravitational constant (m^3/s^2)
        G = 9.81;                % Gravitational acceleration (m/s^2)
        
        % Spacecraft parameters
        INITIAL_MASS = 1200;     % kg
        POWER = 10000;           % Watts
        ISPECIFIC = 3300;        % Specific impulse (s)
        EFFICIENCY = 0.65;       % Propulsion efficiency
        
        % Initial and final orbit conditions
        LEO_ALT = 550e3;         % Initial LEO altitude (m)
        LEO_INC = 28.5 * pi/180; % Initial LEO inclination (rad)
        GEO_ALT = 42164e3;       % GEO altitude (m)
    end
    
    methods (Static)
     function R_RTN_to_ECI = computeRotationMatrix(r, v)
            % Compute rotation matrix from RTN to ECI frame
            
            % R (radial) direction
            R_unit = r / norm(r);
            
            % N (normal) direction
            H_vec = cross(r, v);
            N_unit = H_vec / norm(H_vec);
            
            % T (tangential) direction - completes right-handed system
            T_unit = cross(N_unit, R_unit);
            
            % Rotation matrix from RTN to ECI
            R_RTN_to_ECI = [R_unit, T_unit, N_unit]';
        end
        
        function [QU, varargout] = computeThrustDirection(t, state, params)
            % Compute thrust direction based on guidance scheme
            
            % Extract position and velocity
            r = state(1:3);
            v = state(4:6);
            
            % Compute orbital elements
            [a, e, i, ~, ~, mu] = LowThrustTransfer.orbitalElementsFromState(state);
            
            % Flight path angle
            r_mag = norm(r);
            v_mag = norm(v);
            flightPathAngle = asin(r(3) / r_mag);
            
            % Compute orbital plane and local horizontal frame
            R_RTN_to_ECI = LowThrustTransfer.computeRotationMatrix(r, v);
            
            % Compute optimal in-plane and out-of-plane steering angles
            [optimalInPlaneSteer, optimalOutOfPlaneSteering] = ...
                LowThrustTransfer.computeOptimalSteering(a, e, i, mu, v_mag, flightPathAngle);
            
            % Blend guidance parameters
            [inPlaneSteering, outOfPlaneSteering] = ...
                LowThrustTransfer.blendGuidanceControls(t, optimalInPlaneSteer, optimalOutOfPlaneSteering, params);
            
            % Compute thrust direction unit vector in RTN frame
            QU_RTN = [sin(inPlaneSteering) * cos(outOfPlaneSteering);
                      cos(inPlaneSteering) * cos(outOfPlaneSteering);
                      sin(outOfPlaneSteering)];
            
            % Transform thrust direction to ECI frame
            QU = R_RTN_to_ECI * QU_RTN;
            
            % Optional outputs for debugging
            if nargout > 1
                varargout{1} = inPlaneSteering;
                varargout{2} = outOfPlaneSteering;
                varargout{3} = R_RTN_to_ECI;
            end
        end
        
        function [a_dot, e_dot, i_dot] = computeElementRates(state, thrust, a)
            % Compute rates of change for orbital elements
            r = norm(state(1:3));
            v = norm(state(4:6));
            thrust_mag = norm(thrust);
            
            % Semimajor axis rate
            a_dot = 2 * a^2 * thrust_mag * cos(angle(state(4:6), thrust)) / v;
            
            % Eccentricity rate
            e_dot = 0; % Placeholder - needs detailed implementation
            
            % Inclination rate
            i_dot = 0; % Placeholder - needs detailed implementation
        end
        
     function [dx] = dynamicalEquations(t, x, params)
            % Full dynamical equations of motion
            state = x(1:6);
            mass = x(7);
            
            % Compute thrust direction
            [u_thrust] = LowThrustTransfer.computeThrustDirection(t, state, params);
            
            % Thrust acceleration
            thrust_mag = params.POWER * params.EFFICIENCY / ...
                (params.G * params.ISPECIFIC * norm(state(4:6)));
            thrust = thrust_mag * u_thrust;
            
            % Position derivative
            r = norm(state(1:3));
            position_dot = state(4:6);
            
            % Velocity derivative
            velocity_dot = -params.MU * state(1:3) / r^3 + thrust / mass;
            
            % Mass derivative
            mass_dot = -thrust_mag / (params.G * params.ISPECIFIC);
            
            dx = [position_dot; velocity_dot; mass_dot];
        end
        
        function [results] = performTransfer()
            % Perform complete LEO to GEO transfer
            
            % Initial conditions
            r_leo = (6371e3 + 550e3);  % Earth radius + LEO altitude
            v_leo = sqrt(LowThrustTransfer4.MU / r_leo);
            
            initial_state = [r_leo; 0; 0; 0; v_leo; 0; 1200];
            
            % Guidance parameters (as in the paper)
            params = struct(...
                'Ke0', 0.0, ...
                'Ke1', -1.5e-3, ...
                'Ke2', 0.0, ...
                'Ki0', -0.33, ...
                'Ki1', -4.3e-3);
            
            % ODE solver options
            options = odeset('MaxStep', 100);
            
            % Solve transfer trajectory
            [t, state] = ode45(@(t,x) LowThrustTransfer4.dynamicalEquations(t, x, params), ...
                [0 216.3*24*3600], initial_state, options);
            
            % Post-process results
            results.time = t;
            results.state = state;
        end
    end
end

% Utility functions
function [optimalInPlane, optimalOutOfPlane] = computeOptimalSteering(a, e, i, mu, v, flightPathAngle)
    % Compute optimal steering angles
    
    % Optimal in-plane steering for semimajor axis
    optimalInPlane = 0;
    
    % Optimal in-plane steering for eccentricity
    r = norm(r);
    trueAnomaly = 0; % This is a simplification
    optimalEccentricitySteering = atan2(r * sin(trueAnomaly), ...
        2 * a * (e + cos(trueAnomaly)));
    
    % Optimal out-of-plane steering
    longitude = 0; % Simplification
    optimalOutOfPlane = (pi/2) * sign(cos(longitude));
end

function [inPlaneSteering, outOfPlaneSteering] = blendGuidanceControls(...
    t, optimalInPlaneSteer, optimalOutOfPlaneSteer, params)
    % Blend guidance control laws
    
    % In-plane steering blending
    if t <= 120 * 24 * 3600  % First phase
        Ke0 = 0.0;
        Ke1 = -1.5e-3;
        Ke2 = 0.0;
    else  % Second phase
        Ke0 = -0.18;
        Ke1 = -1.0e-2;
        Ke2 = 7.0e-5;
    end
    
    % Weighting function
    Ge = Ke0 + Ke1 * t + Ke2 * t^2;
    
    % Out-of-plane steering
    if t <= 120 * 24 * 3600  % First phase
        Ki0 = -0.33;
        Ki1 = -4.3e-3;
    else  % Second phase
        Ki0 = -0.85;
        Ki1 = -5.2e-3;
    end
    
    % Yaw steering
    Gi = Ki0 + Ki1 * t;
    
    % Compute final steering angles
    inPlaneSteering = optimalInPlaneSteer * Ge;
    outOfPlaneSteering = optimalOutOfPlaneSteer * Gi;
end

function [a, e, i, omega, OMEGA, mu] = orbitalElementsFromState(state)
    % Convert state vector to orbital elements (placeholder implementation)
    % This is a simplified version and would need a full implementation
    
    r = state(1:3);
    v = state(4:6);
    
    % Placeholder calculations - would need more rigorous implementation
    a = norm(r);
    e = norm(v)^2/mu - 1/norm(r);
    i = 0;
    omega = 0;
    OMEGA = 0;
    mu = 3.986e14;
end

