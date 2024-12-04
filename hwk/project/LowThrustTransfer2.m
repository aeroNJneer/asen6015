% Low-Thrust Orbital Transfer Simulation
% Based on Kluever 2012 "Simple Guidance Scheme for Low-Thrust Orbit Transfers"

classdef LowThrustTransfer
    properties (Constant)
        % Physical constants
        mu = 3.986e14;         % Earth's gravitational parameter (m^3/s^2)
        g0 = 9.81;             % Gravitational acceleration (m/s^2)
    end
    
    methods (Static)
        % Compute optimal in-plane thrust vector for semimajor axis
        function c_a = computeAxialControlVector(v, r, flight_path_angle)
            c_a = [sin(flight_path_angle); cos(flight_path_angle); 0];
        end
        
        % Compute optimal in-plane thrust vector for eccentricity
        function c_e = computeEccentricityControlVector(a, e, r, true_anomaly, flight_path_angle)
            tan_angle = r * sin(true_anomaly) / (2 * a * (e + cos(true_anomaly)));
            c_e = [sin(flight_path_angle + atan(tan_angle)); 
                   cos(flight_path_angle + atan(tan_angle)); 
                   0];
        end
        
        % Compute out-of-plane thrust steering
        function yaw = computeYawSteering(longitude_angle)
            yaw = (pi/2) * cos(longitude_angle);
        end
        
        % Blended thrust direction computation
        function [thrust_unit, pitch_angle] = computeThrustDirection(c_a, c_e, Ge, flight_path_angle)
            % Normalize and blend control vectors
            c_blended = Ge * c_e + c_a;
            c_blended = c_blended / norm(c_blended);
            
            % Extract pitch angle
            pitch_angle = atan2(c_blended(1), c_blended(2));
        end
        
        % Guidance parameter schedule
        function [Ge, Gi] = computeWeightingFunctions(t, phase_params)
            % Quadratic weighting function for eccentricity
            Ge = phase_params.Ke0 + phase_params.Ke1 * t + phase_params.Ke2 * t^2;
            
            % Linear weighting for out-of-plane steering
            Gi = phase_params.Ki0 + phase_params.Ki1 * t;
        end
        
        % Spacecraft dynamics
        function [state_dot] = spacecraft_dynamics(t, state, thrust_params)
            % State vector: [x, y, z, vx, vy, vz, mass]
            r = state(1:3);
            v = state(4:6);
            m = state(7);
            
            % Gravitational acceleration
            r_mag = norm(r);
            a_grav = -thrust_params.mu * r / (r_mag^3);
            
            % Compute current thrust direction
            % This would be determined by the guidance laws
            u_thrust = [1; 0; 0];  % Placeholder, actual direction computed dynamically
            
            % Thrust acceleration
            a_thrust = thrust_params.thrust_acceleration * u_thrust / m;
            
            % Mass rate change
            mdot = -thrust_params.mass_flow_rate;
            
            % State derivative
            state_dot = [
                v;
                a_grav + a_thrust;
                mdot
            ];
        end
        
        % Main transfer simulation
        function results = simulate_transfer()
            % Initial conditions (matching paper's scenario)
            initial_altitude = 550e3;  % 550 km
            initial_inclination = 28.5 * pi/180;
            
            % Target GEO conditions
            target_semimajor_axis = 42164e3;
            
            % Simulation parameters
            thrust_params.power = 10e3;           % 10 kW
            thrust_params.Isp = 3300;             % Specific impulse
            thrust_params.efficiency = 0.65;      % Engine efficiency
            thrust_params.initial_mass = 1200;    % kg
            
            % Compute derived parameters
            thrust_params.mu = LowThrustTransfer.mu;
            thrust_params.thrust_acceleration = ...
                thrust_params.power * thrust_params.efficiency / ...
                (thrust_params.Isp * LowThrustTransfer.g0);
            thrust_params.mass_flow_rate = ...
                thrust_params.thrust_acceleration * thrust_params.initial_mass / ...
                (thrust_params.Isp * LowThrustTransfer.g0);
            
            % Phase 1 guidance parameters (first 120 days)
            phase1_params.Ke0 = 0.0;
            phase1_params.Ke1 = -1.5e-3;
            phase1_params.Ke2 = 0.0;
            phase1_params.Ki0 = -0.33;
            phase1_params.Ki1 = -4.3e-3;
            
            % Phase 2 guidance parameters (remaining transfer)
            phase2_params.Ke0 = -0.18;
            phase2_params.Ke1 = -1.0e-2;
            phase2_params.Ke2 = 7e-5;
            phase2_params.Ki0 = -0.85;
            phase2_params.Ki1 = -5.2e-3;
            
            % Simulation setup
            options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
            tspan = [0, 216.3 * 24 * 3600];  % Transfer time in seconds
            
            % Initial state
            r0 = [r_initial, 0, 0];
            v0 = [0, sqrt(thrust_params.mu / r_initial), 0];
            
            % Solve ODE
            [t, state] = ode45(@(t,y) LowThrustTransfer.spacecraft_dynamics(t,y,thrust_params), ...
                tspan, [r0, v0, thrust_params.initial_mass], options);
            
            % Process and return results
            results.time = t;
            results.state = state;
        end
    end
end

% Example usage
transfer_results = LowThrustTransfer.simulate_transfer();