% Low-Thrust Orbital Transfer Simulation
% Based on Kluever 2012 "Simple Guidance Scheme for Low-Thrust Orbit Transfers"

classdef LowThrustTransfer
    properties (Constant)        % Physical constants
        mu = 3.986e14;         % Earth's gravitational parameter (m^3/s^2)
        g0 = 9.81;             % Gravitational acceleration (m/s^2)
        tgo = 216;                % time to go (simulation time)
        dt = 3600;           % time between iterations 
    end

    
    methods (Static)
        % Compute RTN to ECI transformation matrix
        function R_RTN2ECI = compute_RTN_to_ECI_matrix(r, v)
            % Radial unit vector (R direction)
            r_unit = r / norm(r);
            
            % Velocity unit vector for transverse direction
            v_unit = v / norm(v);
            
            % Normal unit vector (cross product of r and v)
            n_unit = cross(r, v) / norm(cross(r, v));
            
            % Transverse unit vector (completes right-handed system)
            t_unit = cross(n_unit, r_unit);
            
            % Transformation matrix from RTN to ECI
            R_RTN2ECI = [r_unit, t_unit, n_unit]';
        end
        
        % Compute orbital elements
        coes = eci2orb1(mu, r, v);
        
        % Compute flight path and longitude angles
        function [flight_path_angle, longitude_angle] = compute_path_angles(r, v)
            % Flight path angle
            flight_path_angle = atan2(v(1),v(2));
            
            % Compute longitude angle
            coes = eci2orb1(LowThrustTransfer.mu,r,v);
            longitude_angle = coes.omega + coes.true_anomaly;
            if longitude_angle > 2*pi
                longitude_angle = longitude_angle - 2*pi';
            end
        end
        
        % Compute optimal in-plane thrust vector for semimajor axis
        function c_a = computeAxialControlVector(flight_path_angle)
            c_a = [sin(flight_path_angle); cos(flight_path_angle); 0];
        end
        
        % Compute optimal in-plane thrust vector for eccentricity
        function c_e = computeEccentricityControlVector(a, e, r, angle, true_anomaly, flight_path_angle)
            % if e ~0 then true_anomaly passed in will be 0 and angle used
            % will be argument of latitude, otherwise true anomaly will be
            % passed in for both angle and true_anomaly
            phi_e = atan(r * sin(angle) / (2 * a * (e + cos(true_anomaly))));
            c_e = [sin(flight_path_angle + phi_e); 
                   cos(flight_path_angle + phi_e); 
                   0];
            if ~isreal(c_e)
                c_e;
            end
        end
        
        % Compute out-of-plane thrust steering
        function yaw = computeYawSteering(longitude_angle)
            yaw = (pi/2) * cos(longitude_angle);
        end
        
        % Blended thrust direction computation in RTN frame
        function [thrust_unit_RTN, pitch_angle, yaw] = computeThrustDirection(c_a, c_e, Ge, Gi, flight_path_angle, longitude_angle)
            % Normalize and blend control vectors
            c_blended = Ge * c_e + c_a;
            c_blended = c_blended / norm(c_blended);
            
            % Extract pitch angle
            pitch_angle = atan2(c_blended(1), c_blended(2));
            
            % Compute thrust direction in RTN frame with yaw
            yaw = Gi*LowThrustTransfer.computeYawSteering(longitude_angle); 
            if yaw >2*pi
                yaw = 2*pi - yaw;
            elseif yaw < -2*pi
                yaw = 2*pi + yaw;
            end
            thrust_unit_RTN = [
                sin(pitch_angle) * cos(yaw);
                cos(pitch_angle) * cos(yaw);
                sin(yaw)
            ];
        end
        
        % Guidance parameter schedule
        function [Ge, Gi] = computeWeightingFunctions(t, phase_params)
            % Quadratic weighting function for eccentricity
            Ge = phase_params.Ke0 + phase_params.Ke1 * t + phase_params.Ke2 * t^2;
            
            % Linear weighting for out-of-plane steering
            Gi = phase_params.Ki0 + phase_params.Ki1 * t;
        end
        
        % Spacecraft dynamics
        function [state_dot, angles] = spacecraft_dynamics(t, state, thrust_params, guidance_params)
            % State vector: [x, y, z, vx, vy, vz, mass]
            r = state(1:3);
            v = state(4:6);
            m = state(7);
            
            % Compute current orbital state
            [flight_path_angle, longitude_angle] = LowThrustTransfer.compute_path_angles(r, v);
            coes = eci2orb1(LowThrustTransfer.mu,r, v);
            
            % Determine current phase and guidance parameters
            if t <= 120 * 24 * 3600
                phase_params = guidance_params.phase1;
            else
                phase_params = guidance_params.phase2;
            end
            
            % Compute weighting functions
            [Ge, Gi] = LowThrustTransfer.computeWeightingFunctions(t, phase_params);
            
            % Compute control vectors
            c_a = LowThrustTransfer.computeAxialControlVector(flight_path_angle);
            c_e = LowThrustTransfer.computeEccentricityControlVector(coes.a,coes.e,norm(r),coes.true_anomaly, coes.true_anomaly, flight_path_angle);

            % Compute thrust direction in RTN frame
            [thrust_unit_RTN, pitch_angle, yaw] = LowThrustTransfer.computeThrustDirection(c_a, c_e, Ge, Gi, flight_path_angle, longitude_angle);
            
            % Transform thrust direction to ECI frame
            R_RTN2ECI = LowThrustTransfer.compute_RTN_to_ECI_matrix(r, v);
            a_rtn2eci = Euler313(coes.raan, coes.inc, coes.omega, coes.true_anomaly);
            thrust_unit_ECI = R_RTN2ECI' * thrust_unit_RTN;
            
            % Gravitational acceleration
            r_mag = norm(r);
            a_grav = -thrust_params.mu * r / (r_mag^3);
            
            % Thrust acceleration
            a_thrust = thrust_params.thrust_acceleration * thrust_unit_ECI / m;
            

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
        function results = simulate_transfer(r0, v0, target_semi_major_axis, thrust_params, guidance_params)
            % Simulation setup
            options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
            dt = LowThrustTransfer.dt;
            todays = 24*3600;
            tspan = 0:dt:LowThrustTransfer.tgo*todays;  % Transfer time in seconds
           % tspan = 0:dt:Tof
            
            % Solve ODE
            [t, state] = ode45(@(t,y) LowThrustTransfer.spacecraft_dynamics(t, y, thrust_params, guidance_params), ...
                tspan, [r0; v0; thrust_params.initial_mass], options);
            
            % Process and return results
            results.time = t;
            results.state = state;
            
            for i = 1:length(state)
                coes(i) = eci2orb1(LowThrustTransfer.mu,state(i,1:3), state(i,4:6));
            end
            results.a = [coes.a];
            results.ecc = [coes.e];
            results.inc = [coes.inc];
            results.true_anomaly = [coes.true_anomaly];

            % Plot results
            figure;
            subplot(2,1,1);
            plot(t/24/3600 , vecnorm(state(:,1:3), 2, 2), 'b-');
            title('Radius vs Time');
            xlabel('Time (days)');
            ylabel('Altitude (m)');
            
            subplot(2,1,2);
            plot(t/24/3600 , state(:,7), 'r-');
            title('Spacecraft Mass vs Time');
            xlabel('Time (days)');
            ylabel('Mass (kg)');

            figure;
            plot3(state(:,1),state(:,2),state(:,3))
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            
            figure;
            plot(t/24/3600,results.ecc,'DisplayName','Eccentricity')
            yyaxis left
            ylabel('Eccentricity')
            hold on
            plot(t/24/3600, results.inc,'--', 'DisplayName','Inclination')
            yyaxis right
            ylabel('Inclination')
            title('Orbit Eccentricty and Inclination')
            xlabel('Semi-Major Axis (km)')
       %     legend('Location','best')
        end
  
        function [r, v] = orbitalElementsToStateVectors(a, e, i, RAAN, w, TA, mu)
            % a    - semi-major axis (km)
            % e    - eccentricity
            % i    - inclination (degrees)
            % RAAN - right ascension of the ascending node (degrees)
            % w    - argument of periapsis (degrees)
            % TA   - true anomaly (degrees)
            % mu   - standard gravitational parameter (km^3/s^2)
        
            % Convert angles from degrees to radians
            i = deg2rad(i);
            RAAN = deg2rad(RAAN);
            w = deg2rad(w);
            TA = deg2rad(TA);
        
            % Compute the distance (r) and speed (v) in the orbital plane
            p = a * (1 - e^2);
            r_orbital = p / (1 + e * cos(TA));
            v_orbital = sqrt(mu * (2/r_orbital - 1/a));
        
            % Position in the orbital plane
            r_orbital_plane = [r_orbital * cos(TA); r_orbital * sin(TA); 0];
        
            % Velocity in the orbital plane
            v_orbital_plane = [-sin(TA); e + cos(TA); 0] * (sqrt(mu/p));
        
            % Rotation matrices
            R3_W = [cos(w), -sin(w), 0; sin(w), cos(w), 0; 0, 0, 1];
            R1_i = [1, 0, 0; 0, cos(i), -sin(i); 0, sin(i), cos(i)];
            R3_RAAN = [cos(RAAN), -sin(RAAN), 0; sin(RAAN), cos(RAAN), 0; 0, 0, 1];
        
            % Combined rotation matrix
            Q = R3_RAAN * R1_i * R3_W;
        
            % Position and velocity in the inertial frame
            r = Q * r_orbital_plane;
            v = Q * v_orbital_plane;
        end
    end
end

