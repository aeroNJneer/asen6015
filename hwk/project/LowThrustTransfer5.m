% Low-Thrust Orbit Transfer Guidance Simulation
% Based on Kluever (1998) Guidance Scheme for LEO to GEO Transfer

classdef LowThrustTransfer < handle
    properties (Constant)
        % Physical constants
        mu = 3.986004418e14;  % Earth gravitational constant, m^3/s^2
        g0 = 9.80665;         % Standard gravity, m/s^2
        
    end
    properties
        % Spacecraft parameters
        initial_mass
        Isp
        efficiency
        initial_thrust_to_weight
        max_power
        
        % Initial orbital conditions
        r0
        v0
        initial_inclination
        initial_RAAN
        initial_omega
        
        % Target orbital conditions
        target_semimajor_axis
        target_eccentricity
        target_inclination
    end
    
    methods
        function obj = LowThrustTransfer()
            % Constructor with default parameters from the paper
            obj.initial_mass = 1200;           % kg
            obj.Isp = 3300;                    % s
            obj.efficiency = 0.65;             % 65%
            obj.max_power = 10000;             % W
            
            % Initial LEO conditions
            alt = 550e3;                       % 550 km altitude
            obj.r0 = (6378e3 + alt) * [1; 0; 0];
            obj.initial_inclination = deg2rad(28.5);
            obj.initial_RAAN = 0;
            obj.initial_omega = 0;
            
            % Compute initial orbital velocity
            v_circ = sqrt(obj.mu / norm(obj.r0));
            
            % Rotation matrices for initial orbital plane
            Rz_omega = [cos(obj.initial_omega), -sin(obj.initial_omega), 0;
                        sin(obj.initial_omega),  cos(obj.initial_omega), 0;
                        0, 0, 1];
            
            Rx_inc = [1, 0, 0;
                      0,  cos(obj.initial_inclination), -sin(obj.initial_inclination);
                      0,  sin(obj.initial_inclination),  cos(obj.initial_inclination)];
            
            Rz_RAAN = [cos(obj.initial_RAAN), -sin(obj.initial_RAAN), 0;
                        sin(obj.initial_RAAN),  cos(obj.initial_RAAN), 0;
                        0, 0, 1];
            
            % Initial velocity vector
            R_total = Rz_RAAN * Rx_inc * Rz_omega;
            v0_orbit_plane = [0; v_circ; 0];
            obj.v0 = R_total * v0_orbit_plane;
            
            % Target orbital conditions
            obj.target_semimajor_axis = 42164e3;
            obj.target_eccentricity = 0;
            obj.target_inclination = 0;
        end
        
        function [u_RTN, steering_angles] = compute_steering(obj, t, r, v, mass)
            % Compute steering angles based on time-varying guidance parameters
            
            % Guidance parameter phases as described in the paper
            if t <= 120 * 86400  % First 120 days
                Ke0 = 0.0;
                Ke1 = -1.5e-3;
                Ke2 = 0.0;
                Ki0 = -0.33;
                Ki1 = -4.3e-3;
            else  % Remainder of transfer
                Ke0 = -0.18;
                Ke1 = -1.0e-2;
                Ke2 = 7.0e-5;
                Ki0 = -0.85;
                Ki1 = -5.2e-3;
            end
            
            % Compute orbital elements
            oev = eci2orb1(obj.mu, r, v);
            
            % Compute flight path angle
            h_vec = cross(r, v);
            flight_path_angle = asin(dot(r, v) / (norm(r) * norm(v)));
            
            % Compute true anomaly and argument of periapsis
            true_anomaly = acos(dot(r, v) / (norm(r) * norm(v)));
            theta = oev.omega + true_anomaly;
            
            % Optimal in-plane steering (semimajor axis and eccentricity)
            c_a = [sin(0); cos(0); 0];  % Tangent steering
            c_e = [sin(atan2(norm(r)*sin(theta), 2*oev.a*(oev.e+cos(theta)))); 
                   cos(atan2(norm(r)*sin(theta), 2*oev.a*(oev.e+cos(theta)))); 
                   0];
            
            % Blending of in-plane steering
            Ge = Ke0 + Ke1*t/86400 + Ke2*(t/86400)^2;
            c_blended = (Ge*c_e + c_a) / norm(Ge*c_e + c_a);
            
            % Compute pitch and yaw angles
            pitch_angle = atan2(c_blended(1), c_blended(2));
            
            % Yaw steering
            Gi = Ki0 + Ki1*t/(86400);
            yaw_angle = Gi * (pi/2 * cos(theta));
            
            % Thrust direction unit vector in RTN frame
            u_RTN = [sin(pitch_angle)*cos(yaw_angle);
                     cos(pitch_angle)*cos(yaw_angle);
                     sin(yaw_angle)];
            
            steering_angles.pitch = pitch_angle;
            steering_angles.yaw = yaw_angle;
        end
        
        function [r_dot, v_dot, m_dot] = compute_dynamics(obj, t, r, v, mass, u_RTN)
            % Compute spacecraft dynamics in ECI frame
            
            % Rotation from RTN to ECI
            [a, e, i, omega, RAAN, true_anomaly, mu] = cart2orb(r, v);
            R_RTN_to_ECI = compute_RTN_to_ECI_rotation(r, v);
            
            % Compute thrust acceleration in ECI frame
            thrust_acceleration_magnitude = obj.compute_thrust_acceleration(mass);
            u_ECI = R_RTN_to_ECI * u_RTN;
            
            % Compute acceleration due to thrust and gravity
            thrust_accel = thrust_acceleration_magnitude * u_ECI;
            gravity_accel = -mu * r / norm(r)^3;
            
            % State derivatives
            r_dot = v;
            v_dot = gravity_accel + thrust_accel;
            
            % Mass rate (rocket equation)
            m_dot = -thrust_acceleration_magnitude * norm(u_ECI) / (obj.Isp * obj.g0);
        end
        
        function thrust_accel = compute_thrust_acceleration(obj, mass)
            % Compute thrust acceleration based on power and mass
            max_thrust = obj.max_power / (obj.Isp * obj.g0) * obj.efficiency;
            thrust_accel = max_thrust / mass;
        end

    end
end


