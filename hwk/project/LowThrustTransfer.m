% Low-Thrust Orbit Transfer Guidance Scheme
% Based on the guidance method described by Craig A. Kluever
% Implements a guidance scheme for solar electric propulsion (SEP) spacecraft

classdef LowThrustTransfer
    properties
        % Spacecraft and mission parameters
        initial_orbit       % Initial orbit parameters [a, e, i]
        target_orbit        % Target orbit parameters [a, e, i]
        
        % Propulsion parameters
        power               % Input power (W)
        specific_impulse    % Specific impulse (s)
        engine_efficiency   % Engine efficiency
        
        % Guidance parameters
        Ke0                 % Eccentricity guidance parameter
        Ke1                 % Eccentricity guidance parameter
        Ke2                 % Eccentricity guidance parameter
        Ki0                 % Inclination guidance parameter
        Ki1                 % Inclination guidance parameter
    end
    
    methods
        % Constructor
        function obj = LowThrustTransfer(initial_orbit, target_orbit, power, Isp, efficiency)
            obj.initial_orbit = initial_orbit;
            obj.target_orbit = target_orbit;
            obj.power = power;
            obj.specific_impulse = Isp;
            obj.engine_efficiency = efficiency;
        end
        
        % Compute optimal in-plane thrust steering angle for semimajor axis
        function angle = compute_semimajor_thrust_angle(obj, state)
            % Equation (4): Optimal steering for semimajor axis
            angle = 0;
        end
        
        % Compute optimal in-plane thrust steering angle for eccentricity
        function angle = compute_eccentricity_thrust_angle(obj, state)
            % Equation (5): Optimal steering for eccentricity
            r = norm(state(1:3));
            v = norm(state(4:6));
            true_anomaly = compute_true_anomaly(state);
            a = compute_semimajor_axis(state);
            e = compute_eccentricity(state);
            
            angle = atan2(r * sin(true_anomaly), ...
                          2*a*(e + cos(true_anomaly)));
        end
        
        % Compute optimal out-of-plane (yaw) steering angle
        function angle = compute_yaw_steering_angle(obj, state)
            % Equation (7): Yaw steering law
            longitude_angle = compute_longitude_angle(state);
            angle = (pi/2) * cos(longitude_angle);
        end
        
        % Compute guidance weighting functions
        function [Ge, Gi] = compute_guidance_weights(obj, t, phase)
            % Equation (12): Quadratic weighting function
            if phase == 1
                % First phase parameters (0 ≤ t ≤ 120 days)
                Ke0 = 0.0;
                Ke1 = -1.5e-3;
                Ke2 = 0.0;
                Ki0 = -0.33;
                Ki1 = -4.3e-3;
            else
                % Second phase parameters (120 < t ≤ transfer time)
                Ke0 = -0.18;
                Ke1 = -1.0e-2;
                Ke2 = 7.0e-5;
                Ki0 = -0.85;
                Ki1 = -5.2e-3;
            end
            
            % Compute eccentricity weighting
            Ge = Ke0 + Ke1*t + Ke2*t^2;
            
            % Compute inclination weighting
            Gi = Ki0 + Ki1*t;
        end
        
        % Compute thrust direction
        function thrust_direction = compute_thrust_direction(obj, state, t)
            % Compute thrust direction in RTN frame
            flight_path_angle = compute_flight_path_angle(state, mu);
            
            % Determine mission phase
            if t<=120
                phase = 1;
            else
                phase = 2;
            end
            
            % Compute guidance weights
            [Ge, Gi] = obj.compute_guidance_weights(t, phase);
            
            % Compute optimal thrust angles
            angle_semimajor = obj.compute_semimajor_thrust_angle(state);
            angle_eccentricity = obj.compute_eccentricity_thrust_angle(state);
            yaw_angle = obj.compute_yaw_steering_angle(state);
            
            % Blend in-plane thrust direction vectors
            ca = [sin(angle_semimajor); cos(angle_semimajor); 0];
            ce = [sin(angle_eccentricity); cos(angle_eccentricity); 0];
            
            c_blended = (Ge*ca + Ge*ce) / norm(Ge*ca + Ge*ce);
         
            % Compute pitch and yaw steering angles
            pitch_angle = atan2(c_blended(1), c_blended(2));
 
            % scale the amplitude of yaw steering
            yaw_angle = Gi * yaw_angle;
   
            % Compute final thrust unit vector
            thrust_direction = [
                sin(pitch_angle) * cos(yaw_angle);
                cos(pitch_angle) * cos(yaw_angle);
                sin(yaw_angle)
            ];
        end
        
        % Perform low-thrust orbit transfer simulation
        function [trajectory, transfer_time] = perform_transfer(obj)
            % Implement numerical integration of spacecraft trajectory
            % This is a placeholder for a full numerical integration routine
            % In practice, this would use Runge-Kutta or another numerical method
            
            % Initialization
            state = initialize_spacecraft_state(obj.initial_orbit);
            t = 0;
            max_transfer_time = 300;  % days
            % calculate current period and simulate over 100 points
            while ~orbit_transfer_complete(state, obj.target_orbit) && t < max_transfer_time
                % Compute thrust direction
                thrust_dir = obj.compute_thrust_direction(state, t);
                
                % Update spacecraft state (numerical integration)
                [ts, state] = ode45(@(t,state) update_spacecraft_state(state, thrust_dir, t), [0 1], x0);
                t = t + time_step;
            end
            
            trajectory = state;
            transfer_time = t;
        end
    end
end

% Helper functions (would need detailed implementation)
function true_anom = compute_true_anomaly(state, mu)
    % Compute true anomaly from state vector
    r = state(1:3);
    v = state(4:6);
    H = cross(r,v);
    true_anom =  deg2rad(acosd((H^2/(mu*r) -1)/e));
end

function semimajor_axis = compute_semimajor_axis(state, mu)
    % Compute semimajor axis from state vector
    r = norm(state(1:3));
    v = norm(state(4:6));
    semimajor_axis = -mu/(v^2 -(2*mu/r));
end

function ecc = compute_eccentricity(state, mu)
    % Compute orbit eccentricity from state vector
    r = state(1:3);
    v = state(4:6);
    H = cross(r,v);
    Ecc = cross(v,H)/mu - r/norm(r);
    ecc = norm(Ecc);
end

function long_angle = compute_longitude_angle(state)
    % Compute longitude angle from state vector

end

function flight_path = compute_flight_path_angle(state, mu)
    % Compute flight path angle from state vector
    e = compute_eccentricity(state, mu);
    TA  = compute_true_anomaly(state, mu);
    flight_path = atan2((e*sin(TA)),(1+(e*cos(TA))));
end

function complete = orbit_transfer_complete(state, target_orbit)
    % Check if orbit transfer is complete based on target orbit parameters

end

function [t, new_state] = update_spacecraft_state(t, state, thrust_dir)
    % Numerical integration of spacecraft dynamics
    % Implement Runge-Kutta or other numerical integration method
    % Unpack state
    r = state(1:3);
    v = state(4:6);
    mdot = -2*eta*P/(g*Isp)^2;

%     1) Compute the unit vectors ca and ce that define the in-plane optimal
% controls for maximum da=dt and de=dt , respectively. The unit
% vectors are expressed in a local rotating radial-transverse-normal
% (RTN) coordinate frame where the R axis is along the radial direction,
% the T axis is in the orbit plane along the transverse direction,
% and the N axis is normal to the orbit plane (however, the in-plane
% steering unit vectors ca and ce will only have R and T components).
% A general expression for the unit vector in the RTN frame is
% ck D sin ° C Á¤k ; cos ° C Á¤k ; 0
% T
% ; k D a; e (8)
% where Á¤a and Á¤e are computed using the optimal controls (4) and
% (5). The angle ° is the  ight-path angle and is measured from the
% local horizon to the velocity vector.

    % Gravitational acceleration
    g_eci = -mu*r/norm(r)^3;
    
    % State derivatives
    dr = v;
    dv = g_eci + acc*thrust_dir;
    
    dxdt = [dr; dv; mdot];
end

function initial_state = initialize_spacecraft_state(initial_orbit)
    % Convert orbital elements to initial state vector
    [r, v] = oe2rv()
    
end

function [a, e, i, Omega, omega, nu] = rv2oe(r, v, target)
    % Convert state vectors to orbital elements
    global mu 
    h = cross(r, v);
    n = cross([0; 0; 1], h);
    e_vec = cross(v, h)/mu - r/norm(r);
    e = norm(e_vec);
    a = 1/(2/norm(r) - dot(v,v)/mu);
    
    i = acos(h(3)/norm(h));
    Omega = acos(n(1)/norm(n));
    if n(2) < 0
        Omega = 2*pi - Omega;
    end
    
    omega = acos(dot(n,e_vec)/(norm(n)*e));
    if e_vec(3) < 0
        omega = 2*pi - omega;
    end
    % project x onto orbit ellipse to find nu
    vp_eci = norm(target.vp_eci);
    rp_eci = norm(target.rp_eci);
    nu = atan2(dot(r,(target.vp_eci/vp_eci)),dot(r,(target.rp_eci/rp_eci)));
end

function [r, v] = orbital_elements_to_state_vector(a, e, i, RAAN, omega, nu, mu)
    % Convert orbital elements to state vector
    % 
    % Inputs:
    %   a       - semimajor axis (km)
    %   e       - eccentricity (dimensionless)
    %   i       - inclination (radians)
    %   RAAN    - Right Ascension of Ascending Node (radians)
    %   omega   - Argument of periapsis (radians)
    %   nu      - True anomaly (radians)
    %   mu      - Gravitational parameter (km^3/s^2), default is Earth's if not provided
    %
    % Outputs:
    %   r       - Position vector in ECI frame (km)
    %   v       - Velocity vector in ECI frame (km/s)

    % Default to Earth's gravitational parameter if not specified
    if nargin < 7
        mu = 3.986004418e5;  % Earth's standard gravitational parameter (km^3/s^2)
    end

    % Calculate p (semi-latus rectum)
    p = a * (1 - e^2);

    % Calculate position in orbital (perifocal) frame
    r_pqw = [p * cos(nu) / (1 + e * cos(nu));
        p * sin(nu) / (1 + e * cos(nu));
        0];

    % Calculate velocity in orbital (perifocal) frame
    v_pqw = [-sqrt(mu/p) * sin(nu);
        sqrt(mu/p) * (e + cos(nu));
        0];

    DCM = compute_Euler313(i, RAAN, omega, nu);

    % Transform position and velocity to ECI frame
    r = DCM * r_pqw;
    v = DCM * v_pqw;
end

function [a, e, i, RAAN, omega, nu] = state_vector_to_orbital_elements(r, v, mu)
    % Convert state vector to orbital elements
    
    % Default to Earth's gravitational parameter if not specified
    if nargin < 3
        mu = 3.986004418e5;  % Earth's standard gravitational parameter (km^3/s^2)
    end

    % Angular momentum vector
    h = cross(r, v);

    % Eccentricity vector
    e_vec = cross(v, h)/mu - r/norm(r);

    % Eccentricity magnitude
    e = norm(e_vec);

    % Node vector (perpendicular to orbital plane)
    n = cross([0, 0, 1], h);

    % Semimajor axis
    r_mag = norm(r);
    v_mag = norm(v);
    a = 1 / (2/r_mag - v_mag^2 / mu);

    % Inclination
    i = acos(h(3) / norm(h));

    % Right Ascension of Ascending Node
    if norm(n) ~= 0
        RAAN = atan2(n(2), n(1));
        if RAAN < 0
            RAAN = RAAN + 2*pi;
        end
    else
        RAAN = 0;
    end

    % Argument of periapsis
    if norm(n) ~= 0
        omega = atan2(dot(e_vec, cross(n, r/r_mag)), dot(n, e_vec));
        if omega < 0
            omega = omega + 2*pi;
        end
    else
        omega = 0;
    end

    % True anomaly
    nu = atan2(dot(h, cross(r, v)), dot(r, v));
    if nu < 0
        nu = nu + 2*pi;
    end
end


function DCM_COE_to_ECI = compute_Euler313(i, RAAN, w, nu)
    % use Euler 3-1-3 rotation matrix to find state vector
    % rotation matrices
    M1 = @(a) [1 0 0; 0 cosd(a) sind(a); 0 -sind(a) cosd(a)];
    M3 = @(a) [cosd(a) sind(a) 0; -sind(a) cosd(a) 0; 0 0 1];
    
    % DCM  rth to xyz *** theta = trueAnom + arg of periapsis 
    DCM = @(Raan, i, theta) M3(-Raan) * M1(-i) * M3(-(nu+w));
end


function DCM_RTN_to_ECI = compute_RTN_to_ECI_matrix(r, v)
    % r: position vector in ECI frame
    % v: velocity vector in ECI frame
    % Radial unit vector (R) = normalized position vector
    R = r / norm(r);    
    % Angular momentum vector (H) = cross product of r and v
    H = cross(r, v);    
    % Normal unit vector (N) = normalized angular momentum vector
    N = H / norm(H);    
    % Transverse unit vector (T) = cross product of N and R
    T = cross(N, R);
    
    % Direction Cosine Matrix (rows are RTN basis, columns are ECI basis)
    DCM_RTN_to_ECI = [R', T', N']';
end