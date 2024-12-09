
% Main simulation script
function simulate_low_thrust_transfer()
close all;
clear;
clc;
    % Initialize transfer
    transfer = LowThrustTransfer();
    
    % Simulation parameters
    t_start = 0;
    t_end = 250 * 86400;  % 250 days
    dt = 100;  % Time step
    
    % Initialize state
    r = transfer.r0;
    v = transfer.v0;
    mass = transfer.initial_mass;
    
    % Storage for results
    results.time = [];
    results.position = [];
    results.velocity = [];
    results.mass = [];
    results.orbital_elements = [];
    results.steering = [];
    
    % Numerical integration (simple Euler method for demonstration)
    for t = t_start:dt:t_end
        % Compute steering
        [u_RTN, steering_angles] = transfer.compute_steering(t, r, v, mass);
        
        % Compute dynamics
        [r_dot, v_dot, m_dot] = transfer.compute_dynamics(t, r, v, mass, u_RTN);
        
        % Update state
        r = r + r_dot * dt;
        v = v + v_dot * dt;
        mass = mass + m_dot * dt;
        
        % Compute and store orbital elements
        [a, e, i, ~, ~, ~] = cart2orb(r, v);
        
        % Store results
        results.time(end+1) = t;
        results.position(:,end+1) = r;
        results.velocity(:,end+1) = v;
        results.mass(end+1) = mass;
        results.orbital_elements(:,end+1) = [a; e; i];
        results.steering(end+1).pitch = steering_angles.pitch;
        results.steering(end+1).yaw = steering_angles.yaw;
    end
    
    % Plotting results
    plot_results(results);
end

function plot_results(results)
    % Plotting orbital elements
    figure;
    subplot(3,1,1);
    plot(results.time/86400, results.orbital_elements(1,:));
    title('Semimajor Axis');
    xlabel('Time (days)'); ylabel('a (km)');
    
    subplot(3,1,2);
    plot(results.time/86400, results.orbital_elements(2,:));
    title('Eccentricity');
    xlabel('Time (days)'); ylabel('e');
    
    subplot(3,1,3);
    plot(results.time/86400, rad2deg(results.orbital_elements(3,:)));
    title('Inclination');
    xlabel('Time (days)'); ylabel('i (degrees)');
    
    % Plotting steering angles
    figure;
    subplot(2,1,1);
    plot(results.time(1:length(results.steering))/86400, ...
         arrayfun(@(x) x.pitch, results.steering));
    title('Pitch Steering Angle');
    xlabel('Time (days)'); ylabel('Angle (rad)');
    
    subplot(2,1,2);
    plot(results.time(1:length(results.steering))/86400, ...
         arrayfun(@(x) x.yaw, results.steering));
    title('Yaw Steering Angle');
    xlabel('Time (days)'); ylabel('Angle (rad)');
end
