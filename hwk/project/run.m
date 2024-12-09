clear all;
clc;
close all;


% ASEN-6015 FINAL PROJECT
% NICOLE JAMES
LTT = LowThrustTransfer;
% CONSTANTS
%  WATT = kg m^2/s^3  Isp = sec  
% Simulation parameters
thrust_params.power = 10e3;           % 10 kW
thrust_params.Isp = 3300;             % Specific impulse
thrust_params.efficiency = 0.65;      % Engine efficiency
thrust_params.initial_mass = 1200;    % kg
thrust_params.mu = LowThrustTransfer.mu;
thrust_params.initial_tw = 3.4e-2;    % thrust in mN (?)          
thrust_params.thrust_acceleration = thrust_params.initial_tw * LowThrustTransfer.g0 * thrust_params.initial_mass;
thrust_params.mass_flow_rate = 2*thrust_params.efficiency* thrust_params.power/(LowThrustTransfer.g0*thrust_params.Isp)^2; 

% Phase 1 and 2 guidance parameters  
guidance_params.phase1.Ke0 = 0;
guidance_params.phase1.Ke1 = -1.5e-3;
guidance_params.phase1.Ke2 = 0.0;
guidance_params.phase1.Ki0 = -0.33;
guidance_params.phase1.Ki1 = -4.3e-3;

guidance_params.phase2.Ke0 = -0.18;
guidance_params.phase2.Ke1 = -1.0e-2;
guidance_params.phase2.Ke2 = 7e-5;
guidance_params.phase2.Ki0 = -0.85;
guidance_params.phase2.Ki1 = -5.2e-3;

% INITIAL CONDITIONS

initial_altitude = 550e3;  
initial_inclination = deg2rad(28.5);
r_earth = 6378e3;  % Earth radius
r_initial = [r_earth + initial_altitude; 0; 0 ];
v_initial = [0; sqrt(LTT.mu / norm(r_initial)); 0];

% Initial state - assume simulation starts at RAAN
rtn2eci = Euler313(0, initial_inclination, 0, 0);
r0 = rtn2eci * r_initial;
v0 = rtn2eci * v_initial;

% TARGET CONDITIONS

% Target GEO conditions
target_semimajor_axis = 42164e3;
vf = sqrt(LTT.mu/target_semimajor_axis);

results = LTT.simulate_transfer(r0, v0, target_semimajor_axis, thrust_params, guidance_params);
