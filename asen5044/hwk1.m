A MATLAB Script
%% ASEN5044 Assignment 1 Solutions
% Y. Shen
% 31 August 2018
%% Problem 3
% Part B
Ix = 500; % kg mˆ3
Iy = 750; % kg mˆ3
Iz = 1000; % kg mˆ3
p0 = 20; % rad/s
dt = 0.1;
A = [0 0 0; 0 0 p0*(Ix - Iz)/Iy; 0 p0*(Iy - Ix)/Iz 0];
Phi = expm(A*dt);
% Part C
t = 0:0.01:5; % Times (s)
x0 = [0; 0.1; 0]; % Initial conditions (rad/s)
x = zeros(3, length(t)); % States (rad/s)
for i = 1:length(t)
x(:, i) = expm(A*t(i))*x0;
end
figure;
plot(t, x);
xlabel('Time t (s)');
ylabel('Angular Rate Perturbation (rad/s)');
title('Angular Rate Perturbations vs. Time, Using expm');
legend('Roll Perturbation \Deltap', 'Pitch Perturbation \Deltaq', 'Yaw Perturbation \Deltar', 'location', 'best');
grid on;
%% Advanced Question 1
c = p0*(Iz - Ix)/Iy;
b = p0*(Iy - Ix)/Iz;
Phi = @(t) [0 0 0;
0 cos(sqrt(c*b)*t) -sqrt(c/b)*sin(sqrt(c*b)*t);
0 sqrt(b/c)*sin(sqrt(c*b)*t) cos(sqrt(c*b)*t)];
for i = 1:length(t)
x(:, i) = Phi(t(i))*x0;
end
figure;
6
plot(t, x);
xlabel('Time t (s)');
ylabel('Angular Rate Perturbation (rad/s)');
title('Angular Rate Perturbations vs. Time, Using Analytical Solution');
legend('Roll Perturbation \Deltap', 'Pitch Perturbation \Deltaq', 'Yaw Perturbation \Deltar', 'location', 'best');
grid on;