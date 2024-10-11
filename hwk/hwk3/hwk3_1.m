% ASEN-6015 Hwk 3 

% a) gamma, a, b neq 0
% b) a, b = 0
% c) gamma neq 0, a or b new 0

clear
close all
clc
a = 1;              % lower a results in lower J and better response
b = 1               % lower b results in lower J, but lower response
q = 100;            % lower q results in lower J
x0 = [1 -1 1 -1 1]';
dt = 0.01;
tsim = 1:dt:10;
u = zeros(2,length(tsim));
u(:,1) = 1;         % impulse function

cost = @(x,y) y.*x.^2;
% integral J = 1/3 x3^3
A = [5 0 3 0 1;
     3 0 0 -2 0;
     0 -2 4 1 0;
     1 3 -4 1 3;
     0 2 2 0 -1];

B = [0 1; 
     0 2;
     0 0;
     1 3;
     1 1];

C = [0 0 1 0 0];
D = [];
Q = [0 0 0 0 0;
    0 0 0 0 0;
    0 0 q 0 0;
    0 0 0 0 0;
    0 0 0 0 0];

R = [a 0; 0 b];

i=1;
[K,P,poles] = lqr(A,B,Q,R);

H = A-B*K;
eig(H)      % is closed loop system symptotically stable?

S = Q + K'*R*K;
eig(S)      % is S PSD or PD?
sys1 = ss(H,B,C,D);

J(i) = x0'*P*x0       % minimize J by adjusting P/K

[ycalc, tout, xcalc] = initial(sys1,x0,tsim);

figure(1)
plot(tout,xcalc)
title('LQR Controller a=1, b=1, \gamma=100')
xlabel('Time (s)');
ylabel('States');
legend('x1','x2','x3','x4','x5');

x = x0;

for i = 1:length(tsim)
    xdot(:,i) = H*x(:,i) + B*u(:,i);
    x(:,i+1) = x(:,i) + xdot(:,i)*dt;
    y(i) = C * x(:,i);
end

figure(2)
plot(tsim,x(:,1:end-1))
title('LQR Controller States a=1, b=1, \gamma=100')
xlabel('Time (s)');
ylabel('States');
legend('x1','x2','x3','x4','x5');
 
figure(3)
plot(tsim,y)
title('LQR Controller Measurements a=1, b=1, \gamma=100')
xlabel('Time (s)');
ylabel('States');


% Problem 2

L = eye(5);
Rwi = [0.1 0.2 0.1 0.3 0.2];
wt = randn(5)*Rwi'

% use method of adjoints


tbw = 10:-dt:0;

Kf= S;      % K(tf,tf) = S
kdot = @(A,B,Q,R,K) -K*A + K*B*inv(R)*B'*K - Q - K*A';

[tout,kout] = ode45(@(tout,kout) kdot(A,B,Q,R,K),tbw,Kf);

% Integrate the variational equations with control
%[t, state] = ode45(@(t, state) gauss_variational_equations_with_control(t, state, xref, K), tspan, initial_state);

