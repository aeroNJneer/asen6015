% ASEN-6015 Hwk 3 

% a) gamma, a, b neq 0
% b) a, b = 0
% c) gamma neq 0, a or b new 0

clear
close all
clc
a = 0;              % lower a results in lower J and better response
b = 1 ;              % lower b results in lower J, but lower response
q = 100;            % lower q results in lower J
x0 = [1 -1 1 -1 1]';
dt = 0.01;
tsim = 1:dt:10;
u = zeros(2,length(tsim));
%u(:,1) = 1;         % impulse function

cost = @(x,y) y.*x.^2;

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
eig(H)  ;    % is closed loop system symptotically stable?

S = Q + K'*R*K;
eig(S)   ;   % is S PSD or PD?
sys1 = ss(H,B,C,D);

J(i) = x0'*P*x0   ;    % minimize J by adjusting P/K

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
% part a
L = eye(5);
Rwi = [0.1 0.2 0.1 0.3 0.2]';
%wt = randn(5)*Rwi';

% use method of adjoints

tbw = 6:-dt:0;  % integrate backwards at t=6 .. t=0
pdot = @(A,p) -A'*p;
p0 = C';
[tout,pout] = ode45(@(t,p) pdot(H,p),tbw,p0); % find history of adjoints
xp = x0;
Ey=zeros(5,1);
% integrate forward in time
for i = 1:length(tout)
    wt = randn(5)*Rwi;     % generate random noise
    xdotp(:,i) = H*xp(:,i) + L*wt;
    xp(:,i+1) = xp(:,i) + xdotp(:,i)*dt;
    yp(i) = pout(i,:) * xp(:,i);
    dEy = pout(i,:)*L*Rwi*L'*pout(i,:)';
    Ey = Ey + dEy*dt;
end

figure(4) 
plot(tout,pout)
title('Adjoints Evolution')
xlabel('Time (s)')
ylabel('Adjoint System State')
legend('p_1','p_2','p_3','p_4','p_5','Location','best');

tout = flip(tout);
figure(5)
plot(tout,yp);
title('Measurement over Time (Adjoints)')
xlabel('Time (s)')
ylabel('Measurement (x_3)')

% part B
% Full LQR solution to system with NO noise
K0 = zeros(5);

[tkout,kout] = ode45(@(t,k) Kdot(k,A,B,R,Q),tbw,K0); % find history of K
xstar(:,1) = x0;
kout = reshape(kout.',5,5,[]);
tkout = flip(tkout);

for i = 1:length(tkout)
    k = kout(:,:,length(tkout)+1-i);
    ustar = -inv(R)*B'*k;
    xdotstar(:,i) = A*xstar(:,i) + B*ustar*xstar(:,i);
    xstar(:,i+1) = xstar(:,i) + xdotstar(:,i)*dt;
    ystar(i) = C * xstar(:,i);
    % calculate cost ... 
end

figure(6)
plot(tkout,ystar)
title('Measurement over Time (Optimal LQR)')
xlabel('Time (s)')
ylabel('Measurement (x_3)')

figure(7)
plot(tkout,xstar(:,1:end-1))
title('Optimal LQR Controller States a=1, b=1, \gamma=100')
xlabel('Time (s)');
ylabel('States');
legend('x1','x2','x3','x4','x5');
 

function dk = Kdot(K,A,B,R,Q)
    k = reshape(K,size(A));
    kdot = -k*A + k*B*inv(R)*B'*k - Q - A'*k;
    dk = kdot(:);
end
