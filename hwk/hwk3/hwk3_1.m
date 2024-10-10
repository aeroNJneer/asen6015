% ASEN-6015 Hwk 3 

% a) gamma, a, b neq 0
% b) a, b = 0
%c) gamma neq 0, a or b new 0

clear
close all
clc
a = [1; 0; 10];             % lower a results in lower J and better response
b = [0.1; 0; 0.1];          % lower b results in lower J, but lower response
q = [1; 0.1; 0.1];       % lower q results in lower J
x0 = [1 -1 1 -1 1]';
dt = 0.1;
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
sys = ss(A,B,C,[]);

i=1;
for i = 1:length(a)
    % update Q and R
    Q = [0 0 0 0 0;
        0 0 0 0 0;
        0 0 q(i) 0 0;
        0 0 0 0 0;
        0 0 0 0 0];

    R = [a(i) 0;
         0 b(i)];    
    % if R = 0 then cost does not depend on R
    % so minimize J = integral x'Q x dt
    % 
    
    [Pstar, K, L] = icare(A,B,Q,R,[],[],[])
    
    try
        [K,P,poles] = lqr(A,B,Q,R);
    catch
        K = zeros(2,5);
    end
    % Lyupanov equation     H'P + PH = -S
    H = A-B*K;
    eig(H)      % is closed loop system symptotically stable?
    
    S = Q + K'*R*K;
    eig(S)      % is S PSD or PD?
    sys1 = ss(H,B,C,D);
    
    J(i) = x0'*P*x0       % minimize J by adjusting P/K
    
   [ycalc, tout, xcalc] = initial(sys1,x0,tsim)
   figure(i)
   plot(tout,xcalc)
   title('LQR Controller a=1, b=1, \gamma=0.5')
   xlabel('Time (s)');
   ylabel('States');
   legend('x1','x2','x3','x4','x5');
   % figure(i+1)
   % plot(tout,ycalc)
end
min(J)
Stm = expm(A*20);

x = x0;

for i = 1:length(tsim)
    xdot(:,i) = H*x(:,i) + B*u(:,i);
    x(:,i+1) = Stm * x(:,i) + Stm*B*u(:,i)*dt;
    y(i) = C * x(:,i);
end

    % steps in lrq
    % solve ricatti equation for P*, K, L (eigenvalues)
    % [Pstar, K, L] = icare(A,B,Q,R,[],[],[]);
    % chol(Pstar)
    % K = inv(R)*B'*Pstar;
    % H = A-B*K;      % must be asymptotically stable
    % eig(H)          % all negative real parts in eigenvalues
    
% Problem 2
% 
% L = eye(5);
% Rwi = [0:1 0:2 0:1 0:3 0:2];
% wt = rand(5)' * Rwi'

