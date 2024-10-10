% ASEN-6015 Hwk 3 

% a) gamma, a, b neq 0
% b) a, b = 0
%c) gamma neq 0, a or b new 0

clear
close all
a = [1; 0];
b = [1; 1];
gam = [0.5; 0.5; 0.5];
x0 = [1 -1 1 -1 1]';
dt = 0.1;
tsim = 1:dt:10;
u = zeros(2,length(tsim));
u(:,1) = 1;         % impulse function

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

Q = zeros(5,5);
i=1;
%for i = 1:length(a)
    Q(3,3) = gam(i) * x0(3)^2;

    R = [a(i) 0;
         0 b(i)];    
    [~,pd] = chol(R)

    [K,P,poles] = lqr(A,B,Q,R);
    H = A-B*K;
    eig(H)
    
    S = Q + K'*R*K;
    
    sys1 = ss(H,B,C,D);
    
    J = x0'*P*x0;       % minimize J by adjusting P/K
    
   [ycalc, tout, xcalc] = initial(sys1,x0,tsim)
   figure(i)
   plot(tout,xcalc)
   title('LQR Controller a=1, b=1, \gamma=0.5')
   xlabel('Time (s)');
   ylabel('States');
   legend('x1','x2','x3','x4','x5');
   figure(2)
   plot(tout,ycalc)
%end
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

L = eye(5);
Rwi = [0:1 0:2 0:1 0:3 0:2];
wt = rand(5)' * Rwi'