% Midterm #2
close all;
% Q3
clear
omega = 0.045;
dt = 0.5;

A = [0 1 0 0 ; 0 0 0 -omega; 0 0 0 1; 0 omega 0 0 ];

F = [1 (sin(omega*dt))/omega 0 -(1-cos(omega*dt))/omega; 
    0 cos(omega * dt) 0 -sin(omega*dt); 
    0 (1-cos(omega*dt))/omega 1 (sin(omega*dt))/omega;
    0 sin(omega*dt) 0 cos(omega*dt)];

HA = [1 0 0 0; 0 0 1 0];

%ymsmt = importdata("OneDrive - Microsoft/Documents/school/ASEN5044DynAnal/hw6problem5data.csv");
y3b = importdata("midterm2_problem3b.mat");

k=length(y3b);
k=1
RA = [75 7.5; 7.5 75] + [12.5*sin(k/10) 25.5*sin(k/10); 25.5*sin(k/10) 12.5*cos(k/10)];
PArr = []; xhatArr = []; 
sigma = [];

M = [0; 85*cos(pi/4); 0; -82*sin(pi/4)];


x0 = inv(HA'*inv(RA)*HA) * (HA'*inv(RA)*y3b(1,1));

H30 = []; y30=[]; R30 = cell(30);
for i= 1:30
    H30 = [H30;H];
    y30 = [y30;ymsmt(1:3,i)];
    for j=1:30
        if i==j
            R30{i} = [R30{i};R];
        else
            R30{i} = [R30{i};zeros(3)];
        end
    end
end


p0 = RA
P0 = diag(p0);
PArr(:,:,1) = P0;
xhatArr(:,:,1) = M;

H = eye(3);

for k = 1:300
   sigma(:,k) = 2 * sqrt(diag(P0));
   P0 =  F * P0 * F';
   PArr(:,:,k+1) = P0;
   M = F*M;
   xhatArr = [xhatArr M];
end
% 
% % Q4
% 
% f = 0.8;
% q = 10;
% sig = 10;
% p = sig;
% p50 = 50;
% for k = 1:50
%     p(:,k+1) = (f^(2)*p(:,k))+q;
%     p50(:,k+1) = (f^(2)*p50(:,k))+q;
%     
% end
% % figure 
% % hold on;
% %  plot((1:1:51),p)
% %  plot((1:1:51),p50)
% %  title('Covariance Prediction with Various \sigma_0')
% %  set(gca,'FontSize',12);
% %  legend('\sigma^2_0 = 10','\sigma^2_0 = 50')
% %  xlabel('Time Step \Deltat'); ylabel('Covariance (P)')
% %  xlim([1,50])
% % hold off
% 
% % part d
% F = [0.99 0.2;0 -0.76];
% Q = [1 0.37;0.37 2.5];
% n = 200;
% P = 10*eye(2);
% p=zeros(2)
% PqArr(:,:,1) = eye(2);
% PArr(:,:,1) = P;
% 
% for k=1:n
%     PArr(:,:,k+1) = F*PArr(:,:,k)*F' +Q';
%    % PqArr(:,:,k+1) =  (PqArr(:,:,k) + ((F')^k * Q * F^k)) ;
%    % P =  P + ((F')^k * Q * F^k);
% end
% 
% Pinf = dlyap(F,Q)
% 
% for i = 1:n
%     p1(i) = PArr(1,1,i);
%     p2(i) = PArr(1,2,i);
%     p3(i) = PArr(2,2,i);
% end
% 
% % figure(2); % open a new figure
% % set(gcf,'Color','white');
% % k = 1:1:n;
% % 
% % subplot(3,1,1); hold on;
% % plot(k, p1, 'r-', 'LineWidth', 2);
% % plot(k,Pinf(1:1,1:1)*ones(n,1),'b:', 'LineWidth',2)
% % title('Covariance Prediction (Iterative vs dlyap.m)')
% % ylabel('Covariance (P_1_1)')
% %  legend('Iterative','dlyap')
% % subplot(3,1,2); hold on;
% % plot(k,p2);
% % plot(k,Pinf(2:2,1:1)*ones(n,1),'b:', 'LineWidth',2)
% % subplot(3,1,3); hold on;
% % plot(k,p3,'r-', 'LineWidth', 2);
% % plot(k,Pinf(2:2,2:2)*ones(n,1),'b:', 'LineWidth',2)
% % ylabel('Covariance (P_1_2,P_2_1)')
% %  legend('Iterative','dlyap')
% % ylim([5,10])
% % set(gca,'FontSize',12);
% % legend('\sigma^2_0 = 10','\sigma^2_0 = 50')
% % xlabel('Time Step \Deltat'); ylabel('Covariance (P_2_2)')
% % xlim([1,n])
% %  legend('Iterative','dlyap')
% %  hold off
