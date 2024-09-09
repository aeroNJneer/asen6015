clear
syms t
k=398600;
dt=10;
r=6678;
w=sqrt(k/(r^3));

x0=[r+10;-.5;0;w+2.5e-5]; % perturbation initial conditions

tspan = 0:10:6000; % sample every 10 seconds
opts = odeset('RelTol',1e-32,'AbsTol',1e-32);
[t,y] = ode45(@(t,y) myfunc(t,y,k), tspan, x0, opts);

figure;
s=stackedplot(t, y);
xlabel('Time t (s)');
s.DisplayLabels={'Radial Position (km)', 'Radial Velocity (km/sec)', 'Angular Position (rad)', 'Angular Velocity (rad/sec)'};
ax = findobj(s.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
title('Nonlinear State Change vs. Time (ode45)');
grid on;

function dy=myfunc(t,y,a)
    dy=zeros(4,1);

    dy(1)=y(2);
    dy(2)=y(1).*(y(4)^2)-(a./(y(1)^2));
    dy(3)=y(4);
    dy(4)=(-2.*y(2).*(y(4))./y(1));
end

