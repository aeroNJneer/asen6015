

A = [0 1; -2 3];
B = [0; 1];
C = [1 0];
D = [0];

K = [-1 -1];

sys = ss(A, B, C, D);


tspan = [0 10];
x0 = [1; -1]; % Initial conditions
[t, x] = ode45(@sys, tspan, x0);




%% functions
function dy=myfunc(t,y)
init;
a = Gearth;
    dy=zeros(4,1);
    dy(1)=y(2);
    dy(2)=y(1).*(y(4)^2)-(a./(y(1)^2));
    dy(3)=y(4);
    dy(4)=(-2.*y(2).*(y(4))./y(1))
end


function dx = sys(t, x)

    A =;
    B = ;
    C = ;
    K = ;

    u=K*x
    dx = A*x + B*u;


end


