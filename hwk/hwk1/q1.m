init;
w = sym('w'); 
r = sym('r');
dw = sym('dw');
dr = sym('dr');
mu = sym('mu');

A = [0 1 0 0 0 0;
    dw^2+(2*mu/(r^3)) 1 0 2*r*dw 0 0;
    0 0 0 1 0 0;
    (2*dr*dw/(r^2))-dw/(r^2) (-dw/r)+(1/r) 0 (dr/r)+(1/r) 0 0;
    0 0 0 0 0 1;
    0 0 0 0 0 0];
B = [0 0 0; 1 0 0 ; 0 0 0; 0 1 0; 0 0 0; 0 0 0]

Co = [B A*B A^2*B A^3*B]

rank(Co)
