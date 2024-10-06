clear

syms g; syms mu; syms r;syms x; syms y; syms z; syms R;

R = [x;y;z]
r = norm(R)
r_hat = R/r
g = (mu/r^3)*R
dg = jacobian(g,R)
ddg = jacobian(dg,R)