clear

syms g; syms mu; syms r;syms x; syms y; syms z; syms R;

R = [x y z]
g = mu/norm(R)^3 * R
dg = jacobian(g,R)
ddg = jacobian(dg,R)
