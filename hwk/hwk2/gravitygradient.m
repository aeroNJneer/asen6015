clear

syms g; syms mu; syms r;syms x; syms y; syms z; syms R;

R = [r 0 0]
g = mu/R * eye(3)
dg = jacobian(g,R)
ddg = jacobian(dg,R)
