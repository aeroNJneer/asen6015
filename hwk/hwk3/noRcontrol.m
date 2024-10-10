syms x1; syms x2; syms x3; syms x4; syms x5; syms q; syms X; syms Y;

X = [x1;x2;x3;x4;x5];
XT = [x1 x2 x3 x4 x5]
q = [0 0 0 0 0;
    0 0 0 0 0;
    0 0 0.5 0 0;
    0 0 0 0 0;
    0 0 0 0 0]

Y = XT*q*X

