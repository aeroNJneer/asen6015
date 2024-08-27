% quiz 4

F = [1,1,1;1,1,1;1,1,1];
H = [1,2,3];

O = [H;H*F;H*F^2];
eig(O)
[V,d]=eig(O)
