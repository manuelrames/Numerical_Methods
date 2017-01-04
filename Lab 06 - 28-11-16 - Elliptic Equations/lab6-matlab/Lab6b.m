%Lab6b

h = 0.025;
L = 1;
u0 = 0;
uL = 0;
fun = 'x.^0';
Uex = 'x.*(1-x)/2';

Elliptic(L,h,u0,uL,fun,Uex);