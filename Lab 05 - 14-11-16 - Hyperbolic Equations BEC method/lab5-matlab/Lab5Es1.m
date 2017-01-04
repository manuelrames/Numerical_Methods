%Lab5Es1  
alpha = -3;
beta = 3; 
Tf = 1;

a = 1;
h = 0.1;

CFL = 2.0; 
dt = CFL*h/a;

fun1 = inline('cos(pi*x).^4 .* (abs(x) < 0.5)');
inflowCondition = 0;
funExact = inline('cos(pi*(x-a*t)).^4 .* (abs(x - a*t) < 0.5)','x','t','a');
[u] = BEhyperbolic(h, dt, Tf, alpha, beta, inflowCondition, fun1, funExact, a);