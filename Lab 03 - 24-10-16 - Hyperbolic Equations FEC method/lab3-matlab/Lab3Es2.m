%lab3es2
alpha = -3;
beta = 3;
Tf = 4.5;    
a = 0.5;

fun1 = inline('cos(pi*x).^2 .* (abs(x) <= 0.5)', 'x');
conditionInflow = 0;
funExact = inline('cos(pi*(x-a*t)).^2 .* (abs(x - a*t ) <= 0.5)', 'x','t','a');

h = 0.1;
dt = 0.1;
u = FEhyperbolic(h, dt, Tf, alpha, beta, conditionInflow, fun1, funExact, a);