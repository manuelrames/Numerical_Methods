%lab4es3
alpha = -3;
beta = 3;
Tf = 2;
a = 1;

fun1 = inline('cos(pi*x).^4 .* (abs(x) <= 0.5)', 'x');
conditionInflow = 0;
funExact = inline('cos(pi*(x-a*t)).^4 .* (abs(x - a*t) <= 0.5)','x','t','a');

h = 0.1;
dt = 0.05;

u = UWhyperbolic(h, dt, Tf, alpha, beta, conditionInflow, fun1, funExact, a);

%%

% UW method
vUW = abs(a) * h / 2; % numerical viscosity for UW
u = GenericHyperbolic(h, dt, Tf, alpha, beta, conditionInflow, fun1, funExact, a, vUW);

%%

% LW method
vLW = a^2 * dt / 2;   % numerical viscosity for LW
u = GenericHyperbolic(h, dt, Tf, alpha, beta, conditionInflow, fun1, funExact, a, vLW);