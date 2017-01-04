%Lab5Es2
alpha = -3;
beta = 3;
Tf = 2;

a = 1;
h = 0.1;

CFL = 0.5;
dt = CFL*h/a;

vUW = abs(a)*h/2;
vLW = a^2*dt/ 2;

fun1 = inline('cos(pi*x).^4 .* (abs(x) < 0.5)');
inflowCondition = 0;
funExact = inline('cos(pi*(x-a*t)).^4 .* (abs(x - a*t) < 0.5)','x','t','a');

M = 6;
errUW = zeros(1, M);
errLW = zeros(1, M);

for k = 0:M
    [u errUW(k+1)] = GenericHyperbolicErr(h, dt, Tf, alpha, beta, inflowCondition, fun1, funExact, a, vUW);
    [u errLW(k+1)] = GenericHyperbolicErr(h, dt, Tf, alpha, beta, inflowCondition, fun1, funExact, a, vLW);
    h = h/2;
    dt = dt/2;
    vUW = abs(a)*h/2;
    vLW = a^2*dt/2;
end

pUW = log2(errUW(1:end-1)./errUW(2:end))
pLW = log2(errLW(1:end-1)./errLW(2:end))
