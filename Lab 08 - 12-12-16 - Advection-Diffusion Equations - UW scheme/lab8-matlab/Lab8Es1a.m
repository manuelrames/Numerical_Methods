%Lab8es1a

alpha = 0;
beta = 1;
h = 0.1;
mu = 0.01;
a = 1;
uAlpha = 0;
uBeta = 1;
fun = '0.*x';
funUex = inline('(exp(a*x/mu)-1)./(exp(a/mu)-1)','x','mu','a');

[u] = AdvectionDiffusionUW(alpha, beta, h, ...
    uAlpha, uBeta, fun, funUex, mu, a);