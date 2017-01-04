% Lab9Es1b

L = 0.05;
T = 5.e-5;
M = 40; % number of space intervals
h = 2*L/M;
t0 = h^2;
alpha = 3/4;
dt = alpha*h^2;
funU0 = inline('1.0/sqrt(4.0*pi*t)*exp(-x.^2/(4.0*t))', 'x', 't');

FEheatEquation(L, t0, T, h, dt, funU0);