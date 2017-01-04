% Lab9Es1c

L = 0.05;
T = 5.e-5;
M = 40; % number of space intervals
h = 2*L/M;
t0 = h^2;
alpha = 1/2;
dt = alpha*h^2;
funU0 = inline('1.0/sqrt(4.0*pi*t)*exp(-x.^2/(4.0*t))', 'x', 't');

% Divide by 2 both h and dt: unstable solution
h = h/2;
dt = dt/2;

% Stable solution
%h = h/2;
%dt = dt/4;

disp(sprintf('Is  dt <= 1/2 * h^2  fulfilled? %d', dt <= h^2/2));
disp(sprintf('dt = %5.2e', dt));
disp(sprintf('1/2 * h^2 = %5.2e', h^2/2));

FEheatEquation(L, t0, T, h, dt, funU0);