% Input parameters for the function FEgeneric
t0 = 0;
tf = 10;
u0 = 0;
fun = inline('1./(1+t.^2) - 2*y.^2','t','y');

% Solve the problem with dt = 0.95
dt = 0.95; 
[t1, u1] = FEGeneric(t0, tf, u0, dt, fun);

% Solve the problem with dt = 0.9
dt = 0.9; 
[t2, u2] = FEGeneric(t0, tf, u0, dt, fun);

% Solve the problem with dt = 0.2
dt = 0.2; 
[t3, u3] = FEGeneric(t0, tf, u0, dt, fun);

% Exact solution
t = [t0:0.01:tf]; 
y = t./(1+t.^2);

% Plot solutions
plot(t, y, 'k-');
hold on
plot(t1, u1, 'go-');
plot(t2, u2, 'ro-');
plot(t3, u3, 'bo-');
legend('y_{ex}','dt=0.95','dt=0.9','dt=0.2')