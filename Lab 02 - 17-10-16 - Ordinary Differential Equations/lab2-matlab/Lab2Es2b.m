% Input parameters for the function BE
lambda = -2;
t0 = 0;
tf = 12;
u0 = 1;

% Solve the reference problem with dt = 0.05
dt = 0.05;
figure(1)
[u, err] = BE(t0, tf, u0, dt, lambda);
title('BE method - dt = 0.05')

%Solve the problem with dt = 1.2
dt = 1.2;
figure(2)
[u, err] = BE(t0, tf, u0, dt, lambda);
title('BE method - dt = 1.2')