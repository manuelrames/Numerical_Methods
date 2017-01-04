% Evaluation of the order of accuracy of the FE method
t0 = 0;
tf = 12;
u0 = 1;
lambda = -2;

M = 8;
dt = 0.05;
p = accuracyOrder(M, t0, tf, u0, dt, lambda);