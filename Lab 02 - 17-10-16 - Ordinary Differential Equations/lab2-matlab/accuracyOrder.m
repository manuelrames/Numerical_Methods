function p = accuracyOrder(M, t0, tf, u0, dt, lambda)
% Estimation of the order of accuracy p
% of the forward Euler method
% 
%  Sintax:
%   p = orderAccuracy(M, t0, tf, u0, dt, lambda)
%
%  Input:
%       M       number of halved time steps
%       t0      initial time
%       Tf      final time
%       u0      initial condition
%       dt      time step discretization
%       lambda  coefficient of the reference problem
%
%  Output:
%       p       estimated order of accuracy
%

e = zeros(M,1);

for j = 1:M
    [u, e(j)] = FE(t0, tf, u0, dt, lambda);
    dt = dt/2;
end

p = log2( e(1:end-1) ./ e(2:end) )'
