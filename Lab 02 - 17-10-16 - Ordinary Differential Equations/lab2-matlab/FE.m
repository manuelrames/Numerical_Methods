function [u, err] = FE(t0, tf, u0, dt, lambda)
% ----  Solution of the reference problem  ----
%  y'(t) = lambda y(t), with y(t0) = u0
%  by using the forward Euler method
% ----   ----
% Sintax:
%   [u, err] = FE(t0, Tf, u0, dt, lambda)
%
% Input:
%	t0      initial time
%	Tf      final time
%	u0      initial condition
%	dt      time step
%   lambda  coefficient of the reference problem
%
% Output:
%   u       numerical solution 
%           of the reference problem
%   err     maximum absolute error
%           with respect to the exact solution


% Definition of the temporal nodes
Nt = floor((tf-t0)/dt)+1;
t = [t0:dt:tf];

% Initialization of the solution vector u
u = zeros(1,Nt);
u(1) = u0;

% Temporal loop
for i = 1:Nt-1
   u(i+1) = (1+lambda*dt)*u(i); 
end

% Computation of the error
uex = inline('u_0*exp(lambda*(t-t_0))','t','u_0','t_0','lambda');
y = feval(uex, t, u0, t0, lambda);
err = max( abs( y-u ) );

% Auxiliar vector for the representation of the exact solution
tt = [t0:0.01:tf];
% Exact solution
y_uex = feval(uex, tt, u0, t0, lambda);

plot(t, u, 'ro-', tt, y_uex);
legend('u', 'y_{ex}')