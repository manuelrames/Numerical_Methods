function [t, u] = FEgeneric(t0, tf, u0, dt, fun)
% ----  Forward Euler scheme for a generic Cauchy problem  ----
% Sintax:
%   [t, u] = FEgeneric(t0, Tf, u0, dt, lambda)
%
% Input:
%	t0      initial time
%	Tf      final time
%	u0      initial condition
%	dt      time step
%   fun     function f(t,y) of the differential problem
%
% Output:
%   t       time steps
%   u       numerical solution 
%           of the reference problem

% Definition of the temporal nodes
Nt = floor((tf-t0)/dt)+1;
t = [t0:dt:tf];

% Initialization of the solution vector u
u = zeros(1,Nt);
u(1) = u0;

% Temporal loop
for i = 1:Nt-1
   u(i+1) = u(i)+dt*feval(fun,t(i),u(i));
end