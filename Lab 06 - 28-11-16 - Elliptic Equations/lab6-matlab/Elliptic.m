function [u, err] = Elliptic(L, h, u0, uL, fun, Uex)
% ----  Solution of the elliptic equation ----
%   -u'' = f with Dirichlet boundary conditions
%   in the interval (0,L)
% ----   ----
% Sintax:
%   [u, err] = Elliptic(L, h, u0, uL, fun, Uex)
%
% Input:
%   L upper bound of the interval
%   h  space discretization step
%   u0 Dirichlet condition in x=0
%   uL Dirichlet condition in  in x=L
%   fun function describing the term f in the equation
%   Uex exact solution
%
% Output:
%   u numerical solution of the elliptic problem
%   err error w.r.t. the exact solution

% Initialization
M = floor(L/h);            % # of intervals
x = linspace(0, L, M+1)';  % M+1 nodes

u = zeros(M-1, 1);         % Solution vector (only in the internal nodes)

% Construction of the matrix in sparse format
e = ones(M-1, 1);
A = h^(-2) * spdiags([-e, 2*e, -e], [-1 0 1], M-1, M-1);

% Right hand side
f = eval(fun);
F = f(2:end-1); % Known term (only in the internal nodes)

% Correction of the known term with the boundary conditions
F(1) = F(1) + u0/h^2;
F(end) = F(end) + uL/h^2;

% Compute the solution (only in the internal nodes)
u = A\F;

% Append the Dirichlet boundary conditions
u = [u0; u; uL];

% Plot
uex = eval(Uex);
plot(x, u, '-o', x, uex, 'linewidth',2);
legend('Numerical solution', 'Exact solution', 'Location', 'NorthWest');

% Compute the error
err = max(abs(u-uex));