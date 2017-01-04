function [u, err] = AdvectionDiffusionCenteredStabilized(alpha, beta, h, ...
    uAlpha, uBeta, fun, funUex, mu, a, phi)
% ----  Solution of the advection-diffusion equation ----
%   -mu u'' + a u' = fun with Dirichlet boundary conditions in
% the [alpha, beta] interval with second order centered approximation
% with a generic stabilization
% ----   ----
% Sintax:
%   [u, err] = AdvectionDiffusionCenteredStabilized(alpha, beta, h,
%           uAlpha, uBeta, fun, funUex, mu, a, phi)
%
% Input:
%   alpha   left bound of the interval
%   beta    right bound of the interval
%   h       space step
%   uAlpha  boundary Dirichlet condition in x = alpha
%   uBeta   boundary Dirichlet condition in x = beta
%   fun     function describing the term f in the equation
%   funUex  exact solution
%   mu      diffusion coefficient
%   a       advection coefficient
%   phi     phi(Pe) function used to apply the stabilization
%
% Output:
%   u numerical solution
%     of the advection-diffusion problem
%   err error with respect to the exact
%       solution in the infinity norm

Pe = a * h / (2*mu); % Peclet number
disp(sprintf('Peclet = %f', Pe));

muh = mu * ( 1 + feval(phi, Pe) );

% Initialization
M = floor((beta - alpha)/h); % # of intervals
x = linspace(alpha, beta, M+1)'; % M+1 nodes

u = zeros(M-1, 1); % only internal nodes
f = eval(fun); % all nodes

% Construction of the matrix in sparse format
e = ones(M-1, 1);
Amu = 1/h^2 * spdiags([-e 2*e -e], [-1 0 1], M-1, M-1);
Aa = 1/(2*h) * spdiags([-e e], [-1 1], M-1, M-1);
A = muh*Amu + a*Aa;

% Right hand side
F = f(2:end-1); % only internal nodes

% Correction of the known term with the boundary conditions
F(1) = F(1)+ (muh/h^2 + a/(2*h)) * uAlpha;
F(end) = F(end) + (muh/h^2 - a/(2*h)) * uBeta;

% Compute the solution in the internal nodes
u = A\F;

% Append the Dirichlet boundary conditions
u = [uAlpha; u; uBeta];

% Plot
xx = [alpha:0.01:beta]'; % Just for visualization
uex = feval(funUex, xx, mu, a); % Evaluated in 'xx'
plot(x, u,'o-', xx, uex, 'linewidth', 2);
legend('Numerical solution', 'Exact Solution', 'location', 'northwest')

% Error
uex = feval(funUex, x, mu, a); % Evaluated in 'x'
err = max( abs(u-uex) );