function [u, err] = GenericHyperbolicErr(h, dt, Tf, alpha, beta, ...
               Uinflow, funU0, funUex, a, vN)
% --------------------------------------------------
%         Solution of the hyperbolic equation
%          with general numerical viscosity
% --------------------------------------------------
% Sintax:
%   [u, err] = GenericHyperbolicErr(h, dt, Tf, alpha, beta,
%         Uinflow, funU0, funUex, a, vN)
%
% Input:
%   dt          time discretization step
%   h           space discretization step
%   Tf          final time
%   alpha       left bound interval
%   beta        right bound interval
%   Uinflow     inflow condition
%   funU0       initial condition
%   funUex      exact solution
%   a           parameter of the equation representing
%               the velocity of the front
%   vN          numerical viscosity
%
% Output:
%   u           numerical solution
%   err         error

% Initialization
lambda = dt/h;
% disp(sprintf('CFL: %f', a*lambda));

Nh = floor((beta-alpha)/h)+1;
x = linspace(alpha,beta,Nh);

Nt = floor(Tf/dt)+1;

% Initial condition
u0 = feval(funU0, x);
uold = u0;
u = zeros(size(u0));

err = 0;

% Temporal loop
for n = 1:Nt-1;
    t = n*dt;
    u(2:Nh-1) = uold(2:Nh-1) - a * lambda / 2 * ( uold(3:Nh) - uold(1:Nh-2) ) ...
              + vN * dt / h^2 * ( uold(3:Nh) - 2* uold(2:Nh-1) + uold(1:Nh-2) );

    % Boundary condition
    u(1) = Uinflow;                                      % Inflow
    u(Nh) = (1-lambda*a)*uold(Nh) + lambda*a*uold(Nh-1); % Outflow

    uold = u;

    uex = feval(funUex,x,t,a);    
    err = max(err, max(abs(u - uex)));
    
%     plot(x, u, '-o', x, uex,'linewidth',2);
%     ylim([-0.5,1.5])
%     legend('Numerical solution','Exact solution', 'location','southwest');
% 
%     pause(dt);
end
