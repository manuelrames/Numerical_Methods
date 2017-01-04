function [u] = FEhyperbolic(h, dt, Tf, alpha, beta, Uinflow, funU0, funUex, a)
% --------------------------------------------------
%         Solution of the hyperbolic equation 
%    with the forward Euler/centered method (FE/C)
% --------------------------------------------------
% Sintax:
%   [u] = FEhyperbolic(h, dt, Tf, alpha, beta, Uinflow, funU0, funUex, a)
%
% Input:
%	h          space discretization step
%	dt         time discretization step
%   Tf         final time
%   alpha 	   left bound interval
%  	beta 	   right bound interval
%   Uinflow    inflow condition
%   funU0 	   initial condition
%   funUex 	   exact solution
%   a          parameter of the equation representing
%     		   the velocity of the front
%
% Output:
%   u 	   numerical solution

% Initialization
lambda = dt/h;

Nh = floor((beta-alpha)/h)+1;
x = linspace(alpha,beta,Nh);

Nt = floor(Tf/dt)+1;

% Initial condition
u0 = feval(funU0,x);
uold = u0;
u = zeros(size(u0));

% Temporal loop 
for n = 1:Nt-1;
    t = n*dt;
    u(2:(Nh-1)) = uold(2:Nh-1) - a*lambda/2 * ( uold(3:Nh) - uold(1:Nh-2) );
    
    % Boundary condition
    u(1) = Uinflow;                                      % Inflow
    u(Nh) = (1-lambda*a)*uold(Nh) + lambda*a*uold(Nh-1); % Outflow
    
    uold = u;
    
    % Plot
    uex = feval(funUex,x,t,a);
    plot(x, u, '-o', x, uex,'linewidth',2);
    ylim([-5,5])
    legend('Numerical solution','Exact solution', 'location','southwest');
    disp(sprintf('t = %f \n', t));
    pause(dt);
end
