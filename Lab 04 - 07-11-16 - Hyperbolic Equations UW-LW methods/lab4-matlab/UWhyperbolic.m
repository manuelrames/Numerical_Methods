function [u] = UWhyperbolic(h, dt, Tf, alpha, beta, Uinflow, funU0, funUex, a)
% --------------------------------------------------
%         Solution of the hyperbolic equation 
%             with the Upwind method (UW)
% --------------------------------------------------
% Sintax:
%   [u] = UWhyperbolic(h, dt, Tf, alpha, beta, Uinflow, funU0, funUex, a)
%
% Input:
%	dt          time discretization step
%	h           space discretization step
%   Tf          final time
%   alpha       left bound interval
%   beta        right bound interval
%   Uinflow     inflow condition
%   funU0       initial condition
%   funUex      exact solution
%   a           parameter of the equation representing
%               the velocity of the front
%
% Output:
%   u           numerical solution

% Initialization
lambda = dt/h;
disp(sprintf('CFL:  %f', a*lambda));

Nh = floor((beta-alpha)/h)+1;
x = linspace(alpha,beta,Nh);

Nt = floor(Tf/dt)+1;

% Initial condition
u0 = feval(funU0, x);
uold = u0;
u = zeros(size(u0));

% Temporal loop 
for n = 1:Nt-1;
    t = n*dt;
    u(2:Nh) = uold(2:Nh) - a * lambda * ( uold(2:Nh) - uold(1:Nh-1) );
    
    % Boundary condition
    u(1) = Uinflow;     % Inflow
    
    uold = u;
    
    % Plot
    uex = feval(funUex,x,t,a);
    plot(x, u, '-o', x, uex,'linewidth',2);
    ylim([-0.5,1.5])
    legend('Numerical solution','Exact solution', 'location','southwest');
    disp(sprintf('t = %f \n', t));
    pause(dt);
end