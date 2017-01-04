function [u] = BEhyperbolic(h, dt, Tf, alpha, beta, Uinflow, funU0, funUex, a)
% ----------------------------------------------
%      Solution of the hyperbolic equation 
%    with the backward Euler/centered method
% ----------------------------------------------
% Sintax:
%   [u] = BEhyperbolic(h, dt, Tf, alpha, beta, Uinflow, funU0, funUex, a)
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
%
% Output:
%   u           numerical solution

lambda = dt/h;
disp(sprintf('CFL:= %f \n', a*lambda));

Nh = floor((beta-alpha)/h)+1; 
x = linspace(alpha,beta,Nh);

Nt = floor(Tf/dt)+1;

% Initial condition
u0 = feval(funU0, x);
uold = u0;
u = zeros(size(u0));
 
e = ones(Nh-1,1);
A = speye(Nh-1) + a*lambda/2 * spdiags([-e e], [-1 1], Nh-1, Nh-1);

% Implicit UW for the last node
A(Nh-1,Nh-1) = 1 + a*lambda; %A(Nh-1,Nh-1) = A(Nh-1,Nh-1) + a*lambda;
A(Nh-1,Nh-2) = -a*lambda; %A(Nh-1,Nh-2) = 2*A(Nh-1,Nh-2);
 
F = zeros(Nh-1,1);
     
for n = 1:Nt-1
    t = n*dt;
          
    F(2:Nh-1) = uold(3:Nh);
     
    % First row of the system has the rhs modified by the inflow condition
    F(1) = uold(2)+a*lambda/2*Uinflow;    
     
    u(2:Nh) = A\F;
    
    u(1) = Uinflow; 
       
    uold = u;
     
    uex = feval(funUex,x,t,a);
    plot(x, u, '-o', x, uex,'linewidth',2);
    legend('Numerical solution BE/C','Exact solution');	
    pause(dt);
end     
