function FEheatEquation(L, t0, T, h, dt, funU0)
% ----  Solution of the heat equation
%   with homogeneous Dirichlet boundary conditions 
%   with the forward Euler method in the interval [-L,L]
% ----   ----
% Sintax:
%   FEheatEquation(L, t0, T, h, dt, funU0)
%
% Input:
%   L length semi-interval
%   t0 inizial time 
%   T final time
%   h  space step
%   dt time step
%   funU0 initial condition

M = floor(2*L/h);            % number of space intervals
N = floor((T-t0)/dt);        % number of time intervals
x = linspace(-L, L, M+1);    % spatial nodes
tt = linspace(t0, T, N+1);   % temporal nodes

e = ones(M-1, 1);
A = h^(-2)*spdiags([-e 2*e -e], [-1 0 1], M-1, M-1);

% Initialization
u = zeros(M+1, N+1);    % total solution
Uint = zeros(M-1,N+1);  % solution in the internal nodes

% Initial condition
u(:, 1) = feval(funU0, x, t0); 
Uint(:, 1) = u(2:end-1, 1);

% Boundary conditions
u(1, :) = zeros(N+1, 1);
u(end, :) = zeros(N+1, 1);

% FE/C method applied to the internal nodes
for n = 1:N
    Uint(:, n+1) = (eye(M-1) - dt*A) * Uint(:, n);
end

u(2:end-1,:) = Uint;

% Plot the numerical solution
figure(1);    
mesh(tt, x', u);  
title('Numerical solution')
xlabel('t');   
ylabel('x');  
axis([dt T -L L -20 100]);

% Exact solution (approximation)
uex = zeros(M+1, N+1);
for i = 2:M
    for n = 1:N+1
        uex(i,n) = 1/sqrt(4*pi*tt(n))*exp(-x(i)^2/(4*tt(n)));
    end
end
uex(1,:) = zeros(N+1,1);
uex(end,:) = zeros(N+1,1);

% Plot the exact solution
figure(2);    
mesh(tt, x', uex);
title('Exact solution')
xlabel('t');   
ylabel('x');  
axis([dt T -L L -20 100]);