%Lab6c

h = 0.025;
L = 1;
u0 = 0;
uL = 0;

% case a)
fun = 'x.^0';
Uex = 'x.*(1-x)/2';

% case b)
% fun = '4*pi^2*cos(2*pi*x)';
% Uex = 'cos(2*pi*x)-1'

iter = 6;
e = zeros(iter,1);
H = [];

for j = 1:iter
    [U, e(j)] = Elliptic(L,h,u0,uL,fun,Uex);
    H = [H, h];
    h = h/2;
end

format long
e
format
p = log2( e(1:end-1) ./ e(2:end) )'

% Plot of the error as function of h (in logarithmic scale)
loglog(H, e, H, H, H, H.^2)
legend('err', 'H', 'H^2')
xlabel('h')
ylabel('err')
grid on