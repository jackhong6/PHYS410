function maxerror = finitedifference(N)

% Writing the system of equations to be solved
x = linspace(0,1,N+1)'; h = x(2) - x(1);
L = zeros(N+1,N+1);         % preallocate linear operator
for i=2:N
    L(i,i) = -(2*x(i)*(1-x(i))/h^2 + 1);
    L(i,i+1) = x(i)*(1-x(i))/h^2 - 3*x(i)/(2*h);
    L(i,i-1) = x(i)*(1-x(i))/h^2 + 3*x(i)/(2*h);
end
% Boundary conditions
L(1,1) = 1;
L(N+1,N+1) = 1;
f = zeros(N+1,1); f(1) = 1; f(N+1) = 1/2;

% Solving and comparing to actual solution
u = L\f;
usol = (1-x+x.*log(x))./(x-1).^2;
maxerror = max(abs(u-usol));

% Higher resolution for solution
xx = linspace(0,1,1000);
u0 = (1-xx+xx.*log(xx))./(xx-1).^2;

plot(x,u,'.-',xx,u0)



end