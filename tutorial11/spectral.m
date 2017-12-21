N = 10;
[D, x] = cheb(N);
x = x./2 + 1/2;
D = 2.*D;

L = diag(x(1-x))*D^2 - diag(3*x)*D - eye(N+1);
id = eye(N+1);
f = zeros(1,N+1); f(1) = 1/2; f(N+1) = 1;

