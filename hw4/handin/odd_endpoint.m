function [ psi0 ] = odd_endpoint( E )
% Shooting method

hbar2 = 0.076199682;  % eV * m_e * nm^2
m = 1;  % mass of electron
alpha = 500;  % eV nm^-2
beta = 3500;  % eV nm^-4 

x0 = 0.6;  % nm
epsilon = 1e-5;

V = @(x) -alpha*x.^2 + beta*x.^4 + alpha^2/(4*beta);
d2psi_dx2 = @(x, psi) 2*m*(V(x)-E)/hbar2 .* psi;

f = @(x,y) [y(2), d2psi_dx2(x, y(1))];
xspan = [-x0, 0];
y0 = [0, epsilon];
ylist = RK4(f, xspan, y0);

psi0 = ylist(end, 1);
end