% PHYS 410: Homework 4
% Jack Hong, 30935134
% November 27, 2016

clear all; close all;
%% Define constants and functions
hbar2 = 0.076199682;  % eV * m_e * nm^2
m = 1;  % mass of electron
alpha = 500;  % eV nm^-2
beta = 3500;  % eV nm^-4 

x0 = 0.6;  % nm
E = 1;  % eV

epsilon = 1e-5;

V = @(x) -alpha*x.^2 + beta*x.^4 + alpha^2/(4*beta);
d2psi_dx2 = @(x, psi) 2*m*(V(x)-E)/hbar2 .* psi;

%% Solve ODE with RK4.m
f = @(x,y) [y(2), d2psi_dx2(x, y(1))];
xspan = [-x0, x0];
y0 = [0, epsilon];
[ylist, xlist] = RK4(f, xspan, y0);

%% Plot V(x)
figure(1)
x = linspace(-x0, x0);
plot(x, V(x))
title('Potential Well: $V(x)=-\alpha x^2+\beta x^4+\alpha^2/(4\beta)$', 'Interpreter', 'latex')
ylabel('V(x)')
xlabel('x')

%% Plot psi(x)
figure(2)
plot(xlist, ylist(:, 1));
title('Numerical Solution for $\psi(x)$', 'Interpreter', 'latex')
ylabel('$\psi(x)$', 'Interpreter', 'latex')
xlabel('x')

%% Even solution energy eigenvalues
n = 3;
even_energies = zeros(1, n);
even_guesses = [5, 15, 22];
for ind = 1:n
    even_energies(ind) = fzero(@even_endpoint, even_guesses(ind));
end
even_energies

%% Odd solution energy eigenvalues
n = 3;
odd_energies = zeros(1, n);
odd_guesses = [5, 16, 27];
for ind = 1:n
    odd_energies(ind) = fzero(@odd_endpoint, odd_guesses(ind));
end
odd_energies

%% Plot the six eigenstates
figure(3)
x = linspace(-x0, x0);
plot(x, V(x), 'DisplayName', 'Potential')
xlabel('x')
ylabel('$\psi(x)$', 'Interpreter', 'latex')
hold on

figure(4)
plot(x, V(x), 'DisplayName', 'Potential')
xlabel('x')
ylabel('$|\psi(x)|^2$', 'Interpreter', 'latex')
hold on

energies = [even_energies odd_energies];

for E = energies
    d2psi_dx2 = @(x, psi) 2*m*(V(x)-E)/hbar2 .* psi;
    f = @(x,y) [y(2), d2psi_dx2(x, y(1))];
    xspan = [-x0, x0];
    y0 = [0, epsilon];
    [ylist, xlist] = RK4(f, xspan, y0);
    a = trapz(xlist, abs(ylist(:, 1).^2));
    
    figure(3)
    plot(xlist, ylist(:, 1)./sqrt(a)+E, 'DisplayName', ['E=' num2str(E)])
    
    figure(4)
    plot(xlist, ylist(:, 1).^2./a+E, 'DisplayName', ['E=' num2str(E)])
end

figure(3); ylim([0, 30]); legend(gca, 'show');
figure(4); ylim([0, 35]); legend(gca, 'show');

