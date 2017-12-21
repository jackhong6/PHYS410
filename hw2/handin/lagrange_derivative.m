% PHYS 410: Homework 2, Question 1
% Jack Hong, 30935134

function [ ] = lagrange_derivative()
% Return the derivative of f at x_star using a 7-points centered formula.
%    - Required weights are derived using the idea of differentiating the Lagrange
%      Interpolation of f using the 7 points around x_star.
%    - value of h optimized to minimize the combination of discretization
%      and round off errors.

TE = @(h) 2181.29/140 .* h.^6 + 11/6.*eps./h;  % Total error upper bound

dTE_dh =  @(h) 2181.29*6*h^5/140 - 11*eps/6/h^2; % derivative of total error with respect to h
h_optimal = fzero(dTE_dh, 0.001); 

f = @(x) sin(x.^2);
f_prime = @(x) 2.*x.*cos(x.^2);

lagrange = @(x0,h) (3/4*(f(x0+h)-f(x0-h))./h) - (3/20*(f(x0+2*h)-f(x0-2*h))./h)...
                 + (1/60*(f(x0+3*h)-f(x0-3*h))./h);

% ============================ Plot derivative ============================
x = linspace(0,1, 100);
l_pnts = lagrange(x,h_optimal);
f_prime_pnts = f_prime(x);

figure(1);
subplot(2,1,1);
plot(x, l_pnts, x, f_prime_pnts); % only see one line since approximation is good.
title('Derivative of f(x) = sin(x^2)');
xlabel('x');
ylabel('f''(x)');
legend('Lagrange', 'f''(x) = 2xcos(x^2)');

subplot(2,1,2);
plot(x, l_pnts - f_prime_pnts);
title('Residuals between Lagrange approximation and f''(x)')
hold off;
% =========================================================================

% ===================== Plot errors as function of h ======================
figure(2)
num_pnts = 1000;
hs = linspace(1e-4, 1e-2, num_pnts);
x_sample = 0.5; % arbitrary point to compare errors.
errors = abs(lagrange(x_sample,hs) - f_prime(x_sample));

loglog(hs, TE(hs)); hold on;
scatter(hs, errors);
title('Errors as a function of h')
legend('Theoretical error upper bound', 'Error');
xlabel('h');
ylabel('Errors');


end

