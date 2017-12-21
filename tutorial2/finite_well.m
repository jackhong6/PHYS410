% Function to find the finite square well's bound states' energies.
% We define this script as a function so that the parameters of
% interest (width, mass, depth) are its input arguments.
% This function can be called from the command window.

function [E_even1,E_even2,E_odd1,E_odd2] = finite_well(a,m,v0)
% Save as finite_well.m in MATLAB.

hbar2 = calculate_hbar2;      % in units of eV*m_e*nm^2

% Defining functions
alpha = @(E) sqrt(2*m*E/hbar2);        % called as alpha(E)
beta = @(E) sqrt(2*m*(v0-E)/hbar2);    % called as beta(E)

f_even = @(E) beta(E).*cos(alpha(E).*a) - alpha(E).*sin(alpha(E).*a);
f_odd = @(E) alpha(E).*cos(alpha(E).*a) + beta(E).*sin(alpha(E).*a);
g_even = @(E) alpha(E).*tan(alpha(E)*a) - beta(E);
g_odd = @(E) alpha(E).*cot(alpha(E)*a) + beta(E);

% Plot function to approximately find roots
E = linspace(0,v0,100);
figure;
subplot(2,1,1);
plot(E,f_even(E),E,g_even(E),E,zeros(1,100));
ylim([-20,20]);
legend('f_{even}', 'g_{even}');

subplot(2,1,2);
plot(E,f_odd(E),E,g_odd(E),E,zeros(1,100));
ylim([-20,20]);
legend('f_{odd}', 'g_{odd}');

% Newton's method for root-finding
alpha_prime = @(E) sqrt(m/(2*hbar2*E));
beta_prime = @(E) -sqrt(m/(2*hbar2*(v0-E)));
f_even_prime = @(E) beta_prime(E).*cos(alpha(E).*a) - a.*beta(E).*sin(alpha(E).*a).*alpha_prime(E) ...
                    - alpha_prime(E).*sin(alpha(E).*a)-a.*alpha(E).*cos(alpha(E).*a).*alpha_prime(E);
f_odd_prime = @(E) alpha_prime(E).*(cos(a.*alpha(E)) - a.*alpha(E).*sin(a.*alpha(E))) ...
                   + beta_prime(E).*sin(alpha(E).*a) + a.*beta(E).*cos(alpha(E).*a).*alpha_prime(E);
                
E_even1 = newton_root(f_even, f_even_prime, 0.8);
E_even2 = newton_root(f_even, f_even_prime, 6.1);
E_odd1 = newton_root(f_odd,f_odd_prime,2.9);
E_odd2 = newton_root(f_odd,f_odd_prime,9.8);

end