function E_bound = single_well()
% PHYS 410: Project 1, Problem 1; Solution to exercise 2.16
% 
% Function to find the bound state energies of the single well potential using
% the propagator method discussed in section 2.13 of the textbook. User
% changeable parameters are placed at the top of the code.
%
% Author: Jack Hong, 30935134
% Last Modified: September 30, 2016 

% ---------------------- Define user parameters: --------------------------
b = 1;          % Left limit of the square well in nm
c = 1.6;        % Right limit of the square well in nm
v0 = 10;        % The depth of the well in eV
m = 1;          % The mass of the particle in m_e (mass of electron)

plotfig = 1;    % If plotfig == 1, the script will generate a plot of the function.
                % Otherwise, it won't.
        
tol = 1e-8;     % Tolerance for when a root is "found"
max_iter = 100; % Maximum number of iterations for running the Newton method.

% ----------------Define constants and coefficients: ----------------------
hbar2 = calculate_hbar2;   % approx.  0.0762 eV*m_e*nm^2
alpha = @(E) sqrt(2*m*E/hbar2);
beta = @(E) sqrt(2*m*(v0-E)/hbar2);
w = c - b;
aw = @(E) alpha(E).*w;  % used to simplify expressions

% Derivatives of alpha and beta with respect to energy. (Found by hand)
alpha_prime = @(E) sqrt(m/(2*hbar2*E));
beta_prime = @(E) -sqrt(m/(2*hbar2*(v0-E)));

% Define the function that we want to find the root of
% f is obtained by arbitrarily setting psi(b)=1 and psi'(b)=beta and using
% the P_allowed matrix to propagate the solution to c.
f = @(E) -alpha(E).*sin(aw(E)) + 2.*beta(E).*cos(aw(E)) ...
         + (beta(E).^2)./alpha(E).*sin(aw(E));
% The derivative of f with respect to energy. (Found by hand)
f_prime = @(E) -alpha_prime(E).*sin(aw(E)) - w.*alpha_prime(E).*alpha(E).*cos(aw(E)) ...
               + 2.*beta_prime(E).*cos(aw(E)) - 2.*w.*alpha_prime(E).*beta(E).*sin(aw(E))...
               + (2.*beta_prime(E).*alpha(E)-beta(E).^2.*alpha_prime(E))/(alpha(E).^2).*sin(aw(E))...
               + (beta(E).^2/alpha(E)).*w.*alpha_prime(E).*cos(aw(E));

% Find the roots using the Newton method.

% Initialize array to store the energies found. Indices in E_bound
% correspond to indices in init_guesses
init_guesses = get_guesses(f, 0, v0, 1,4);

E_bound = zeros(1,length(init_guesses));
index = 1; % Index of E_bound to store the energy found
for E_init = init_guesses
    [root, num_iter] = newton_root(f, f_prime, E_init, tol, max_iter);
    if num_iter < max_iter % If the root was actually found by newton_root,
        E_bound(index) = root; % assign the root found to the slot in E_bound
    else % otherwise, if the root was not found,
        E_bound(index) = NaN; % assign NaN to the spot in E_bound
    end
    index = index + 1;
end

% Plot used to estimate the roots of the function
if plotfig == 1
    E = linspace(0,v0,1000);
    
    figure(1);
    subplot(2,1,1);
    plot(E, f(E));
    line([0 v0],[0 0],'Color',[0 0 0]);
    title({'f(E) = \psi''(c) + \beta\psi(c)';...
        'Roots are bound state energies for a single square well';...
        'v0 = 10eV; width = 0.6nm'});
    ylabel('\bf{f(E)}');
    xlabel('\bf{E}');
    grid minor;
    
    subplot(2,1,2);
    plot(E, f(E));
    ylim([-3,3]);
    line([0 v0],[0 0],'Color',[0 0 0]);
    title({'ZOOMED IN OF ABOVE: f(E) = \psi''(c) + \beta\psi(c)';...
        'Roots are bound state energies for a single square well';...
        'v0 = 10eV; width = 0.6nm'});
    ylabel('\bf{f(E)}');
    xlabel('\bf{E}');
    grid minor;
end

end