function [ E_bound ] = n_square_wells(n)
% PHYS 410: Project 1, Problem 3; Solution to exercise 2.19
% 
% Function to find the bound state energies of the double well potential using
% the propagator method discussed in section 2.13 of the textbook. The
% double wells are assumed to be identical. All user parameters are placed
% at the top of the code.
%
% Author: Jack Hong, 30935134
% Last Modified: September 30, 2016 

% ---------------------- Define user parameters: --------------------------
w = 0.6;        % The width of the wells in nm
s = 0.1;        % The separation between the two wells in nm
v0 = 10;        % The depth of the wells in eV
m = 1;          % The mass of the particle in m_e (mass of electron)
if nargin == 1
    N = n;
else
    N = 5;          % The number of wells
end

plotfig = 1;    % If plotfig == 1, a plot of the f(E) will be generated.
plot_evidence_for_band_structure = 1;
                        
tol = 1e-8;     % Tolerance for when a root is "found"
max_iter = 100; % Maximum number of iterations for running the Newton method.

% -------------Define the functions we need to find the roots -------------
f = get_n_wells_func(w,s,v0,m,N);
f_prime = get_f_prime(f,0.01,v0-0.01);

% ----------- Plot used to estimate the roots of the function -------------
if plotfig == 1
    E = linspace(0.01,v0-0.01,1000);
    
    figure(1);
    if plot_evidence_for_band_structure
        subplot(2,1,1)
    end
    plot(E, f(E));
    line([0 v0],[0 0],'Color',[0 0 0]);
    title({'f(E)';...
        'Roots are bound state energies for the N-square well potential';
        ['N=',num2str(N),'; v0=',num2str(v0),'; width=',num2str(w), '; separation=',num2str(s)]});
    ylabel('\bf{f(E)}');
    xlabel('\bf{E}');
    ylim([-10 , 10]);
    grid minor;
end

% ---------------- Find all roots using Newton's method -------------------
E_bound = find_roots(f, f_prime, 0.01,v0-0.01, N, tol, max_iter);

if plot_evidence_for_band_structure
    figure(1);
    subplot(2,1,2);
    for E = E_bound
        line([0 v0], [E E])
    end
    title({'Visual Evidence for band structure of the N-square well';
        ['N=',num2str(N),'; v0=',num2str(v0),'; width=',num2str(w), '; separation=',num2str(s)]})
    ylabel('Bound Energy (eV)')
end

end


