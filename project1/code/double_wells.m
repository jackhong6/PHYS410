function E_bound = double_wells()
% PHYS 410: Project 1, Problem 2; Solution to exercise 2.18
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
s = 0.2;        % The separation between the two wells in nm
v0 = 10;        % The depth of the wells in eV
m = 1;          % The mass of the particle in m_e (mass of electron)

plotfig = 1;    % If plotfig == 1, a plot of the f(E) will be generated.
                
plotsfig = 1;   % If plotsfig == 1, a plot of E2-E1 vs s will be generated.
        
tol = 1e-8;     % Tolerance for when a root is "found"
max_iter = 100; % Maximum number of iterations for running the Newton method.

% -------------Define the functions we need to find the roots -------------

% ----------------- Old method using symbolic functions -------------------
% Using symbolic functions will take 2 to 3 minutes to find the roots.
% f = get_double_well_func(w,m,v0,s);
% f_prime = diff(f);
% -------------------------------------------------------------------------
f = get_n_wells_func(w,s,v0,m,2);
f_prime = get_f_prime(f,0.01,v0-0.01);

% ----------- Plot used to estimate the roots of the function -------------
if plotfig == 1
    E = linspace(0.1,v0-0.1,1000);
    
    figure(1);
    plot(E, f(E));
    ylim([-3,3]);
    line([0 v0],[0 0],'Color',[0 0 0]);
    title({'f(E)';'Roots are bound state energies for a double square well potential';...
        ['v0=',num2str(v0),'; width=',num2str(w), '; separation=',num2str(s)]});
    ylabel('\bf{f(E)}');
    xlabel('\bf{E}');
    grid minor;
end

% ---------------- Find all roots using Newton's method -------------------
E_bound = find_roots(f, f_prime, 0.01,v0-0.01, 2, tol, max_iter);

% ---------------------- Plot E1-E0 vs separation -------------------------
if plotsfig == 1
    s = 0.05:0.05:0.4;
    delta_E = zeros(1,length(s));
    index = 1;
    for si = s
        f = get_n_wells_func(w,si,v0,m,2);
        f_prime = get_f_prime(f,0.01,v0-0.01);
        roots = find_roots(f,f_prime, 0.01,v0-0.01, 2, tol, max_iter, 2);
        delta_E(index) = roots(2)-roots(1);
        index = index + 1;
    end
    
    figure(2);
    plot(s, delta_E);
    title('(E_1 - E_0) vs Separation distance');
    ylabel('\bf{(E_1 - E_0) (eV)}');
    xlabel('\bf{Separation distance (nm)}');
end

end





