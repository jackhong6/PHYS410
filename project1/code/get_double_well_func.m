function f = get_double_well_func( w, m, v0, s)
% Return a symfun, f(E), for a double square well with width, w (in nm),
% mass, m (in m_e), depth, v0 (in eV), and separation, s (in nm).
%
% Author: Jack Hong, 30935134
% Last Modified: Sept. 30, 2016

hbar2 = calculate_hbar2;   % approx.  0.0762 eV*m_e*nm^2
syms E d
alpha = symfun( sqrt(2*m*E/hbar2), E );
beta = symfun( sqrt(2*m*(v0-E)/hbar2), E );

psi_init = symfun([1; beta(E)], E);

P_allowed = symfun([cos(alpha(E).*d), 1./alpha(E).*sin(alpha(E).*d); ...
                    -alpha(E).*sin(alpha(E).*d), cos(alpha(E).*d)], [E,d]);
                
P_forbidden = symfun([cosh(beta(E).*d), 1./beta(E).*sinh(beta(E).*d); ...
                      beta(E).*sinh(beta(E).*d) cosh(beta(E).*d)], [E,d]);
                
psi_vector = P_allowed(E,w)*P_forbidden(E,s)*P_allowed(E,w)*psi_init(E); 
psi = symfun(psi_vector(1), E);
psi_prime = symfun(psi_vector(2), E);

f = symfun(psi_prime(E)+beta(E).*psi(E), E);

end

