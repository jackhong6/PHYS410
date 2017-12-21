function [ f ] = get_n_wells_func( w,s,v0,m,N )
% Return a function f(E), for N square wells with width, w (in nm),
% separation, s (in nm), depth, v0 (in eV), and particle mass, m (in m_e).
%
% Author: Jack Hong, 30935134
% Last Modified: Sept. 30, 2016

hbar2 = calculate_hbar2;   % approx.  0.0762 eV*m_e*nm^2

alpha = @(E) sqrt(2*m*E/hbar2);
beta = @(E) sqrt(2*m*(v0-E)/hbar2);

P_allowed = @(E,d) [cos(alpha(E).*d), 1./alpha(E).*sin(alpha(E).*d); ...
                    -alpha(E).*sin(alpha(E).*d), cos(alpha(E).*d)];
                
P_forbidden = @(E,d) [cosh(beta(E).*d), 1./beta(E).*sinh(beta(E).*d); ...
                      beta(E).*sinh(beta(E).*d), cosh(beta(E).*d)];
                  
function psi_v = propagate(E)    
    psi_v = P_allowed(E,w) * [1; beta(E)];
    for i = 1:N-1
        psi_v = feval(P_allowed,E,w) * feval(P_forbidden,E,s) *  psi_v;
    end
end

function f_evals = get_f(E)
    f_evals = zeros(1, length(E));
    index = 1;
    for Ei = E
        f_evals(index) = [beta(Ei) 1] * propagate(Ei);
        index = index + 1;
    end
end

f = @(E) get_f(E);

end

