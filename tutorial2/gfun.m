% Defining fucntions, way 2
% Can be located either after finite_well.m, or
% as its own file in the same directory.
function [g_even,g_odd] = gfun(E,a,m,v0,hbar2)

alpha = sqrt(2*m*E/hbar2);
beta  = sqrt(2*m*(V0-E)/hbar2);

g_even = alpha.*tan(alpha*a) - beta;
g_odd = alpha.*cot(alpha*a) + beta;

end