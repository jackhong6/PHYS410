function [I,mesh,fevals] = adaptivequad(a,b,func,tol)
% Evaluates adaptative Simpson quadrature.

%% Initialization (crude)
fa = func(a); fb = func(b); fab2 = func((a+b)/2);
fevals = 3;     % counts function evaluations

Iab = srule(a,b,fa,fab2,fb);    % approximation to integral
mesh = [a, (a+b)/2, b]; % first mesh

%% Adaptative evaluation
maxlevel = 10;  % security feature if tolerance never met

[I,mesh,fevals] = quade(a,b,tol,fa,fab2,fb,Iab,maxlevel,fevals,mesh,func);

mesh = sort(mesh);  % sort mesh in ascending order 

x = linspace(a,b,2000);
plot(x,func(x),'-',mesh,zeros(size(mesh)),'xr')

end

function [I,mesh,fevals] = quade(a,b,tol,fa,fab2,fb,Iab,level,fevals,mesh,func)
% Adaptive evaluation of the integral. Recursive.

%% Refine the mesh for the first time
fa3b1 = func((3*a+b)/4); fa1b3 = func((a+3*b)/4);
mesh(fevals+1) = (3*a+b)/4; mesh(fevals+2) = (a+3*b)/4; fevals = fevals+2;

ab2 = (a+b)/2;
Ileft = srule(a,ab2,fa,fa3b1,fab2);
Iright = srule(ab2,b,fab2,fa1b3,fb);

%% Check answer against tolerance - otherwise, recursion

if abs(Iab - (Ileft + Iright) ) < 15*tol 
    I = Ileft + Iright;
else
    % must recursively apply the same logic
    % Step 1: check if max recursion reached
    
    if level==0
        display('Max recursion level reached.')
        I = Ileft + Iright;
    else
        
        [Ileft,mesh,fevals] = quade(a,ab2,tol/2,fa,fa3b1,fab2,...
            Ileft,level-1,fevals,mesh,func);
        
        [Iright,mesh,fevals] = quade(ab2,b,tol/2,fab2,fa1b3,fb,...
            Iright,level-1,fevals,mesh,func);
        
        I = Ileft + Iright;
    end
end

end

function val = srule(a,b,fa,fab2,fb)
% evaluates basic Simpson rule

val = (b-a)/6 * (fa + 4*fab2 + fb);

end