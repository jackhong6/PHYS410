% Newton-Raphson method for finding roots of functions of a single variable.
% Use as starting point, add comments and elaborations of basic code as needed.
% Currently insufficiently documented, also good to add more safeguards in case method diverges, and some measures of performance.
function root=newton(F, Fprime,guess)
iter=1;
maxiter=1000;
x=guess;
accuracy=1;
tolerance=1e-12;
while accuracy>tolerance && iter<maxiter;
    x=x-feval(F,x)/feval(Fprime,x);
    accuracy=abs(feval(F,x));
end
root=x;
end

function out=F(x)
out=cos(x)-x;
end

function out=Fprime(x)
out=-sin(x)-1;
end