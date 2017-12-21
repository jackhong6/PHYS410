function [ root, num_iter ] = newton_root( f, f_prime, guess, tol, max_iter )
% Use the Newton-Raphson method to find the root of f
%
% Author: Jack Hong, 30935134
% Last modified: Sept. 29, 2016

root = guess;
num_iter = 0;
error = abs(feval(f,guess));

while error > tol && num_iter <= max_iter
    root = root - feval(f,root)/feval(f_prime, root);
    error = abs(feval(f,root));
    num_iter = num_iter + 1;
end

end

