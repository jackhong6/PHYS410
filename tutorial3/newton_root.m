function [ root, num_iter ] = newton_root( f, f_prime, z0, tol, max_iter )
% Use the Newton-Raphson formula to find the root of f_z

root = z0;
num_iter = 0;
error = abs(feval(f,z0));

while error > tol && num_iter < max_iter
    root = root - feval(f,root)/feval(f_prime, root);
    error = abs(feval(f,root));
    num_iter = num_iter + 1;
end


end

