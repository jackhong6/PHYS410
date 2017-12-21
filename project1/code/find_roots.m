function [ roots ] = find_roots( f, f_prime, a, b, N, tol, max_iter, n)
% Return the first n roots of the function f on the interval [a,b]
% using Newton's method. If n is not specified, return all roots in [a,b].
% If there are fewer than n roots, return all the roots found.
% N is the number of wells, used to determine the search increment for
% getting the guesses.
% 
% Author: Jack Hong, 30935134
% Last Modified: September 29, 2016

if nargin == 8
    guesses = get_guesses(f,a,b,N,n);
else
    guesses = get_guesses(f,a,b,N);
end

roots = zeros(1,length(guesses));
index = 1; % Index of E_bound to store the energy found
for guess = guesses
    [root, num_iter] = newton_root(f, f_prime, guess, tol, max_iter);
    if num_iter < max_iter % If the root was actually found by newton_root,
        roots(index) = root; % assign the root found to the slot in E_bound
    else % otherwise, if the root was not found,
        roots(index) = NaN; % assign NaN to the spot in E_bound
    end
    index = index + 1;
end


end

