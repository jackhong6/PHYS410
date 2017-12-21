function [ guesses ] = get_guesses( f, a, b, N, n )
% Use the diff and sign functions to get the first n guesses for 
% the roots of f in the interval [a,b]. If there are fewer than n roots, 
% then return all the roots found. N is the number of wells, which 
% determines the search increment of the algorithm.
%
% f is assumed to be a continuous function that takes exactly one scalar 
% input and outputs exactly one scalar output.
%
% Author: Jack Hong, 30935134
% Last Modified: Sept. 30, 2016

incr = 0.01/N; % The search increment for a system of N wells
x = a:incr:b;

if nargin == 5 % if the caller specified n
    find_root = abs(diff(sign(f(x))))/2;
    guesses = x( find(find_root==1,n) );
else
    guesses = x( abs(diff(sign(f(x)))) > 0 );
end



