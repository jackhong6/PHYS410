% PHYS 410: Homework 2, Question 2
% Jack Hong, 30935134
% October 16, 2016

function [ result ] = romberg_integration( f, a, b)
% Return the Romberg integration of function f over the interval [a,b]
% using the the composite trapezoid formula. Start with error of order h^2,
% where h=b-a, then use subintervals of length h*2^(-m) etc. until all
% errors up to order h^6 are eliminated.

m = 0;  % start with h=b-a (m=0 -> h=b-a)

    function approx = trapezoid(m)
    %Return the integral calculated using the composite trapzoid rule with 2^m intervals.
        num_intervals = 2^m;
        h = abs(b-a) / num_intervals;
        x = a+h : h : b-h;   % x1, x2, x3, ... (exclude x0 and xN)
        approx = h*( f(a)/2 + sum( f(x) ) + f(b)/2 ); % trapezoid rule
    end

% From section 4.5 of the textbook.
T0 = @(m) trapezoid(m);                   % error = O(h^2)
T1 = @(m) ( 4*T0(m+1) - T0(m) ) / 3;      % error = O(h^4)
T2 = @(m) ( 16*T1(m+2) - T1(m+1) ) / 15;  % error = O(h^6)

result = T2(m);

end