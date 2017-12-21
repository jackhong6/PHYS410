% PHYS 410: Homework 2, Question 1
% Jack Hong, 30935134
% October 17, 2016

function [ Lp_xstar ] = num_diff( f, x_star )
% Return the derivative of f at x_star using a 7-points centered formula.
%    - Required weights are derived using the idea of differentiating the Lagrange
%      Interpolation of f using the 7 points around x_star.
%    - value of h optimized to minimize the combination of discretization
%      and round off errors.

h = 0.1; % optimal h
n = [-3*h, -2*h, -h, 0, h, 2*h, 3*h];
x_array = x_star + n;   % the 7 points used for Lagrange interpolation
f_array = f(x_array);   % the 7 values of f at xi

    function lagrange_poly = get_lagrange_poly(jj)
        % Return the jj_th Lagrange basis polynomial as a function handle.
        m = [1:jj-1, jj+1: 7];  % 1 <= m <= 7; m!=jj
        lagrange_poly = @(y) prod( (y-x_array(m)) ./ (x_array(jj)-x_array(m)) );
    end

    el1 = get_lagrange_poly(1);
    el2 = get_lagrange_poly(2);
    el3 = get_lagrange_poly(3);
    el4 = get_lagrange_poly(4);
    el5 = get_lagrange_poly(5);
    el6 = get_lagrange_poly(6);
    el7 = get_lagrange_poly(7);
    
    function el_array = get_el_array(x)
        % Return an array of l_j(x) - essentially just [0,0,0,1,0,0,0] when
        % called with x_star
        el_array = [el1(x), el2(x), el3(x), el4(x), el5(x), el6(x), el7(x)];
    end

    function elp_array = get_elp_array(x)
        % Return an array of 7 values l_j'(x), the derivative of the 
        % seven Lagrange polynomials at x_star
        m = [ [2:7]; [1,3:7]; [1:2,4:7]; [1:3,5:7]; [1:4,6:7]; [1:5,7]; [1:6] ];
        w = 1./(x-x_array(m));
        elp_array = get_el_array(x) .* sum( 1./(x-x_array(m)),2 );
    end
            
Lp_xstar = sum( f_array .* get_elp_array(x_star), 2 );



end

