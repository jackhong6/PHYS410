function f_prime = get_f_prime(f,a,b)
    h = 1e-4;
    x = a:h:b;
    diff_est = diff(f(x))/h;
    f_prime = @(E) diff_est(floor((E-a)/h));
end
