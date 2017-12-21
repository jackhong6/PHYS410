function E_star = find_root( f,f_prime,E0 )
tolerance = 1e-8; maxiter = 100;
iter = 0; accuracy = 1; change = 1;    % initialize parameters

while accuracy > tolerance && change > tolerance && iter < maxiter
    Enew = E0 - feval(f,E0)/feval(f_prime,E0);
    accuracy = abs( feval(f,Enew) );
    change = abs(Enew - E0);
    E0 = Enew; iter = iter + 1;
end

E_star = E0;

end

