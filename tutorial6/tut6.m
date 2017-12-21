tol = 1e-5;
err = @(Is, I_left, I_right) 1/15 * (I_left + I_right - Is);

[I, mesh, fevals] = quade(a,b,tol,fa,fab2,fb,I_S,mesh,fevals,f);
mesh = sort(mesh);


while err(Is, I_left, I_right) > tol
    
end