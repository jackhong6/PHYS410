clear

f_z = @(z) z^3 - 1;
f_z_prime = @(z) 3*z^2;

xx = linspace(-1,1,500);
[xx, yy] = ndgrid(xx, xx);

for ii = 1:500
    for jj = 1:500
        [roots(ii,jj), n_iters(ii,jj)] = newton_root(f_z, f_z_prime, xx(ii) + yy(jj)*i, 1e-8, 1e6);
    end
end