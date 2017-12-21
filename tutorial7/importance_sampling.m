function [I1,I2] = importance_sampling(k)

f1 = @(x) x.^(-1/3) + x/10;

f2 = @(x) 3*(1 + x.^2/10)/2;

%% Uniform sampling vs importance sampling

x0 = rand(k,1);

I1 = mean(f1(x0));
I2 = mean(f2(x0));


end