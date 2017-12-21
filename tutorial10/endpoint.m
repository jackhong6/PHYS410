function [ v1 ] = endpoint(C)
% Return v(1)
    f = @(t, y) [y(2), -exp(y(1))];
    tspan = [0, 1];
    y0 = [0, C];
    
    y = RK4(f, tspan, y0);
    v1 = y(end,1);
end