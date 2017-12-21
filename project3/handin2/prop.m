function n = prop( y, t, A )
% For use with RK45.m
omega=2/3; nu=1/2; m=1; g=1;
n=[y(2); A/m*sin(omega*t)-nu/m*y(2)-g*sin(y(1))];

end

