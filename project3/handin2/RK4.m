% Explicit Fourth-Order Runge-Kutta with constant step size

function [y,t] = RK4(f,tspan,y0)

% f: function evaluating the time derivatives; y'(t) = f(t,y)
% Note that f must accept TWO arguments: t and y(t). For systems of N first  
% order ODEs, f and y should both output row vectors with N components.
% y0: row vector of initial conditions
% tspan = [t0 tf] is the time interval over which we solve the ODE

%% Parameters

h = 1e-2;
t0 = tspan(1); tf = tspan(2);

%% Initilization

iter = round((tf-t0)/h);        % number of time steps
y = zeros(iter+1,length(y0));   % preallocating
y(1,:) = y0;                    % y(t0) = y0
    
for i=1:iter
    k1 = feval(f, t0 + (i-1)*h,       y(i,:)            );
    k2 = feval(f, t0 + (i-1)*h + h/2, y(i,:) + (h/2)*k1 );
    k3 = feval(f, t0 + (i-1)*h + h/2, y(i,:) + (h/2)*k2 );
    k4 = feval(f, t0 + (i-1)*h + h,   y(i,:) + h*k3     );
    y(i+1,:) = y(i,:) + (h/6)*(k1+2*k2+2*k3+k4);
end

t = t0 + h*(0:iter)';

end