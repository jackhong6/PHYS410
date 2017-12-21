function shooting

% Rewriting the BVP
f = @(x,u) [u(2), 3*u(2)/(1-x) + u(1)/(x*(1-x))];

% Finding C0
C0 = fzero(@endpoint,-1/6);     % lim_{x->1} u'(x) = -1/6 

epsilon = 1e-9;                 % need cutoff at x=1
[u,x] = RK4(f,[1-epsilon 0],[1/2 C0]);

% Plot numerical solution against analytic solution
xx = linspace(0,1,1000);
u0 = (1-xx+xx.*log(xx))./(xx-1).^2;

plot(x,u(:,1),'.',xx,u0); 
ylim([0.49 1.01]); xlim([0 1]);

end


function u1 = endpoint(C)

f = @(x,u) [u(2), 3*u(2)/(1-x) + u(1)/(x*(1-x))];

epsilon = 1e-9;                  % need cutoff at x=1
u = RK4(f,[1-epsilon 0],[1/2 C]);

u1 = u(end,1) - 1;               % so that u(0) = 1

end

function [y,t] = RK4(f,tspan,y0)

% f: function evaluating the time derivatives; y'(t) = f(t,y)
% Note that f must accept TWO arguments: t and y(t). For systems of N first  
% order ODEs, f and y should both output row vectors with N components.
% y0: row vector of initial conditions
% tspan = [t0 tf] is the time interval over which we solve the ODE

%% Parameters

h = -1e-2;
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