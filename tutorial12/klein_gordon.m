function [ ] = klein_gordon( g, s, p )
% PHYS 410 tutorial 12

L = 8; c = 1;
dx = 0.05;
x = -L:dx:L;
dt = sqrt(s)*dx/c;
t = 0:dt:7;

phi1 = @(x) exp(-5*x.^2);
phi2 = @(x) sin(2*pi*x/L);

psi1 = @(x) 0;
psi2 = @(x) 0;

u1 = zeros(length(t), length(x));
u1(1, :) = phi1(x);
u1(2, :) = dt*psi1(x) + 0.5*(phi1(circshift(x,-1))*(dx/dt)^2 ...
            + phi1(x).*(2-dt^2*(2/dx^2+g*phi1(x).^2+1)) ...
            + phi1(circshift(x,1))*(dt/dx)^2);

u2 = zeros(length(t), length(x));
u2(1, :) = phi2(x);
u2(2, :) = dt*psi2(x) + 0.5*(phi2(circshift(x,-1))*(dx/dt)^2 ...
            + phi2(x).*(2-dt^2*(2/dx^2+g*phi2(x).^2+1)) ...
            + phi2(circshift(x,1))*(dt/dx)^2);


for k = 3:length(t)
    u1(k, :) = circshift(u1(k-1, :), -1) ...
        + u1(k-1, :).*(2-dt^2*(2/dx^2+g*u1(k-1, :).^2+1)) ...
        + circshift(u1(k-1, :), 1)*(dt/dx)^2-u1(k-2, :);
    u2(k, :) = circshift(u2(k-1, :), -1) ...
        + u2(k-1, :).*(2-dt^2*(2/dx^2+g*u2(k-1, :).^2+1)) ...
        + circshift(u2(k-1, :), 1)*(dt/dx)^2-u2(k-2, :);
end

for r = 1:length(t)
    if p==1, plot(x, u1(r, :)); end
    if p==2, plot(x, u2(r, :)); end
    ylim([-1, 1])
    xlim([-8, 8])
    pause(0.03)
end


end

