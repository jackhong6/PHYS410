function [] = tutorial12(s, p)
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
u1(2, :) = dt*psi1(x) ...
    + s/2*(phi1(circshift(x,-1))+phi1(circshift(x,1))) ...
    + (1-s)*phi1(x);

u2 = zeros(length(t), length(x));
u2(1, :) = phi2(x);
u2(2, :) = dt*psi2(x) ...
    + s/2*(phi2(circshift(x,-1))+phi2(circshift(x,1))) ...
    + (1-s)*phi2(x);

for k = 3:length(t)
    u1(k, :) = s*(circshift(u1(k-1, :), -1) + circshift(u1(k-1, :), 1)) ...
        + 2*u1(k-1, :)*(1-s) - u1(k-2, :);
    u2(k, :) = s*(circshift(u2(k-1, :), -1) + circshift(u2(k-1, :), 1)) ...
        + 2*u2(k-1, :)*(1-s) - u2(k-2, :);
end

for r = 1:length(t)
    if p==1, plot(x, u1(r, :)); end
    if p==2, plot(x, u2(r, :)); end
    ylim([-1, 1])
    xlim([-8, 8])
    pause(0.03)
end



end