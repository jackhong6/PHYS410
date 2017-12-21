% PHYS 410: Project 3 ? Driven damped pendulum motion
% Jack Hong, 30935134
% November 22, 2016

clear all; close all;

%% Define constants
l = 1;  % m
m = 1;  % kg
g = 1;  % m/s^2

%% Plot unforced motion of pendulum with nu=1,5,10.

nus = [1, 3, 5, 10];  % Frequency of force

figure 
subplot(2, 1, 1)
ylabel('$\theta(t)\ (radians)$', 'Interpreter', 'latex')
xlabel('Time\ (s)', 'Interpreter', 'latex')
hold on

subplot(2, 1, 2)
xlabel('$\theta(t)\ (radians)$', 'Interpreter', 'latex')
ylabel('$v(t)\ (s^{-1})$', 'Interpreter', 'latex')
hold on

for nu = nus
    f = @(t, y) [y(2), -nu/m*y(2)-g*sin(y(1))];  % Return [v, alpha]
    tspan = [0, 60];  % Time interval interval over which we solve the ODE
    y0 = [0.2, 0];  % Initial conditions of the pendulum. [theta, omega]
    [ylist, t] = RK4(f, tspan, y0);
    ylist(:,1) = mod(ylist(:,1)+pi, 2*pi);

    subplot(2,1,1)
    plot(t, ylist(:, 1), 'DisplayName', ['$\nu=' num2str(nu) '$'])
    
    subplot(2,1,2)
    plot(ylist(:, 1), ylist(:, 2), 'DisplayName', ['$\nu=' num2str(nu) '$'])
end

subplot(2, 1, 1)
legend1 = legend(gca, 'show');

subplot(2, 1, 2);
legend2 = legend(gca, 'show');

set(legend1, 'Interpreter', 'latex'); set(legend2, 'Interpreter', 'latex')

%% Plot forced motion of pendulum with nu=1/2, omega=2/3, A=0.5, 1.2
nu = 1/2;
omega = 2/3;

A = 0.5;
f = @(t, y) [y(2), A/m*sin(omega*t)-nu/m*y(2)-g*sin(y(1))];  % Return [v, alpha]
tspan = [0, 120];  % Time interval interval over which we solve the ODE
y0 = [0.2, 0];  % Initial conditions of the pendulum. [theta, v]
[ylist, t] = RK4(f, tspan, y0);
ylist(:,1) = mod(ylist(:,1)+pi, 2*pi);

figure
subplot(2, 1, 1)
plot(t, ylist(:, 1), 'Displayname', ['A=' num2str(A)])
ylabel('$\theta(t)\ (radians)$', 'Interpreter', 'latex')
xlabel('Time\ (s)', 'Interpreter', 'latex')
title(['A=' num2str(A)])

subplot(2, 1, 2)
plot(ylist(:, 1), ylist(:, 2), 'Displayname', ['A=' num2str(A)])
xlabel('$\theta(t)\ (radians)$', 'Interpreter', 'latex')
ylabel('$v(t)\ (s^{-1})$', 'Interpreter', 'latex')

% ===============================

A = 1.2;
y0 = [0.2; 0];  % Initial conditions of the pendulum. [theta, v]
[ylist, t] = RK45(y0, 300, 1.2);
ylist(1,:) = mod(ylist(1,:)+pi, 2*pi);

figure
subplot(2, 1, 1)
plot(t, ylist(1, :), 'Displayname', ['A=' num2str(A)])


subplot(2, 1, 2)
plot(ylist(1, :), ylist(2, :), 'Displayname', ['A=' num2str(A)])
xlabel('$\theta(t)\ (radians)$', 'Interpreter', 'latex')
ylabel('$v(t)\ (s^{-1})$', 'Interpreter', 'latex')

%% Plot theta(t) for A=1.35, 1.44, 1.465, and their Poincare sections
As = [1.35, 1.44, 1.465];

for A = As
    y0 = [0.2; 0];  % Initial conditions of the pendulum. [theta, v]
    [ylist, t] = RK45(y0, 120, A);
    ylist(1,:) = mod(ylist(1,:)+pi, 2*pi);
    
    figure
    plot(t, ylist(1, :))
    ylabel('$\theta(t)\ (radians)$', 'Interpreter', 'latex')
    xlabel('Time\ (s)', 'Interpreter', 'latex')
    title(['A=' num2str(A)])
end

%% Plot Poincare section
As = [1.35, 1.44, 1.465];
omega = 2/3;

figure

for A = As
    y0 = [0.2; 0];  % Initial conditions of the pendulum. [theta, v]
    [ylist, t] = RK45(y0, 3000, A);
    ylist(1,:) = mod(ylist(1,:)+pi, 2*pi);
    
    tol = 1e-4;
    n = t.*(omega/(2*pi));
    poincare_pnts = abs(round(n)-n) < tol;
    %t(poincare_pnts) .* (omega/2/pi)
    scatter(ylist(1, poincare_pnts), ylist(2, poincare_pnts), 'DisplayName', ['A=' num2str(A)])
    hold on
end
ylim([-2, -1.2])
legend(gca, 'show')
xlabel('$\theta(t)\ (radians)$', 'Interpreter', 'latex')
ylabel('$v(t)\ (s^{-1})$', 'Interpreter', 'latex')
title('Poincare Section')

