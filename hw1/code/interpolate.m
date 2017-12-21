function [ ] = interpolate(  )
% PHYS 410: Homework 1
%
% Jack Hong, 30935134
% Last modified: Ocotober 3, 2016

plot_fig = false;
plot_spline = true;
plot_barycentric = true;
plot_chebyshev = true;

plot_overlay = true;  % Overlay plot of f(x) over the interpolate plot
plot_residuals = true; % Plot residuals of f(x) and the interpolate plot

n = 1000; % Number of points used to plot
N = [10, 20, 30];

domain = [1, 2]; % if changed, must also change the chebyshev points function to match the new domain

f = @(x) exp(3.*x).*sin(20.*x.^2) ./ (1+2.*x.^2);  % 1<x<2
x = linspace(domain(1),domain(2),n);

chebyshev = @(ii,N) 1.5 + 0.5.*cos((2.*ii-1).*pi./(2.*N)); % must change to match domain

if plot_fig
    figure;
    plot(x,f(x))
    title('Function to interpolate');
    ylabel('$$f(x) = \frac{e^{3x}sin(20x^2)}{1+2x^2}, 1<x<2$$',...
        'interpreter','latex','fontsize',14);
    xlabel('x','fontsize',14);
end

if plot_spline
    figure;
    subplot(length(N)+1, 1+plot_residuals, 1);
    plot(x, f(x), 'Color', 'r');
    title('Function to interpolate: $$f(x) = \frac{e^{3x}sin(20x^2)}{1+2x^2}, 1<x<2$$',...
        'interpreter','latex','fontsize',14')
    plot_num = 2 + plot_residuals;
    for Ni = N
        inter_pnts = linspace(domain(1),domain(2),Ni);
        subplot(length(N)+1,1+plot_residuals,plot_num);
        y = spline(inter_pnts, f(inter_pnts),x);
        plot(x,y);
        title(['Cubic Spline Interpolation with N = ', num2str(Ni)]);
        ylim([-50,50])
        if plot_overlay
            hold on
            plot(x, f(x));
            ylim([-50,50])
        end
        if plot_residuals
            subplot(length(N)+1,1+plot_residuals,plot_num+1);
            plot(x, abs(f(x)-y))
            title('Residuals')
            ylim([0,100]);
        end
        plot_num = plot_num + 1 + plot_residuals;
    end
end

if plot_barycentric
    figure;
    subplot(length(N)+1, 1+plot_residuals, 1);
    plot(x, f(x), 'Color', 'r');
    title('Function to interpolate: $$f(x) = \frac{e^{3x}sin(20x^2)}{1+2x^2}, 1<x<2$$',...
        'interpreter','latex','fontsize',14')
    plot_num = 2 + plot_residuals;
    for Ni = N
        inter_pnts = linspace(domain(1),domain(2),Ni);
        subplot(length(N)+1,1+plot_residuals,plot_num);
        y = barycentric(inter_pnts, f(inter_pnts), x);
        plot(x,y);
        ylim([-50,50])
        title(['Barycentric Interpolation with N = ', num2str(Ni)]);
        if plot_overlay
            hold on
            plot(x, f(x));
            ylim([-50,50])
        end
        if plot_residuals
            subplot(length(N)+1,1+plot_residuals,plot_num+1);
            plot(x, abs(f(x)-y))
            title('Residuals')
            ylim([0,100]);
        end
        plot_num = plot_num + 1 + plot_residuals;
    end
end

if plot_chebyshev
    figure;
    subplot(length(N)+1, 1+plot_residuals, 1);
    plot(x, f(x), 'Color', 'r');
    title('Function to interpolate: $$f(x) = \frac{e^{3x}sin(20x^2)}{1+2x^2}, 1<x<2$$',...
        'interpreter','latex','fontsize',14')
    plot_num = 2 + plot_residuals;
    for Ni = N
        inter_pnts = chebyshev(1:Ni,Ni);
        subplot(length(N)+1,1+plot_residuals,plot_num);
        y = barycentric(inter_pnts, f(inter_pnts), x);
        plot(x,y);
        title(['Barycentric Interpolation with Chebyshev Points with N = ', num2str(Ni)]);
        ylim([-50,50])
        if plot_overlay
            hold on
            plot(x, f(x));
            ylim([-50,50])
        end
        if plot_residuals
            subplot(length(N)+1,1+plot_residuals,plot_num+1);
            plot(x, abs(f(x)-y))
            title('Residuals')
            ylim([0,100]);
        end
        plot_num = plot_num + 1 + plot_residuals;
    end
end


end

