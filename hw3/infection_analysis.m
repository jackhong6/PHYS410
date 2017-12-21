% PHYS 410: Homework 3
% Jack Hong, 30935134
% October 28, 2016
% NOTE: May take around 30 minutes to run

function [ lopts, mopts, hopts ] = infection_analysis()
% Analysis of 2D random walk simulation (infection.m)
%    Plot points near optimal points obtained from prelim_analysis.m

tic
nruns = 50;  % number of times to run the simulation at each point
npnts = 31;

    function [costs, errors] = get_costs(density, pvacs)
        costs = zeros(1, npnts);
        errors = zeros(1, npnts);
        for k = 1:npnts
            pvac = pvacs(k);
            runs = zeros(1, nruns);
            for m = 1:nruns
                runs(m) = infection(density, pvac, false);
            end
            costs(k) = mean(runs);
            errors(k) = std(runs);  % use standard deviation as errors
        end
    end

lpvacs = linspace(0, 0.3, npnts);
mpvacs = linspace(0.25, 0.55, npnts);
hpvacs = linspace(0.45, 0.75, npnts); 


[lden_costs, lden_errs] = get_costs(0.25, lpvacs);
[mden_costs, mden_errs] = get_costs(0.35, mpvacs);
[hden_costs, hden_errs] = get_costs(0.5, hpvacs);

[lcost, ind] = min(lden_costs);
lerror = lden_errs(ind);
lopts = lpvacs( lden_costs-lden_errs < lcost+lerror );

[mcost, ind] = min(mden_costs);
merror = mden_errs(ind);
mopts = mpvacs( mden_costs-mden_errs < mcost+merror );

[hcost, ind] = min(hden_costs);
herror = hden_errs(ind);
hopts = hpvacs( hden_costs-hden_errs < hcost+herror );


figure(1)
errorbar(lpvacs, lden_costs, lden_errs, '.')
xlabel('Percent of Population Vaccinated')
ylabel('Total Costs = Vaccines + Cures ($)')
title('Density = 0.25')
grid minor

figure(2)
errorbar(mpvacs, mden_costs, mden_errs, '.')
xlabel('Percent of Population Vaccinated')
ylabel('Total Costs = Vaccines + Cures ($)')  
title('Density = 0.35')
grid minor

figure(3)
errorbar(hpvacs, hden_costs, hden_errs, '.')    
xlabel('Percent of Population Vaccinated')
ylabel('Total Costs = Vaccines + Cures ($)')
title('Density = 0.5')
grid minor

running_time = toc;

disp(horzcat('Infection analysis completed in ', ...
    num2str(running_time),' seconds.'));
end

