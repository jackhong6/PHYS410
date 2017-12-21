% PHYS 410: Homework 3
% Jack Hong, 30935134
% October 28, 2016
% NOTE: May take a couple of minutes to run

function [ low_opt, low_cost, low_err, ...
           mid_opt, mid_cost, mid_err, ...
           high_opt, high_cost, high_err ] = prelim_analysis()
%    Preliminary analysis of 2D random walk simulation (infection.m)
%    Return apporx optimal fraction of vaccinated population and average cost of
%    the optimal point for low (0.25), mid (0.35), and high (0.5) density.
%    Optimal is simply defined as where the total average cost is at its
%    minimum, without taking into account the standard deviation of that
%    average.
tic
plot = true;
npnts = 21;
nruns = 10;
pvacs = linspace( 0, 1, npnts );  % pvacs = 0, 0.1, 0.2, 0.3, ..., 1.0

    function [costs, errors] = get_costs(density)
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

[lden_costs, lden_errs] = get_costs(0.25);
[mden_costs, mden_errs] = get_costs(0.35);
[hden_costs, hden_errs] = get_costs(0.5);


[low_cost, ind] = min(lden_costs);
low_opt = pvacs(ind);
low_err = lden_errs(ind);

[mid_cost, ind] = min(mden_costs);
mid_opt = pvacs(ind);
mid_err = mden_errs(ind);

[high_cost, ind] = min(hden_costs);
high_opt = pvacs(ind);
high_err = hden_errs(ind);


if plot
    figure(1)
    errorbar(pvacs, lden_costs, lden_errs, '.')
    xlabel('Percent of Population Vaccinated')
    ylabel('Total Costs = Vaccines + Cures ($)')
    title('Density = 0.25')
    
    figure(2)
    errorbar(pvacs, mden_costs, mden_errs, '.')
    xlabel('Percent of Population Vaccinated')
    ylabel('Total Costs = Vaccines + Cures ($)')  
    title('Density = 0.35')
    
    figure(3)
    errorbar(pvacs, hden_costs, hden_errs, '.')    
    xlabel('Percent of Population Vaccinated')
    ylabel('Total Costs = Vaccines + Cures ($)')
    title('Density = 0.5')
end
running_time = toc;

disp(horzcat('Infection analysis completed in ', ...
    num2str(running_time),' seconds.'));
end

