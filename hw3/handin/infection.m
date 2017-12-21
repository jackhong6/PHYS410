% PHYS 410: Homework 3
% Jack Hong, 30935134
% October 28, 2016

function [cost, ncured, V, iter] = infection(density, pvac, plot)
% 2D random walk simulation of virus infection.
% pvac -- percent of people initially vaccinated.
% density -- density of the population.
% plot -- flag for whether or not to plot 2D visualization

%tic

% Definitions =============================================================
N = 128;          % grid size: N*N
T = 15;           % hours remaining before cure
Nitermax = 5000;  % max iterations

prob = 0.66;      % probability of infection

cost_per_vac = 30;    % in dollars
cost_per_cure = 100;  % in dollars

M = floor(density*N^2);     % number of walkers
V = floor(pvac*M);          % number vaccinated (initially immune)
% =========================================================================

infect = zeros(M,1);    
% infect takes the values: 0 (healthy); 1 (sick); 2 (immunized)
hours = T*ones(M,1);     % hours gets initialized to T upon infection and 
% decreases by 1 every iteration until hours==0, which implies recovery

% Initialize the walkers' starting positions randomly
x = randi(N,M,1);
y = randi(N,M,1);

% Vaccination campaign for V walkers
infect(randperm(M,V)) = 2;

% Infect a random walker
if pvac == 1    % 100% vaccination rate
    nsick = 0;  % so no one is sick
else
    while true  % infect a random person who is not immunized
        z = randi(M);
        if infect(z) ~= 2
            infect(z) = 1;
            break;
        end
    end
    nsick = 1;
end

nimmunized = V;     % only vaccinated people immunized initially
iter = 1;           % about to start first iteration

while nsick > 0 && iter < Nitermax
    
    % Count infection countdowns and cures before taking a random step
    hours(infect==1) = hours(infect==1) - 1;
    
    % Is the sick person cured?
    infect(hours==0) = 2;     
    nimmunized = sum(infect==2);
    
    % Random walk
    rand_walk = rand(M,1);
    walk_x = rand_walk < 0.5;
    walk_y = rand_walk >= 0.5;
    
        % Walk in the x direction
        deltax = 2*floor(2*rand(size(x))) - 1;  % random directions
        deltax(infect==3) = 0;  % the dead don't move
        deltax(walk_y) = 0;     % don't walk in y!
        xold = x;
        xnew = x + deltax; 
        % Periodic boundary conditions
        xnew(xnew==0) = N;
        xnew(xnew==N+1) = 1;
        x = xnew;
        
        % Walk in the y direction
        deltay = 2*floor(2*rand(size(x))) - 1;  % random directions
        deltay(infect==3) = 0;  % the dead don't move
        deltay(walk_x) = 0;     % don't walk in x!
        yold = y;
        ynew = y + deltay;
        % Periodic boundary conditions
        ynew(ynew==0) = N;
        ynew(ynew==N+1) = 1;
        y = ynew;
    
    % Contamination step 
    H_coord = zeros(M,2);     % because we need to keep track of positions
    pos_H = infect==0;
    H_coord(pos_H,:) = [x(pos_H), y(pos_H)];  % positions of healthy people
    
    pos_I = infect==1;
    I_coord = [x(pos_I), y(pos_I); xold(pos_I), yold(pos_I)];  
    % positions of sick people, before and after random step
    
    indH = find( ismember(H_coord,I_coord,'rows') == 1);
    % indH are the "j-positions" of the healty walkers who are
    % getting infected (infect(j) = 0 -> 1 for j in indH)
    
    infection_status = rand(size(indH)) < prob;  % true = turns into zombie
    infect(indH) = infection_status;
    nsick = sum(infect==1);
    
    % Increment
    iter = iter + 1;
    if iter == Nitermax
        disp('Maximum number of iterations reached; increase Nitermax')
    end
    
    % Visualization in 2D
    if plot == 1
        colorgrid = 4*ones(N);
        map = [ 0, 1, 0;    % healthy = green
                1, 0, 0;    % sick = red
                0, 0, 1;    % immunized = blue
                1, 1, 1];   % empty = white
        for i=1:M
            colorgrid(x(i),y(i)) = infect(i);
        end
        g = pcolor(colorgrid); colormap(map); axis square; 
        set(g,'LineStyle','none'); 
        set(gca,'XDir','normal')
        drawnow
    end

end
                            
ncured = nimmunized - V;
cost = V * cost_per_vac + ncured * cost_per_cure;

%running_time = toc;

%disp(horzcat('Simulation for zombie apocalypse completed in ', ...
%    num2str(running_time),' seconds.'));
end