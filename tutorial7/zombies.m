function [nalive, ndeaths, iter] = zombies()

tic 
% -----------------
% Definitions
% -----------------

N = 100;        % grid size: N*N
T = 24;         % hours remaining before death
Nitermax = 5000;    % max iterations

density = 0.4;  % density of population

M = floor(density*N^2);     % number of walkers

plot_flag = 1;

% -----------------
% Initialization
% -----------------
x = randi(N,M,1);
y = randi(N,M,1);

infect = zeros(1,M);
hours = T*ones(1,M);

ndeaths = 0;
nzombies = 1;
iter = 1;
% ...............
%   Question 1
% ...............
infect(1) = 1;

% -----------------
% Random walk
% -----------------
while nzombies > 0 && iter < Nitermax
    
    % Count infections and deaths before taking a random step
    hours(infect==1) = hours(infect==1)-1;
    infect(hours<=0) = 2;
    ndeaths = sum(infect==2);
    
    % Random walk
    rand_walk = rand(M,1);
    walk_x = rand_walk < 0.5;
    walk_y = rand_walk >= 0.5;
    
    % Walk in the x direction
        deltax = 2*floor(2*rand(size(x))) - 1;  % random directions
        deltax(infect==2) = 0;  % the dead don't move
        deltax(walk_y) = 0;     % don't walk in y!
        xold = x;
        xnew = x + deltax; 
        xnew(xnew==0) = N;
        xnew(xnew==N+1) = 1;
        x = xnew;
    % Walk in the y direction
        deltay = 2*floor(2*rand(size(y))) - 1;  % random directions
        deltay(infect==2) = 0;  % the dead don't move
        deltay(walk_x) = 0;     % don't walk in y!
        yold = y;
        ynew = y + deltay; 
        ynew(ynew==0) = N;
        ynew(ynew==N+1) = 1;
        y = ynew;
    
    % Contamination 
    H_coord = zeros(M,2);     % because we need to keep track of position
    pos_H = infect==0;
    H_coord(pos_H,:) = [x(pos_H), y(pos_H)];  % positions of healthy people
    
    pos_I = infect==1;
    I_coord = [xold(pos_I), yold(pos_I); xnew(pos_I), ynew(pos_I)];
    
    indH = find(ismember(H_coord,I_coord,'rows') == 1);
    % indH are the location, in H_coord, of the healty walkers who are
    % getting infected
    
    infection_status = 1;
    infect(indH) = infection_status;
    nzombies = sum(infect==1);
    
    % Increment
    iter = iter + 1;
    if iter == Nitermax
        display('Maximum number of iterations reached; increase Nitermax')
    end
    
    % Visualization in 2D
    if plot_flag == 1
        colorgrid = 3*ones(N+1);
        
        % Ensure correct colors are set when all values of colorgrid are
        % the same.
        colorgrid(1,N+1) = 0;
        colorgrid(2,N+1) = 1;
        colorgrid(3,N+1) = 2;
        colorgrid(4,N+1) = 3;

        map = [ 0, 1, 0;    % healthy = green
                1, 0, 0;    % zombies = red
                0, 0, 0;    % dead = black
                1, 1, 1];   % empty = white
            
        for i=1:M
            colorgrid(x(i),y(i)) = infect(i);
        end
        g = pcolor(colorgrid); 
        colormap(map); 
        axis square; 
        set(g,'LineStyle','none');
        drawnow
    end

end
                            
nalive = M - nzombies - ndeaths; % alive and not immunized

running_time = toc;

display(horzcat('Simulation for zombie apocalypse completed in ', ...
    num2str(running_time),' seconds.'));

end