% PHYS 410: Project 2
% Jack Hong, 30935134
% November 11, 2016
%
% Calculate the average energy and magnetization for the 2D Ising model
% T=Temperature, N=lattice side length, J=Ising coupling, t=number of steps
function [E,M,Cv,X,ElistF,MlistF]=ising2D(T,N,J,t,plot_flag)

%%  Initial configuration
grid = sign(.5-rand(N));  % Random initial configuration

%%  Initiation
Elist=zeros(t,1); Mlist=zeros(t,1);
lr_neighbours = grid.*circshift(grid,1,1);
ud_neighbours = grid.*circshift(grid,1,2);
Energy = -J*( sum(lr_neighbours(:))+sum(ud_neighbours(:)) );  %initial Energy
Magnet=sum(grid(:));      %initial magnetization
trials_x = randi(N,t,1);  %cheaper to generate all at once
trials_y = randi(N,t,1);

%%  Metropolis algorithm
for ii=1:t
    x = trials_x(ii); 
    y = trials_y(ii); 
    
    if x~=1; left_ind=x-1; else, left_ind=N; end
    if x~=N; right_ind=x+1; else, right_ind=1; end
    if y~=1; top_ind=y-1; else, top_ind=N; end
    if y~=N; down_ind=y+1; else, down_ind=1; end
    
    left = grid(left_ind,y); right = grid(right_ind,y);
    up = grid(x, top_ind); down = grid(x, down_ind);
    dE=2*J*grid(x,y)*(left+right+up+down);  % change in energy
    
    % Acceptance test (including the case dE<0).
    p = exp(-dE/T);
    if rand <= p
        grid(x,y) = -1*grid(x,y); 
        Energy=Energy+dE;
        Magnet=Magnet+2*grid(x,y);
    end
    
    % Update energy and magnetization.
    Mlist(ii) = Magnet;
    Elist(ii) = Energy;
    
    %Refresh display of spin configuration every N trials.
    %if mod(ii,N) && plot_flag==1, pcolor(grid); drawnow;end
end

%% Display time series of energy and magnetization
Elist(Elist==0)=[];Mlist(Mlist==0)=[];
Mlist=abs(Mlist);
Mlist=Mlist/N; Elist=Elist/N;    %normalize.
if plot_flag==1
    figure; 
    subplot(2,1,1)
    plot(Elist)
    subplot(2,1,2)
    plot(Mlist)
end

%%  Magnetization and energy density
% Eliminate all configurations before thermalization.
MlistF = Mlist; ElistF = Elist;
Mlist(1:2e4)=[]; Elist(1:2e4)=[];

% Average over post thermalization configurations.
M = mean(Mlist);
E = mean(Elist);
Cv = var(Elist)/T^2;
X = var(Mlist)/T;
end

