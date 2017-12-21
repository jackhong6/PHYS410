% Function to calculate average energy and magnetization for the 1D Ising model
% T = Temperature, N = linear lattice size, J = Ising coupling.
function [E,M]=ising1D(T,N,J,plot_flag)
%%  Initial configuration
grid = sign(.5-rand(N,1)); % Random initial configuration
%%  Initiation
t =1e4*N;% Number of steps
Elist=zeros(t,1); Mlist=zeros(t,1);
Energy = -J*sum(grid.*circshift(grid,1));%initial Energy
Magnet=sum(grid); %initial magnetization
trials = randi(N,t,1); %cheaper to generate all at once
%%  Metropolis algorithm
for i=1:t,
    s=trials(i); 
    if s~=1; left=grid(s-1);else left=grid(N);end
    if s~=N; right=grid(s+1);else right=grid(1);end
    dE=2*J* grid(s)*(left+right);  % change in energy
    p = exp(-dE/T);
    % Acceptance test (including the case dE<0).
    if rand <= p,
        grid(s) = -1*grid(s); 
        Energy=Energy+dE;
        Magnet=Magnet+2*grid(s);
    end
    % Update energy and magnetization.
    Mlist(i) =Magnet;
    Elist(i) =Energy;
    %Refresh display of spin configuration every N trials.
             %if mod(i,N)==0 && plot_flag==1;
             %    bar(grid); drawnow;
             %end
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
%Mlist(1:50*N)=[]; Elist(1:50*N)=[];
% Average over post thermalization configurations.
M=sum(Mlist)/numel(Mlist);
E=sum(Elist)/numel(Elist);
end

