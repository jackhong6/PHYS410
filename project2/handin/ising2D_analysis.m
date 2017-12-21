%% Clear all variables and close all figures
clear all; close all;

%% Plot thermalization energy vs time
figure
[E1,M1,Cv1,X1,Elist1,Mlist1] = ising2D(3, 20, 1, 2e4, false);
[E2,M2,Cv2,X2,Elist2,Mlist2] = ising2D(3, 40, 1, 2e4, false);
plot(Elist1); hold on;
plot(Elist2);
legend('N=20','N=40')
title('Total Energy as a Function of Time for T=3, J=1, and N=20,40')
ylabel('Total Normalized Energy')
xlabel('Time')


%% Plot E, M, Cv, X as a function of temperature for N = 10
make_plots(10);

%% Plot E, M, Cv, X as a function of temperature for N = 20
make_plots(20);

%% Plot E, M, Cv, X as a function of temperature for N = 50
make_plots(50);

%% Curie Temperature
figure
temps = 2.2:0.01:2.4;
ntemps = size(temps,2);

Es=zeros(ntemps,1); Ms=zeros(ntemps,1); Cvs=zeros(ntemps,1); Xs=zeros(ntemps,1);
for k = 1:ntemps
    [Es(k),Ms(k),Cvs(k),Xs(k),Elist,Mlist] = ising2D(temps(k), 50, 1, 5e4*50^2, false);
end

subplot(2,1,1)
plot(temps, Cvs)
ylabel('Specfic Heat')
xlabel('Temperature')

subplot(2,1,2)
plot(temps, Xs)
ylabel('Magnetic Susceptibility')
xlabel('Temperature')