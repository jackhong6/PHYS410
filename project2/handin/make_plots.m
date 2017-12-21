function fig = make_plots(N)
fig = figure;
temps = 2:0.02:2.6;
ntemps = size(temps,2);

Es=zeros(ntemps,1); Ms=zeros(ntemps,1); Cvs=zeros(ntemps,1); Xs=zeros(ntemps,1);
for k = 1:ntemps
    [Es(k),Ms(k),Cvs(k),Xs(k),Elist,Mlist] = ising2D(temps(k), N, 1, 5e4*N^2, false);
end
subplot(2,2,1)
plot(temps, Es)
title('Energy')
xlabel('Temperature')
legend(['N = ' num2str(N)], 'Location', 'best')
xlim([2, 2.6])

subplot(2,2,2)
plot(temps, Ms)
title('Magnetization')
xlabel('Temperature')
legend(['N = ' num2str(N)], 'Location', 'best')
xlim([2, 2.6])

subplot(2,2,3)
plot(temps, Cvs)
title('Heat Capacity')
xlabel('Temperature')
legend(['N = ' num2str(N)], 'Location', 'best')
xlim([2, 2.6])

subplot(2,2,4)
plot(temps, Xs)
title('Magnetic Susceptibility')
xlabel('Temperature')
legend(['N = ' num2str(N)], 'Location', 'best')
xlim([2, 2.6])
end