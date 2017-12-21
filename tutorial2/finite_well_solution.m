function [root_even, root_odd] = finite_well_solution(a,m,V0)

%% Definitions

hbar2 = 0.076199682;     % in units of eV*m_e*nm^2

alpha = @(E) sqrt(2*m*E/hbar2);  
beta = @(E) sqrt(2*m*(V0-E)/hbar2);

f_even = @(E) beta(E).*cos(alpha(E)*a) - alpha(E).*sin(alpha(E)*a);
f_odd = @(E) alpha(E).*cos(alpha(E)*a) + beta(E).*sin(alpha(E)*a);

g_even = @(E) alpha(E).*tan(alpha(E)*a) - beta(E);
g_odd = @(E) alpha(E).*cot(alpha(E)*a) + beta(E);

%% Question 2: Where are the bound states?

E0 = linspace(0,V0,1000);

plotfig=1;
if plotfig==1
    figure(1)
    
    subplot(2,1,1)
    plot(E0,f_even(E0),E0,g_even(E0)); ylim([-max(f_even(E0)) max(f_even(E0))])
    line([0 V0],[0 0],'Color',[0 0 0])
    title('Finding the even roots'); xlabel('Energy'); ylabel('Transcendental equation');
    legend('f(E)','g(E)')
    
    subplot(2,1,2)
    plot(E0,f_odd(E0),E0,g_odd(E0)); ylim([-max(f_even(E0)) max(f_even(E0))])
    line([0 V0],[0 0],'Color',[0 0 0])
    title('Finding the odd roots'); xlabel('Energy'); ylabel('Transcendental equation');
    legend('f(E)','g(E)')
end

% How many bound states are there?
n_even = sum(abs(diff(sign(f_even(E0))))/2);
n_odd = sum(floor(abs(diff(sign(f_odd(E0))))/2));

%% Question 3: Lowest eigenvalues

%---------------------------
% Step 1: Calculate the derivatives (by hand)
%---------------------------

fprime_even = @(E) (-m/hbar2)*(1./beta(E)+a).*cos(alpha(E)*a) - ...
                    (m/hbar2)*(beta(E)*a+1)./alpha(E).*sin(alpha(E)*a);
fprime_odd = @(E) (m/hbar2)*(1+beta(E)*a)./alpha(E).*cos(alpha(E)*a) - ...
                    (m/hbar2)*(1./beta(E)+a).*sin(alpha(E)*a);
       
%---------------------------
% Step 2: Find roots using Newton method
%---------------------------

iter=1; maxiter = 1000; accuracy=1; tolerance=1e-8;    % parameters

% Find initial guess

find_even = abs(diff(sign(f_even(E0))))/2;  % logical array
find_odd = abs(diff(sign(f_odd(E0))))/2;    % logical array
n1 = min(find(find_even==1));
n2 = min(find(find_odd==1));

if isempty(n2)==1 % there is no odd bound state
    E0_even_initial = E0(n1); E0_odd_initial = 0;
else
    E0_even_initial = E0(n1); E0_odd_initial = E0(n2);  % initial guesses
end

% Newton method: for the even root

change = 1;
E0_even1 = E0_even_initial;
while accuracy > tolerance && iter < maxiter && change > tolerance
    E0_even = E0_even1 - f_even(E0_even1)/fprime_even(E0_even1);
    accuracy = abs(f_even(E0_even));
%     change = abs(E0_even - E0_even1);
    iter = iter + 1;
    E0_even1 = E0_even;
end
root_even = E0_even;

% Newton method: for the odd root

E0_odd1 = E0_odd_initial;
if E0_odd1 == 0 
    root_odd = 0; % there is no odd bound state
else
    iter=1; accuracy=1; change=1;  % reinitialize iter and accuracy
        while accuracy > tolerance && iter < maxiter && change > tolerance
            E0_odd = E0_odd1 - f_odd(E0_odd1)/fprime_odd(E0_odd1);
            accuracy = abs(f_odd(E0_odd));
            change = abs(E0_odd - E0_odd1);
            iter = iter + 1;
            E0_odd1 = E0_odd;
        end
    root_odd = E0_odd;
end 

%---------------------------
% Step 3: Verification of accuracy
%---------------------------

epsilon = tolerance/2;
% fe1 and fe2 should be of different signs
fe1 = f_even(root_even + epsilon);
fe2 = f_even(root_even - epsilon);
sign_even = sign(fe1*fe2)<0;    % true if 1
% fo1 and fo2 should be of different signs
fo1 = f_odd(root_odd + epsilon);
fo2 = f_odd(root_odd - epsilon);
sign_odd = sign(fo1*fo2)<0;     % true if 1

%---------------------------
% Step 4: Plotting the eigenvalues and eigenvectors
%---------------------------

% Define the domain in a piecewise fashion
x0 = 3*a;
x1 = linspace(-x0,-a,100);
x2 = linspace(-a,a,200);
x3 = linspace(a,x0,100);
x = [x1 x2 x3];

% Define the even wavefunction
C_even = cos(alpha(root_even)*a)*exp(beta(root_even)*a);
psi_even1 = C_even*exp(beta(root_even)*x1) + root_even;
psi_even2 = cos(alpha(root_even)*x2) + root_even;
psi_even3 = C_even*exp(-beta(root_even)*x3) + root_even;
psi_even = [psi_even1 psi_even2 psi_even3];

% Define the odd wavefunction
C_odd = sin(alpha(root_odd)*a)*exp(beta(root_odd)*a);
psi_odd1 = -C_odd*exp(beta(root_odd)*x1) + root_odd;
psi_odd2 = sin(alpha(root_odd)*x2) + root_odd;
psi_odd3 = C_odd*exp(-beta(root_odd)*x3) + root_odd;
psi_odd = [psi_odd1 psi_odd2 psi_odd3];

% Plot the wavefunctions
plotfig=1;
if plotfig==1
    figure(2)
    
    hold on
    plot(x,psi_even,'-b',x,psi_odd,'-r')
    plot(x1,V0*ones(size(x1)),'-k',x2,zeros(size(x2)),'-k',x3,V0*ones(size(x3)),'-k'); 
    ylim([-0.5 V0+0.5])
    hold off

    line([-a -a],[0 V0],'Color',[0 0 0])
    line([a a],[0 V0],'Color',[0 0 0])

    title('Wavefunctions of the two lowest eigenvalues for the finite well')
    xlabel('Position'); ylabel('Energy');
    legend('Even eigenstate','Odd eigenstate')
end

%% Answers for other questions

% Question 4.2: V0* = 1.04 eV
% Question 4.3: a* = 0.0968 nm
% Question 5.1b: Edelta = V0-m*(2*a*V0)^2/(2*hbar2);


end