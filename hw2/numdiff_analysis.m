function numdiff_analysis

%% Absolute error 

E = @(h) abs(1 - (exp(h)-exp(-h))./(2*h) );

%% Roundoff error

% E_RO = abs(D_h - \bar{D}_h)
%      = abs( (f_h - f_-h)/(2h) - (f_h + e_h - f_-h - e_-h)/(2h) )
%      = abs( e_h - e_-h )/(2h)
% Using triangle inequality:
% E_R0 <= abs(e_h)/(2h) + abs(e_-h)/(2h)   % abs(e_h) ~= abs(e_-h)
%      <= epsilon/h                        % abs(e(x)) <= epsilon for all x

%% Interplay between truncation and roundoff

M = 1; % truncation; M arbitrary, as long as it yields a proper error bound
epsilon = eps;  % roundoff ~ machine precision 
ErrorBound = @(h) h.^2*M/6 + epsilon./h;

%% Plotting results

h = 10.^linspace(-1,-9);
loglog(h,E(h),'.',h,ErrorBound(h),'-')
title('Absolute error and error bound for numerical differentiation')
xlabel('Discretization size h'); ylabel('Absolute error E(h)')

end