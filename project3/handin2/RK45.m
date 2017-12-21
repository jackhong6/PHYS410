% Runge-Kutta-Fehlberg with adaptive step size.
function [dynamicallist,times]=RK45(initial,tfinal, A)
% prop: function evaluating the time derivatives.
% Initial: vector of initial conditions.
% tfinal is the final time.
% Produces the values of the dynamical variables and the 
% corresponding times.
%% Parameters
hmin=1e-6;
hmax=1e-2;
abstol=1e-4;
reltol=1e-4;
%% Initiation
iter=round(tfinal/hmax);
times=zeros(1,iter+1); %vector of times when solution is obtained.
dynamicallist=zeros(length(initial),iter+1);dynamicallist(:,1)=initial;
h=sqrt(hmax*hmin);i=1;
%% Constants
a=zeros(6,6);
a(2,1)=1/4;
a(3,1)=3/32;a(3,2)=9/32;
a(4,1)=1932/2197;a(4,2)=-7200/2197; a(4,3)=7296/2197;
a(5,1)=439/216;a(5,2)=-8;a(5,3)=3680/513; a(5,4)=-845/4104;
a(6,1)=-8/27;a(6,2)=2;a(6,3)=-3544/2565;a(6,4)=1859/4104;a(6,5)=-11/40;
b(1)=16/135;b(2)=0;b(3)=6656/12825;b(4)=28561/56430;b(5)=-9/50;b(6)=2/55;
c(1)=0;c(2)=1/4;c(3)=3/8;c(4)=12/13;c(5)=1;c(6)=1/2;
r(1)=1/360;r(2)=0;r(3)=-128/4275;r(4)=-2197/75240;r(5)=1/50;r(6)=2/55;
%% Main loop
k=zeros(6,length(initial));
while (times(i) < tfinal )
    for j=1:6;
        ak=h*(a*k)';
        k(j,:)= prop(initial+ak(:,j), times(i)+c(j)*h, A);
    end
    incr=(b*k)';
    error=norm(r*k)/(abstol+reltol*(norm(incr)));
    if (error<1)
        initial = initial + h*incr;
        % Updating lists
        dynamicallist(:,i+1)=initial;
        times(i+1) = times(i)+h; i = i + 1;
        rej=0;
    else rej=1;% Rejected step
    end;
    % adjusting step size
    if rej==1; facmax=1; else facmax=4;end
    facmin=.2;
    h=h*min(max((.38/error )^(1/5),facmin),facmax);
    if(h>hmax); h=hmax;end;
    if(times(i)+h>tfinal),h=tfinal-times(i);
    elseif (h<hmin),disp ('Solution requires step size smaller than minimum'); return;end;
end
%Pruning the final lists

dynamicallist(:, i+1:end)=[];
times(i+1:end)=[];
end