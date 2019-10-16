%% SF2520 HT19-1 Applied Numerical Methods
%% Computer exercise 3 part 4 Built-in Matlab solvers
%% Partial differential equation of parabolic type
%% Armin Salmasi
%% The task is to compare two built-in Matlab solvers: the explicit method
%% ode23 and the implicit method ode23s, suitable for stiff problems.
clear all; close all; clc;
%% initialization
N = [10 20 40];  % number of grid points
tfinal = 2; %final time
timespan = [0 tfinal]; % time span
yfinal = 1;  % dianl xi
%% Solvers
%% Loop over number of timesteps
for i=1:size(N,2)
    %% stepsize discretization
    hy = (yfinal-0)/N(i);
    y = 0:hy:yfinal; 
    %% sparse initial value matrix
    u0 = sparse(zeros(N(i),1));
    %% sparse A matrix
    diag1 = ones(1,N(i));
    diag2 = ones(1,N(i))*(-2);
    diag3 = ones(1,N(i));
    diag4 = 2;
    A = sparse((diag(diag2(1:end),0)+diag(diag1(1:end-1),-1)+diag(diag3(1:end-1),1)));
    %% Boundary condition
    A(N(i),N(i)-1) = diag4;
    %% ode23
    tic
    options = odeset('RelTol',1e-3, 'Stats', 'on');
    [t,u] = ode23(@(t,u) f_ut(t,u,N(i),A), timespan, u0, options);
    clockode23(i) = toc;
    ntstpsode23(i) = length(t);
    tmaxode23(i) = max(diff(t));
    %% index of t==1
    lenu = length(u(:,1));
    idx = find(t>0.99 & t<1.01);
    [t0, n] = min(abs(t(idx)-ones(length(t(idx)),1)));
    idx = idx(n);
    %% adding initial values to the solution
    utmp = [ones(1, idx) zeros(1, lenu-idx)]';
    u=[utmp u];
    figure
         surf(y,t,u)
         xlabel('\xi');
         ylabel('\tau');
         zlabel('u(\tau,\xi)')
         title(['ode23, N=' num2str(N(i))])
    %% ode23s
    tic
    options = odeset('RelTol',1e-3, 'Stats', 'on');
    [t,u] = ode23s(@(t,u) f_ut(t,u,N(i),A), timespan, u0, options);
    clockode23s(i) = toc;
    ntstepsode23s(i) = length(t);
    tmaxode23s(i) = max(diff(t));
    %% index of t==1 
    lenu = length(u(:,1));
    idx = find(t>0.99 & t<1.01);
    [t0, n] = min(abs(t(idx)-ones(length(t(idx)),1)));
    idx = idx(n);
    %% adding initial values to the solution
    utmp = [ones(1, idx) zeros(1, lenu-idx)]';
    u=[utmp u];
    figure
        surf(y,t,u)
        title(['ode23s, N=' num2str(N(i))])
        xlabel('\xi');
        ylabel('\tau');
        zlabel('u(\tau,\xi)')
end
%% results
tstpsode23s = ntstepsode23s'
cputime23s = clockode23s'
dtmax23s = tmaxode23s'
tstpsode23 = ntstpsode23'
cputimeode23 = clockode23'
dtmaxode23 = tmaxode23'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function decleration
function [m]=alp(t)
%% time dependent left handside boundary condition
if (t>1)
    m=0;
else
    m=1;
end
end
function ut=f_ut(t,u,Nx,A)
%% du/dt=Au+b(t)
hx = 1/Nx;
b = zeros(Nx,1);
b(1,1) = alp(t);
ut = (1/hx^2)*(A*u+b);
end





