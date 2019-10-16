%% SF2520 HT19-1 Applied Numerical Methods
%% Computer exercise 3 part 5 Optimized MATLAB solvers
%% Partial differential equation of parabolic type
%% Armin Salmasi
%% The task is to adjusting the the options for  ’Jacobian’ and ’Jpattern’
%% to take advantage of the fact that matrices are tridiagonal. 
%% This should increase the efficiency of ode23s. N=10,20,40
clear all; close all; clc;
%% initialization
N=[10, 20, 40];  % number of grid points
tfinal = 2; %final time
timespan = [0 tfinal]; % time span
yfinal = 1;  %final xi
%% Loop over number of timesteps
for i=1:size(N,2)
    %% stepsize discretization
    hy = (yfinal-0)/N(i);
    y = 0:hy:1;
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
    %% JPattern test
    %% set JPattern in options
    jpat = @(p) spdiags([[ones(p-1, 1); 0] -2*ones(p, 1) [0; ones(p-1,1)]], -1:1, p, p);
    tic
        options = odeset('JPattern',jpat(N(i)),'RelTol',1e-3);
        [t,u] = ode23s(@(t,u) f_ut(t,u,N(i),A), timespan, u0, options);
        clockJPattern(i) = toc;
        ntstpsJPattern(i) = length(t);
        tmaxJPattern(i) = max(diff(t));
    %% find indeices when t=1
    lenu = length(u(:,1));
    idx = find(t>0.99 & t<1.01);
    [t0, n] = min(abs(t(idx)-ones(length(t(idx)),1)));
    idx = idx(n);
    utemp = [ones(1, idx) zeros(1, lenu-idx)]';
    u=[utemp u];
    %% surf => JPattern
    if N(i)==10
        figure;
        set(gcf,'units','normalized','outerposition',[0 0 0.4 0.7]);
        surf(y,t,u);
        title(['ode23s, N = ' num2str(N) ' , JPattern'])
        ylabel('\tau','Fontsize',15);
        xlabel('\xi','Fontsize',15);
        zlabel('u(\tau,\xi)','Fontsize',15)
        ax = gca;
        ax.FontSize = 15; 
        grid on; box on;
        saveas(gcf,['ex3_5_1.pdf']);
    end
    %% Jacobiabn test
    %% set JPattern in options
    jpat = @(p) spdiags([[ones(p-1, 1); 0] -2*ones(p, 1) [0; ones(p-1,1)]], -1:1, p, p);
    tic
        options = odeset('Jacobian',jpat(N(i)),'RelTol',1e-3);
        [t,u] = ode23s(@(t,u) f_ut(t,u,N(i),A), timespan, u0, options);
        clockJacobian(i) = toc;
        ntstpsJacobian(i) = length(t);
        tmaxJacobian(i) = max(diff(t));
    %% find index of t=1
    lenu = length(u(:,1));
    idx = find(t>0.99 & t<1.01);
    [t0, n] = min(abs(t(idx)-ones(length(t(idx)),1)));
    idx=idx(n);
    utemp = [ones(1, idx) zeros(1, lenu-idx)]';
    u=[utemp u];
    %% surf =>  Jacobian
    if N(i)==10
        figure;
        set(gcf,'units','normalized','outerposition',[0 0 0.4 0.7]);
        surf(y,t,u)
        title(['ode23s, N=' num2str(N) ' , Jacobian'])
        ylabel('\tau','Fontsize',15);   
        xlabel('\xi','Fontsize',15); 
        zlabel('u(\tau,\xi)','Fontsize',15)
        ax = gca;
        ax.FontSize = 15; 
        grid on; box on;
        saveas(gcf,['ex3_5_2.pdf']);
    end
end
tstpsJacobian = ntstpsJacobian'
cputimeJacobian = clockJacobian'
dtmaxJacobian = tmaxJacobian'
tstpsJPattern = ntstpsJPattern'
cputimeJPattern = clockJPattern'
dtmaxJPattern = tmaxJPattern'
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
function ut=f_ut(t,u,Ny,A)
    %% du/dt=Au+b(t)
    hy = 1/Ny;
    b = zeros(Ny,1);
    b(1,1) = alp(t);
    ut = (1/hy^2)*(A*u+b);
end


