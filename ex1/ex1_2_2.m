%% SF2520 HT19-1 Applied Numerical Methods
%% Computer exercise 1 
%% Nomerical solutionof Initial value probles
%% Part 2: Stability investigation of a Runge-Kutta method
%% Problem 2: Adaptive stepsize experiment using Matlab functions
%% Armin Salmasi 

clear variables; close all; clc;
diary output12.txt 

%% initializing
rt=[1e-3 1e-4 1e-5 1e-6];   % relative tolerances
k1 = 0.04;  % initial values of the constant k1..k3
k2 = 1e4;   % initial values of the constant k1..k3
k3 = 3e7;   % initial values of the constant k1..k3
x0 = [1 0 0];   % initial values of x1..x3
tspan_ode23 = [0 1];    % time limits for ode23
tspan_ode23s = [0 1000];    % time limits for ode23s

%% IVP solver
%% ode23 non stiff solver
for i = 1:length(rt)
    disp(["ode23: for relative tolerance " ,num2str(rt(i))])
    %% solve for x1..x3 with ode23,results are [x,t]
    opts = odeset('RelTol',rt(i),'Stats','on'); % ode23 options
    [t,x] = ode23(@(t,x)... %solver  
           [-k1*x(1)+k2*x(2)*x(3); ...
           k1*x(1)-k2*x(2)*x(3)-k3*x(2)^2; k3*x(2)^2],...
           tspan_ode23,x0,opts);      
    disp('####################')
    %% plot stepsize vs time for RelTol 10^-6
    if rt(i) == 1e-6 
        figure('units','normalized','outerposition',[0 0 0.4 0.7]);
        plot(t(1:end-1),diff(t),'r','LineWidth',1);
        legend({['relative tolerances: ',num2str(rt(i))]});
        title('ode23');
        xlabel('time','Fontsize',15);
        ylabel('stepsize','Fontsize',15);
        set(gca,'FontSize',15);
        grid on; box on;
        %saveas(gcf,'fig_ex122_ode23_1.pdf');
    end
end
%% plot x1,x2 and x3 vs time - ode23
figure('units','normalized','outerposition',[0 0 0.4 0.7])
plot(t,x(:,1),t,x(:,2),t,x(:,3),'LineWidth',2);
title('ode23');
xlabel('time','Fontsize',15);
ylabel('x','Fontsize',15);
legend('x_1','x_2','x_3');
set(gca,'FontSize',15);
grid on; box on;
%saveas(gcf,'fig_ex122_ode23_2.pdf');

%% IVP solver
%% ode23s stiff solver
for i = 1:length(rt)
    disp(["ode23s: for relative tolerance " ,num2str(rt(i))])
    %% solve for x1..x3 with ode23,results are [x,t]
    opts = odeset('RelTol',rt(i),'Stats','on');  % ode23s options
    [t,x] = ode23s(@(t,x)...% solver
           [-k1*x(1)+k2*x(2)*x(3); ...
           k1*x(1)-k2*x(2)*x(3)-k3*x(2)^2; k3*x(2)^2],...
           tspan_ode23s,x0,opts);      
    disp('####################')
    %% plot stepsize vs time for RelTol 10^-6
    if rt(i) == 1e-6 
        figure('units','normalized','outerposition',[0 0 0.4 0.7]);
        plot(t(1:end-1),diff(t),'r','LineWidth',2);
        legend({['relative tolerances: ',num2str(rt(i))]});
        title('ode23s');
        xlabel('time','Fontsize',15);
        ylabel('stepsize','Fontsize',15);
        set(gca,'FontSize',15);
        grid on; box on;
        %saveas(gcf,'fig_ex122_ode23s_1.pdf');
    end
end
%% plot for x1,x2 and x3 as a function of time - ode23s
figure('units','normalized','outerposition',[0 0 0.4 0.7])
plot(t,x(:,1),t,x(:,2),t,x(:,3),'LineWidth',2);  
title('ode23s');
xlabel('time','Fontsize',15);
ylabel('x','Fontsize',15);
legend('x_1','x_2','x_3');
set(gca,'FontSize',15);
grid on; box on;
%saveas(gcf,'fig_ex122_ode23s_2.pdf');
diary off