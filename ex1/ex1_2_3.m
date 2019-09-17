clear all, close all, clc

tspan=[0 1000]; %time span
x=[1 0 0]; %Initial values of the variables x1, x2 and x3

%Initial values of the constants
k1=0.04;
k2=1e4;
k3=3e7;

R=[1e-3 1e-4 1e-5 1e-6]; %Values of the relative tolerances
figure 
hold on 
for i=1:length(R)
    options=odeset('RelTol',R(i), 'Stats', 'on'); %options for ode23s
    [t,x] = ode23s(@(t,x) [-k1*x(1)+k2*x(2)*x(3); k1*x(1)-k2*x(2)*x(3)- ...
        k3*x(2)^2; k3*x(2)^2], tspan, [1 0 0], options);                                   %solving for x1, ...x2 and x3 using ode23s
    
    disp('ode23s stats:') %display ode23 status on screen with successful steps taken
    
    
    %figure(i)
    plot(t(1:end-1), diff(t))
    %     title(['Stepsize (h) as a function of time (t) for RelTol=' num2str(R(i))], 'Fontsize',24)
    xlabel('time(t)', 'Fontsize',16)
    ylabel('stepsize(h)', 'Fontsize',16)
    grid on
    grid minor
    %     tic, ode23s(@(t,x) [-k1*x(1)+k2*x(2)*x(3); k1*x(1)-k2*x(2)*x(3)-k3*x(2)^2; k3*x(2)^2], tspan, [1 0 0], options),toc
end

figure(i+2)
plot(t,x(:,1),t,x(:,2),t,x(:,3));                               %plot for x1, x2 and x3 as a function of time
xlabel('time(t)', 'Fontsize',16)
ylabel('x', 'Fontsize',16)
legend('x1','x2','x3')

