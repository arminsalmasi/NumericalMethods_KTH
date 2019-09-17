%% SF2520 HT19-1 Applied Numerical Methods
%% Computer exercise 1
%% Nomerical solutionof Initial value probles
%% Part 3: Parameter study of solutions of an ODE-system
%% Problem 1: Particle flow past a cylinder
%% Armin Salmasi

clear variables; close all; clc;
diary output131.txt

%% intialization
y0 = [ -1.0 -1.6 0.2 0.6 1.0 1.6 3]; % vertical positions of particles
x0 = -4;
r = 2; % radius of cylinder
ts = 0; % start time
te = 10; % end time
N=50; % number of steps
h=te/N; % timestep
t=ts:h:te; %time span

%% calculations
for i = 1:length(y0) % loop over particle positions
    x=[x0; y0(i)];   % x y coorinates of the particle    
    f= @(t,x)...    % right handside function of the ODE    
        [1-r^2*(x(1)^2-x(2)^2)/(x(1)^2+x(2)^2)^2; ...
        -2*x(1)*x(2)*r^2/(x(1)^2+x(2)^2)^2];
    %% Runge-Kutta Method 4k
    for tstp=1:length(t)-1 % loop over time
        %             k1=f(t(k),x(:,k));
        %             k2=f(t(k)+h, x(:,k)+h.*k1);
        %             k3=f(t(k)+h/2, x(:,k)+h.*k1/4+h.*k2/4);
        %             x(:,k+1)=x(:,k)+h/6.*(k1+k2+4.*k3);
        k1 = h * f(t(tstp), x(:,tstp));
        k2 = h * f(t(tstp)+0.5*h, x(:,tstp)+ 0.5*k1);
        k3 = h * f(t(tstp)+0.5*h, x(:,tstp)+ 0.5*k2);
        k4 = h * f(t(tstp)+0.5*h, x(:,tstp)+ k3);
        x(:,tstp+1)=x(:,tstp)+(k1+2*k2+2*k3+k4)/6;
    end
    %% plot x y coordinates of particles
    figure(1);hold on;
    set(gcf,'units','normalized','outerposition',[0 0 0.4 0.7])
    plot(x(1,:),x(2,:),'LineWidth',1);
    xlabel('x','Fontsize',15);
    ylabel('y','Fontsize',15);
    set(gca,'FontSize',15);
    grid on; box on;
    leg(i)=({['y = ',num2str(y0(i))]});
end
ang = linspace(0,360,721);
%% plot cylinder
plot(r*cosd(ang),r*sind(ang),'k-','LineWidth',1)
leg(i+1)=({['Cylinder, r = ', num2str(r)]});
grid on;
legend(leg);
%saveas(gcf,'fig_ex131.pdf');
    axis equal

diary off