%% SF2520 HT19-1 Applied Numerical Methods
%% Computer exercise 1
%% Nomerical solutionof Initial value probles
%% Part 3: Parameter study of solutions of an ODE-system
%% Problem 2: Motion of a particle (ballistics)
%% Armin Salmasi

clear variables; close all; clc;
diary output132.txt 

%% initialization
tspan=[0 15]; % timespan
k=[0.020 0.065]; % drag coefficients
y0=1.5; % initial height
x0 = 0; % initial distance
v0=20; % initial velocity
ang=[30 45 60]; % ballistics angles in degree
tg = [];
%% calcualtion
for i= 1:length(k) %loop over drag coefficients
    figure
    hold on
    for j= 1:length(ang) % loop over angles 
        %% Ode solver IVP ode23 
        it=[x0 y0 v0*cosd(ang(j)) v0*sind(ang(j))]; % initial coordinates and velocities
        opts=odeset('RelTol',1e-6, 'Stats', 'on');
        f = @(t,x)... % right handside function of the ODE
            [x(3); x(4); -k(i)*x(3)*sqrt(x(3)^2+x(4)^2); ...
            -9.81-k(i)*x(4)*sqrt(x(3)^2+x(4)^2);]; 
        [t,x] = ode45(f, tspan, it , opts); % t,x = time and results 
        disp('%%%%%%%%%%%%%%');
        %% cut negetive values from the time particle hit the ground 
        indices = find((x(:,2))<0); % check if y<0
        x(indices,:) = []; %cut negetive values 
        tg = [tg min(indices)]; % time when particle is on the ground
        %% plot Particle Trajectory ( y vs x )
        plot(x(:,1),x(:,2),'LineWidth',1)
        set(gcf,'units','normalized','outerposition',[0 0 0.4 0.7])
        xlabel('x','Fontsize',15);
        ylabel('y','Fontsize',15);
        set(gca,'FontSize',15);
        grid on; box on;
        leg(j)=({['angele = ',num2str(ang(j))]});
    end
    title(['drag coefficients = ' num2str(k(i))])
    legend(leg)
    saveas(gcf,['fig_ex132_' num2str(i) '.pdf']);

end
tg
diary off





















