%% SF2520 HT19-1 Applied Numerical Methods
%% Computer exercise 1 
%% Nomerical solutionof Initial value probles
%% Part 2: Stability investigation of a Runge-Kutta method
%% Problem 1: Constant stepsize experiment
%% Armin Salmasi 

clear variables;close all;clc;

%% initializing
N = [125 250 500 1000 2000]; % number of steps
k1 = 0.04; %value constant k1
k2 = 1e4; %value constant k2
k3 = 3e7; %value constant k3
ts = 0; % start time
te = 1; % end time
y1 = N*0; % x1 at t=1  ( all number of steps N) 
y2 = N*0; % x2 at t=1  ( all number of steps N) 
y3 = N*0; % x3 at t=1  ( all number of steps N) 
h_accum = N*0; % save stepsizes

%% calcualtion
for i = 1:length(N) % loop over number of steps in N
    h = te/N(i);    % calculate stepsize 
    h_accum = [h_accum h];  % accumulate stepsizes
    t = ts:h:te;    % time span
    x = zeros(3,length(t)); % ODE variables x1 x2 x3 
    x(1:3,1) = [1;0;0]; % initial values of x1,x2,x3
    f= @(t,x) ...   % right handside function of the ODE
       [-k1*x(1)+k2*x(2)*x(3); k1*x(1)-k2*x(2)*x(3)-k3*x(2)^2; k3*x(2)^2];   
    %% Runge-Kutta Method
    for k = 1:length(t)-1
        k_1 = f(t(k),x(:,k));
        k_2 = f(t(k)+h,x(:,k)+h.*k_1);
        k_3 = f(t(k)+h/2,x(:,k)+h.*k_1/4+h.*k_2/4);
        x(:,k+1) = x(:,k)+h/6.*(k_1+k_2+4.*k_3); 
    end
    y1(i) = abs(x(1,end)); % x1 at t=1 for different N
    y2(i) = abs(x(2,end)); % x2 at t=1 for different N
    y3(i) = abs(x(3,end)); % x3 at t=1 for different N
    %% plot log-log t vs x % x(1,:) is complex for N<=500
    figure('units','normalized','outerposition',[0 0 0.5 0.5])
        subplot(1,3,1)
            loglog(t,abs(x(1,:)),"r",'LineWidth',2) 
            xlabel('log(time)','Fontsize',15)
            ylabel('log(x_1)','Fontsize',15)
            set(gca,'FontSize',15)
            box on; grid on
        subplot(1,3,2)
            loglog(t,(x(2,:)),"r",'LineWidth',2)
            xlabel('log(time)','Fontsize',15)
            ylabel('log(x_2)','Fontsize',15)
            set(gca,'FontSize',15)
            box on; grid on
        subplot(1,3,3)
            loglog(t,(x(3,:)),"r",'LineWidth',2)
            xlabel('log(time)','Fontsize',15)
            ylabel('log(x_3)','Fontsize',15)
            set(gca,'FontSize',15)
            box on; grid on
        h=gcf;
        set(h,'PaperOrientation','landscape');
        %saveas(gcf,['fig_ex121_' num2str(i) '.pdf'])
end
%% error evaluation
eN1=abs(y1-y1(end)) %error of x1
eN2=abs(y2-y2(end)) %error of x2
eN3=abs(y3-y3(end)) %error of x3
