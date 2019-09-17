%% SF2520 HT19-1 Applied Numerical Methods
%% Computer exercise 1 
%% Nomerical solutionof Initial value probles
%% Part1 : Accuracy of Runge-Kutta method
%% Armin Salmasi 

clear variables; close all; clc;

%% initialization
f = @(u) [u(2);-u(1)-u(2)*(u(1)^2-1)];  % right handside function of the ODE    
N = [10 20 40 80 160 320];  % number of steps
h = 1./N;   % timestep
ts = 0; % start time
te = 1; % end time
yN = N*0;   % results at t=1 for different N

%% calcualtion
for i = 1:length(N) % loop over step sizes
    t = ts:h(i):te; % time span
    u = zeros(2,length(t)); % results at each timestep 
    u(1:2,1) = [1; 0];  % initial values at time zero
    %% Runge-Kutta method     
    for j=1:(length(t)-1)   % loop over timespan                          
        k_1 = f(u(:,j)); 
        k_2 = f(u(:,j)+h(i).*k_1);
        k_3 = f(u(:,j)+0.25.*h(i).*k_1+0.25.*h(i).*k_2);
        u(:,j+1) = u(:,j) + (h(i)/6).*(k_1+k_2+4.*k_3);
    end    
    yN(i) = u(1,length(t)) ; % result at t=1 for different N
end
%% error evaluation / calculation of hte slop
eN = abs(yN - yN(length(N))) 
slope = log(eN(1)/eN(5))/log(h(1)/h(5)); % Calculate the slope

%% plot loglog timestep vs error
figure('units','normalized','outerposition',[0 0 0.4 0.7])
loglog(h,eN , "r-o",'LineWidth',2)
grid on
box on
xlabel('log(timestep)','Fontsize',15)
ylabel('log(error)','Fontsize',15)
set(gca,'FontSize',15)
an = annotation('textbox', [.48 .35 .1 .1], 'String', ...
    ['slope at t=1 is ',num2str(slope)]); an.FontSize = 15;
%saveas(gcf,'fig_ex11.pdf')