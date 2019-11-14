%% SF2520 HT19-2 Applied Numerical Methods
%% Computer exercise 5
%% Numerical experiments with hyperbolic PDE problems
%% Armin Salmasi 
%% Part 1
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear variables;clc;
g_sin=@(t,T) sin(2*pi*t/T); 
g_sq=@(t,T) square(t*2*T*pi);
%% disceritization
N = 100; % number of z grids
T = 1; % end time
a = 1; % a coefficient
Lz=2; % domain size
h = Lz/N; % step in z (dz) 
sigm = 0.9; %[0.1 0.5 0.9 1 1.1]% sigma =a*dt/h
dt = sigm*h/a; % calculate dt from sigma
ts = dt:dt:2*T+dt; % time disceritization
zs = h:h:Lz;  % z grid disceritization 
figure(1)
for funi=[1,2] % i=1 square function/ i=2 sin function
    % A matrix
    A_LxF = full(gallery('tridiag',N,0.5*(1+sigm),0,0.5*(1-sigm)));
    if a>0
        A_UPW = full(gallery('tridiag',N,sigm,1-sigm,0));
    else
        A_UPW = full(gallery('tridiag',N,0,1-sigm,-sigm));
    end    
    A_LxW = full(gallery('tridiag',N,0.5*(sigm+sigm^2),1-sigm^2,0.5*(sigm^2-sigm)));
    % initial value
    u_LxF(:,1) = zeros(N,1);
    u_UPW(:,1) = zeros(N,1);
    u_LxW(:,1) = zeros(N,1);
    b_LxF=u_LxF;
    b_UPW=u_UPW;
    b_LxW=u_LxW;
    % itteration over time 
    for t =1:length(ts)
        % left BC 
        if funi==1
            b_LxF(1)=0.5*(1+sigm)*g_sq(ts(t),T);
            if a>0
                b_UPW(1) = sigm*g_sq(ts(t),T);
            else
                b_UPW(end) = sigm*g_sq(ts(t),T);
            end    
            b_LxW(1) = 0.5*sigm*(sigm+1)*g_sq(ts(t),T);
        else
            b_LxF(1)=0.5*(1+sigm)*g_sin(ts(t),T);
            if a>0
                b_UPW(1) = sigm*g_sin(ts(t),T);
            else
                b_UPW(end) = sigm*g_sin(ts(t),T);
            end    
            b_LxW(1) = 0.5*sigm*(sigm+1)*g_sin(ts(t),T);
        end
        % calculate the solution
        u_LxF(:,t+1)=A_LxF*u_LxF(:,t)+b_LxF;
        u_UPW(:,t+1)=A_UPW*u_UPW(:,t)+b_UPW;
        u_LxW(:,t+1)=A_LxW*u_LxW(:,t)+b_LxW;
        % right bc extrapolation for LxW and LxF
        if t == length(ts)
           u_LxF(end,t+1) =  2*u_LxF(end-1,t+1) - u_LxF(end-2,t+1);
           u_LxW(end,t+1) =  2*u_LxW(end-1,t+1) - u_LxW(end-2,t+1);
        end
    end
    % Plotting
    subplot(1,2,funi)
        hold on;
        box on;
        grid on;
        plot(zs,u_LxW(:,end));
        plot(zs,u_UPW(:,end));
        plot(zs,u_LxF(:,end));
        xlabel('x');
        ylabel('u');
        legend({'lax wendroff' ,'upwind', 'lax friedrich'});
    clear u_LxF u_LxW u_UPW ;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
presskey()
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SF2520 HT19-2 Applied Numerical Methods
% Computer exercise 5
% Numerical experiments with hyperbolic PDE problems
% Armin Salmasi
% part 2
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables;clc;
U_0t=@(t,Tc,Th)((t<=0.25&&0<=t)*(Tc+(Th-Tc)*sin(2*pi*t))+ ...
    (t<=2&&0.25<t)*Th + ...
    (Th+Tc*sin(4*pi*(t-2)))*(t>2));   
% disceritization
h = 0.1; % step in z (dz)
dt = 0.05; % step in t
v = 1; % velocity of the fluid
Lz = 3; % domain size
T=5;
ts = dt:dt:2*T+dt; % time disceritization
zs = 0:h:Lz;  % z grid disceritization
k = 0.2 ;%heat exchange parameter
uc= 10; % T cool
uh = 100; % Thot
alp=k*dt 
gam = alp-(alp^2)/2
sigm = v*dt/h
N=length(zs);
% A matrix
if v>0
    A_UPW = full(gallery('tridiag',N,sigm,1-sigm-alp,0));
end
A_LxW = full(gallery('tridiag',N,0.5*(sigm*(1-alp)+sigm^2),...
                                1-sigm^2-gam,...
                                0.5*(-sigm*(1-alp)+sigm^2)));
% initial value
u_UPW(:,1) = ones(N,1)*uc;
u_LxW(:,1) = ones(N,1)*uc;
b_UPW=ones(N,1)*alp;
b_LxW=ones(N,1)*gam*uc;
% itteration over time
for t =1:length(ts)
    % left BC
    if v>0
        b_UPW(1) = b_UPW(1)+sigm*U_0t(t,uc,uh);
    end
    b_LxW(1)=(b_LxW(1)+ 0.5*(sigm*(1-alp)+sigm^2)*U_0t(t,uc,uh));
    % calculate the solution
    u_UPW(:,t+1)=A_UPW*u_UPW(:,t)+b_UPW;
    u_LxW(:,t+1)=A_LxW*u_LxW(:,t)+b_LxW;
    % right bc extrapolation for LxW
    if t == length(ts)
        u_LxW(end,t+1) =  2*u_LxW(end-1,t+1) - u_LxW(end-2,t+1);
    end
end
%% Plotting
% part a
figure(2)
subplot(2,2,1)
box on;
grid on;
surf([0 ts],zs,u_UPW);
title('upwind')
xlabel('t')
ylabel('x')
zlabel('u')
subplot(2,2,2)
box on;
grid on;
surf([0 ts],zs,u_LxW);
title('lax friedrich')
xlabel('t')
ylabel('x')
zlabel('u')
% part b
% find values ot t=2.5 and t=5
idx1=find(ts==2.5);
idx2=find(ts==5);
subplot(2,2,3)
hold on
box on;
grid on;
plot(zs,u_UPW(:,idx1));
plot(zs,u_LxW(:,idx1));
title('t=2.5')
xlabel('x')
ylabel('u')
legend({'upwind','lax friedrich'})
subplot(2,2,4)
hold on
box on;
grid on;
plot(zs,u_UPW(:,idx2));
plot(zs,u_LxW(:,idx2));
title('t=5')
xlabel('x')
ylabel('u')
legend({'upwind','lax friedrich'})

%% %%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=presskey()
%% wait for enter key
    disp("press enter to continue to part 2")
    currkey=0;
    while currkey~=1
        pause; % wait for a keypress
        currkey=get(gcf,'CurrentKey'); 
        if strcmp(currkey, 'return') % You also want to use strcmp here.
            currkey=1 % Error was here; the "==" should be "="
        else
            currkey=0 % Error was here; the "==" should be "="
        end
    end
end