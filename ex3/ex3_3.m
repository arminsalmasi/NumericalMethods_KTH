
%% SF2520 HT19-1 Applied Numerical Methods
%% Computer exercise 3 part 3 Solution and visualization
%% Partial differential equation of parabolic type
%% Armin Salmasi
%% Stability of Euler’s explicit method.
%% Experiment with different values of the discretization step space and time and study of stability.
close all ;clear all; clc;
%% initialization
%% discritization of time and space
Nx=5;     % spatial discritization, number of grids
hx=1/Nx
y=0:hx:1; 
%%%%%%%%%%%%%%%%%%
Nt=50;    % time discritization, number of timesteps
ht=1/Nt
t=0:ht:2;
stability=ht/hx^2
% Nx=5 Nt=50   stable
% Nx=5 Nt=40   unstable
%% A matrix %% Sparse matrix
% diagonals for matrix A:
diag1 =ones(Nx-1,1);
diag2 =-ones(Nx,1)*2;
diag3 =ones(Nx-1,1);
A = sparse((diag(diag1,-1)+diag(diag2)+diag(diag3,1)));
%% B matrix
B = sparse(zeros(Nx,1));
%% U and u matrices
u0 = sparse(zeros(Nx,1)); % initial solution
U=[];   % matrix U of all solutions
uk = u0;  % solution after each iteration
%% Boundary condition:
A(Nx,Nx-1)=2;
A(Nx,Nx)=-2;
%% T matrix
T = (1/hx^2).*A;
%% Euler's explicit method
for ti = t
    U=[U uk];
    B(1)=alp(ti);
    uk_tmp=uk+ht*(T*uk+(1/hx^2).*B);
    uk=uk_tmp;
end

%% include initial and boundary conditions 
%% indices of t=1
idx = find(t==1);
u0x = [ones(1,idx) zeros(1,length(t)-idx)];
U = [u0x; U];
%% Post processing 
figure;
    set(gcf,'units','normalized','outerposition',[0 0 0.4 0.7]);
    surf(t,y,U);
    xlabel('\tau','Fontsize',20);
    ylabel('\xi','Fontsize',20);
    zlabel('u(\tau,\xi)','Fontsize',20);
    title(['h_{\xi} = ' num2str(hx) ', h_{\tau} = ' num2str(ht)]);
    ax = gca;
    ax.FontSize = 20; 
    grid on; box on;
    saveas(gcf,['ex3_3_1.pdf']);
%% find solutions at t=[0.5,1,1.5,2]:
ts = [0.5 1 1.5 2];
for i = 1:length(ts)
    idx2 = find(t==ts(i));
    u_tmp(i,:) = U(:,idx2);
end
figure;
    hold on;
    set(gcf,'units','normalized','outerposition',[0 0 0.4 0.7]);
    for i = 1:length(ts)
        plot(y,u_tmp(i,:),'LineWidth',2);
    end
    xlabel('\xi','Fontsize',20);   
    ylabel('u(\xi, \tau)','Fontsize',20);
    legend('\tau = 0.5','\tau = 1','\tau = 1.5','\tau = 2','Fontsize',18);
    ax = gca;
    ax.FontSize = 20; 
    grid on; box on;
    saveas(gcf,['ex3_3_2.pdf']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% function declrations
function a = alp(t)
%% time dependent left handside boundary condition
if t > 1
    a = 0;
else
    a = 1;
end
end
