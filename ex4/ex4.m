%% SF2520 HT19-1 Applied Numerical Methods
%% Computer exercise 4
%% Numerical Solution of a elliptic PDE problem
%% Armin Salmasi 

%% 1-a
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear variables;clc;
% domain
Lx = 5;
Ly = 2;
h = 0.1;
M = Ly/h+1;   % number of points on y axis
N = Lx/h-1;   % number of points on x axis
y = linspace(0,Ly,M);  % y grid
x = linspace(0,Lx,N);  % x grid
%% initial values
f = 0;  % source term
b = zeros(M*N,1)+f*h^2;
%% A matrix
SN = full(gallery('tridiag',N,-1,2,-1));
SM = full(gallery('tridiag',M,-1,2,-1));
A = -(kron(-eye(M),SN)+kron(SM,-eye(N)));
%% Neumann BC
for j = 1:N
    A(j,N+j) = A(j,N+j)*2;
end
for j = M*N:-1:M*N-N+1
    A(j,j-N) = A(j,j-N)*2;
end
%% Dirichlet BC
DirLeft = 20;
DirRight = 100;
b(1:N:end) = b(1:N:end)+DirLeft;
b(N:N:end) = b(N:N:end)+DirRight;
%% Solver
C = A\b;
%% Plotting
figure
CC = reshape(C,N,M);
[Y,X]  =  meshgrid(x,y);
surf(Y,X,CC');
view(0,90)
xlabel('X');
ylabel('Y');
title('h = 0.1, f = 0');
c = colorbar;
c.Label.String  =  'T K';
%% find T(2.5,1)
xi = find(x == 2.5);
yi = find(y == 1);
['h = ' num2str(h), ', T(2.5,1) = ', num2str(CC(xi,yi))]
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1-c
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
% domain
Lx = 5;
Ly = 2;
H = [0.1,0.05,0.025];
for h  =  H
    M = Ly/h+1;   % number of points on y axis
    N = Lx/h-1;   % number of points on x axis
    y = linspace(0,Ly,M);  % y grid
    x = linspace(0,Lx,N);  % x grid
    %%  source term
    f = @(x,y)(50+400*exp(-(x-1)^2-2*(y-1.5)^2));
    k = 1;
    for i = 1:M
        for j = 1:N
            b(k,1) = f(x(1,j),y(1,i))*h^2;
            k = k+1;
        end
    end
    %% A matrix
    SN = (gallery('tridiag',N,-1,2,-1));
    SM = (gallery('tridiag',M,-1,2,-1));
    A = -(kron(-eye(M),SN)+kron(SM,-eye(N)));
    %% Neumann BC
    for j = 1:N
        A(j,N+j) = A(j,N+j)*2;
    end
    for j = M*N:-1:M*N-N+1
        A(j,j-N) = A(j,j-N)*2;
    end
    %% Dirichlet BC
    DirLeft = 20;
    DirRight = 100;
    b(1:N:end) = b(1:N:end)+DirLeft;
    b(N:N:end) = b(N:N:end)+DirRight;
    %% Solver
    C = A\b;
    %% Plotting
    CC = reshape(C,N,M);  
    if h == 0.05
        [Y,X]  =  meshgrid(x,y);
        figure;
        mesh(Y,X,CC');
        xlabel('X');
        ylabel('Y');
        zlabel('T');
        title('h = 0.05, f = 50+400exp(-(x-1)^2-2(y-1.5)^2)')
        c = colorbar;
        c.Label.String  =  'T K';
        figure;
        contour(Y,X,CC',20);
        xlabel('X');
        ylabel('Y');
        zlabel('T');
        title('h = 0.05, f = 50+400exp(-(x-1)^2-2(y-1.5)^2)');
        c = colorbar;
        c.Label.String  =  'T K';
        figure;
        meshc(Y,X,CC');
        xlabel('X');
        ylabel('Y');
        zlabel('T');
        title('h = 0.05, f = 50+400exp(-(x-1)^2-2(y-1.5)^2)')
        c = colorbar;
        c.Label.String  =  'T K';
    end
    %% find T(2.5,1)
    xi = find(x == 2.5);
    yi = find(y == 1);
    ['h = ' num2str(h), ', T(2.5,1) = ', num2str(CC(xi,yi))]
end