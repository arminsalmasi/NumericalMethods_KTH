%% SF2520 HT19-1 Applied Numerical Methods
%% Computer exercise 2 
%% Numerical Solution of a Boundary Value Problem
%% Armin Salmasi 

clear variables; close all; clc;

%% Part1: v=0 with N = 10,20,40,80 in. Is it second order?
L = 10;   % L domain size (pipe size)
aso = 1;   % a start of the source domain  
bso = 3;     % b end of the source domain
Q0_src = 50;     % amplitude of the source 
kp = 0.5;         % heat conduction coefficient
alpha0 = 10;        % heat sink coefficient
ro = 1;            %fluid density
cf = 1;        % heat capacity of the fluid
Tout = 300;         % heat sink constant 
T0 = 400;           % inlet temperature
v = 0; % fluid velocity 
N = [10, 20, 40, 80]; % number of timesteps to be taken
figure 
hold on
%% loop over total number of timesteps N
for i = 1:size(N,2) 
    dz = L / N(i);  % %lenght of the finit stepsize in spacial domain
    z = 0:dz:L; % descritisized spacial domain    
    av = sqrt((v^2*ro^2*cf^2)/4 + alpha0^2)...
        -v*ro*cf/2; % heat sink quefficient 
    %% quefficient (A) matrix terms
    c = -kp-v*ro*cf*dz/2;  
    d = 2*kp;
    e = -kp+v*ro*cf*dz/2;
    E_u = e * diag(ones(size(z,2)-2,1),1); % upper diagonal term
    E_m = d * eye(size(z,2)-1); %diagonal term
    E_d=  c * diag(ones(size(z,2)-2,1),-1); % lower diagonal term
    A = E_u+E_d+E_m; % quefficient matrix
    %% Robin conditions, last line of the quefficietnt matrix
    A(end,end) = d-e*dz*av/kp;
    A(end,end-1) = c+e;
    %% boundary conditions + heat source terms 
    B =dz^2.*(q_calc(Q0_src,aso,bso,z(2:end)))'; % source term
    B(1,1) = B(1,1)-c*T0; % dirichlet condition
    B(end,1) = B(end,1)-dz*e*av*Tout/kp; % Robin condition   %% Or 2*dz*e*av*Tout/kp ?? 
    %% solve the linear system of equations 
    RES = double(A\B);
    %% Plot the results
    plot(z(1:end),[T0 RES(1:end)']);
    leg(i) = {['N=' num2str(N(i)) ', dt=' num2str(dz)]};
end
legend(leg);
title(['v=' num2str(v) ]);
xlabel('distance');
ylabel('temperature');
grid on; box on;

%% Part2: v=.01,0.5,1,10 with N = 40 in. Is it second order?
v = [0.1,0.5,1,10]; % fluid velocity 
N = 40; % number of timesteps to be taken
clear leg 
figure 
hold on
dz = L/N;  % %lenght of the finit stepsize in spacial domain
z = 0:dz:L; % descritisized spacial domain   
%% loop over velocity of liquid v
for i = 1:size(v,2) % loop over total number of timesteps
    av = sqrt((v(i)^2*ro^2*cf^2)/4 + alpha0^2)...
        -v(i)*ro*cf/2; % heat sink quefficient 
    %% quefficient (A) matrix terms
    c = -kp-v(i)*ro*cf*dz/2;  
    d = 2*kp;
    e = -kp+v(i)*ro*cf*dz/2;
    E_u = e * diag(ones(size(z,2)-2,1),1); % upper diagonal term
    E_m = d * eye(size(z,2)-1); %diagonal term
    E_d=  c * diag(ones(size(z,2)-2,1),-1); % lower diagonal term
    A = E_u+E_d+E_m; % quefficient matrix
    %% Robin conditions, last line of the quefficietnt matrix
    A(end,end) = d-e*dz*av/kp;
    A(end,end-1) = c+e;
    %% boundary conditions + heat source terms 
    B =dz^2.*(q_calc(Q0_src,aso,bso,z(2:end)))'; % source term
    B(1,1) = B(1,1)-c*T0; % dirichlet condition
    B(end,1) = B(end,1)-dz*e*av*Tout/kp; % Robin condition   %% Or 2*dz*e*av*Tout/kp ?? 
    %% solve the linear system of equations 
    RES = double(A\B);
    %% Plot the results
    plot(z(1:end),[T0 RES(1:end)']);
    leg(i) = {['v=' num2str(v(i))]};
end
legend(leg);
title(['N=' num2str(N) ' dz=' num2str(dz) ]);
xlabel('distance');
ylabel('temperature');
grid on; box on;

%% Part3: v=10 with N=[10,20,40].Pronaounced oscilation for higher h.
v = 10; % fluid velocity 
N = [10,20,40]; % number of timesteps to be taken
clear leg
figure 
hold on
%% loop over total number of timesteps
for i = 1:size(N,2);
    dz = L / N(i);  % %lenght of the finit stepsize in spacial domain
    z = 0:dz:L; % descritisized spacial domain    
    av = sqrt((v^2*ro^2*cf^2)/4 + alpha0^2)...
        -v*ro*cf/2; % heat sink quefficient 
    %% quefficient (A) matrix terms
    c = -kp-v*ro*cf*dz/2;  
    d = 2*kp;
    e = -kp+v*ro*cf*dz/2;
    E_u = e * diag(ones(size(z,2)-2,1),1); % upper diagonal term
    E_m = d * eye(size(z,2)-1); %diagonal term
    E_d=  c * diag(ones(size(z,2)-2,1),-1); % lower diagonal term
    A = E_u+E_d+E_m; % quefficient matrix
    %% Robin conditions, last line of the quefficietnt matrix
    A(end,end) = d-e*dz*av/kp;
    A(end,end-1) = c+e;
    %% boundary conditions + heat source terms   
    B =dz^2.*(q_calc(Q0_src,aso,bso,z(2:end)))'; % source term
    B(1,1) = B(1,1)-c*T0; % dirichlet condition
    B(end,1) = B(end,1)-dz*e*av*Tout/kp; % Robin condition   %% Or 2*dz*e*av*Tout/kp ?? 
    %% solve the linear system of equations 
    RES = double(A\B);
    %% Plot the results
    plot(z(1:end),[T0 RES(1:end)']);
    leg(i) = {['N=' num2str(N(i)) ', dt=' num2str(dz)]};
end
legend(leg);
title(['v=' num2str(v) ]);
xlabel('distance');
ylabel('temperature');
grid on; box on;

%% Part4: v=10 with N=[10,20,40]. Upwind difference 
v = 10; % fluid velocity 
N = [10,20,40]; % number of timesteps to be taken
clear leg
figure 
hold on
%% loop over total number of timesteps
for i = 1:size(N,2)
    dz = L / N(i);  % %lenght of the finit stepsize in spacial domain
    z = 0:dz:L; % descritisized spacial domain    
    av = sqrt((v^2*ro^2*cf^2)/4 + alpha0^2)...
        -v*ro*cf/2; % heat sink quefficient 
    %% quefficient (A) matrix terms
    c = -kp-v*ro*cf*dz;  
    d = 2*kp+v*ro*cf*dz;
    e = -kp;
    E_u = e * diag(ones(size(z,2)-2,1),1); % upper diagonal term
    E_m = d * eye(size(z,2)-1); %diagonal term
    E_d=  c * diag(ones(size(z,2)-2,1),-1); % lower diagonal term
    A = E_u+E_d+E_m; % quefficient matrix
    %% Robin conditions, last line of the quefficietnt matrix
    A(end,end) = d-e*dz*av/kp;
    A(end,end-1) = c+e;
    %% boundary conditions + heat source terms 
    B =dz^2.*(q_calc(Q0_src,aso,bso,z(2:end)))'; % source term
    B(1,1) = B(1,1)-c*T0; % dirichlet condition
    B(end,1) = B(end,1)-dz*e*av*Tout/kp; % Robin condition   %% Or 2*dz*e*av*Tout/kp ?? 
    %% solve the linear system of equations 
    RES = double(A\B);
    %% Plot the results
    plot(z(1:end),[T0 RES(1:end)']);
    leg(i) = {['N=' num2str(N(i)) ', dt=' num2str(dz)]};
end
legend(leg);
title(['v=' num2str(v) ', Upwind difference' ]);
xlabel('distance');
ylabel('temperature');
grid on; box on;


function out = q_calc(A0, min_value, max_value, list)
count=1;
for  element = list
    if  (element>=min_value) && (element<=max_value)
        out(count)= A0*sin((element-min_value)*pi/(max_value-min_value));
        count=count+1;
    else
        out(count) = 0;
        count=count+1;
    end
end
end