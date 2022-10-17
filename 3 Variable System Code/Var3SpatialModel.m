%% Gentrification Model: 3-Variable Spatial System: Euler's Method

%{
Authors: 
Juan G. Restrepo
Nancy Rodriguez
Jonathan D. Shaw

Date modified: 07/27/2022
%} 

%housekeeping: 
clear; close all; clc; 

warning('off')

tic

%% Constants/Initial Conditions

n = 20; 
tspan = [0 10000]; 
dt = 0.01; 
t = (0:dt:tspan(end))'; 

%struct of parameters: 
params = struct('dt', dt, 't', t, 'z', 0.01, 'tau_a', 1, 'tau_n', 50, 'tau_c', 90, 'epsilon', 0.1, 'un', ones(n,1),'t_tol', 1e-2, 'A_rel_tol', 1e-4, 'A_abs_tol', 1e-4);
                
%generating R matrix: 
R = zeros(n,n); v = ones(n-1,1); 
R = R + diag(v,1) + diag(v,-1); 
R(end,1) = 1; R(1,end) = 1;
params.R = R; 

%perturbation testing IC: 
%{
aeq = 1/n; 
delta_a = 10e-10; 

A0 = [aeq + delta_a; aeq-delta_a; aeq+delta_a/2; aeq-delta_a/4; aeq-delta_a/4]; 
N0 = sigma(A0,params.z,params.epsilon); 
C0 = N0; 
%}

%regular simulation, random IC: 
%
A0 = rand(n,1); A0 = A0./sum(A0); 
N0 = rand(n,1); 
C0 = rand(n,1);
%}

%regular simulation, nonrandom IC: 
%{
A0 = 2:2:(2*n); A0 = (A0/sum(A0))'; 
N0 = sigma(A0,params.z,params.epsilon); 
C0 = N0;
%}

%% Nearest Neighbor Coupling

%simulation:
[a_mat, n_mat, c_mat] = Euler3VarSimSpatialRing(n,A0,N0,C0,params);

%plotting artist densities:
figure
hold on
for i=1:n
plot(t,a_mat(:,i), 'linewidth', 1.3)
Legend{i} = ['Region ', num2str(i)];
end
xlabel('Time','fontsize',12,'fontweight', 'bold')
ylabel('Artist Densities','fontsize',12,'fontweight', 'bold')
title('Artist Densities vs Time','fontsize',12,'fontweight', 'bold')
legend(Legend,'fontsize', 12,'fontweight', 'bold')
grid on; grid minor
hold off

%{
%selecting segment for analysis:
t_range = 200; 
t_start_idx = floor(0.75*length(t)); 
t_seg = t(t_start_idx:(t_start_idx+t_range/dt));
a_mat_seg = a_mat(t_start_idx:(t_start_idx + t_range/dt),:); 


%plotting segment
figure
hold on
for i=1:n
plot(t_seg,a_mat_seg(:,i), 'linewidth', 1.3)
Legend{i} = ['Region ', num2str(i)];
end
xlabel('Time','fontsize',12,'fontweight', 'bold')
ylabel('Artist Densities','fontsize',12,'fontweight', 'bold')
title('Artist Densities vs Time','fontsize',12,'fontweight', 'bold')
legend(Legend,'fontsize', 12,'fontweight', 'bold')
grid on; grid minor
hold off
%}

%spatial plot: 
PlotSpatialRing(n,t,a_mat)

%getting chaos info:
GetChaosInfo(n,t,A0,N0,C0,a_mat,n_mat,c_mat,params, 'SpatialRing')

toc
