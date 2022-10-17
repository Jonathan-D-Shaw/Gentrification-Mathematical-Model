%% Gentrification Model: 3-Variable System: Euler's Method

%{
Authors: 
Juan G. Restrepo
Nancy Rodriguez
Jonathan D. Shaw

Date modified: 08/23/2022
%} 

%housekeeping: 
clear; close all; clc; 
%openfig('OscSolutionPhasePlot.fig');
warning('off')

tic

%% Constants/Initial Conditions

n = 10; 
tspan = [0 400]; 
dt = 0.01; 
t = (0:dt:tspan(end))'; 


%struct of parameters: 
params = struct('dt', dt, 't', t, 'z', 0.01, 'tau_a', 1, 'tau_n', 20, 'tau_c', 20, 'epsilon', 0.1, 'un', ones(n,1),'t_tol', 12, 'A_rel_tol', 1e-2, 'A_abs_tol', 1e-4);
                

%perturbation testing IC: 
%
delta_A = 1e-10; 
delta_N = 1e-10; 
delta_C = 1e-10; 

A0 = 1/n * ones(n,1); 
N0 = sigma(A0, params.z, params.epsilon);
C0 = N0; 

if mod(n,2) %testing for odd n
    for i = 1:n-1
        A0(i) = A0(i) + (-1)^(i-1)*delta_A;
        N0(i) = N0(i) + (-1)^(i-1)*delta_N;
        C0(i) = C0(i) + (-1)^(i-1)*delta_C;
    end
else
    for i = 1:n
        A0(i) = A0(i) + (-1)^i*delta_A;
        N0(i) = N0(i) + (-1)^i*delta_N;
        C0(i) = C0(i) + (-1)^i*delta_C;
    end
    
end
%}

%regular simulation, random IC: 
%{
A0 = rand(n,1); A0 = A0./sum(A0); 
N0 = rand(n,1); 
C0 = rand(n,1);
%}

%regular simulation, nonrandom IC: 
%{
A0 = linspace(0,1,n); A0 = (A0/sum(A0))'; 
N0 = sigma(A0,params.z,params.epsilon); 
C0 = N0;
%}

%% Euler's Method Simulation: 
%
%simulation:
[a_mat, n_mat, c_mat] = Euler3VarSim(n,A0,N0,C0,params);

%plotting:
figure
hold on
for i=1:n
plot(t,a_mat(:,i), 'linewidth', 1.3)
Legend{i} = ['Region ', num2str(i)];
end
xline(350, 'r', 'linewidth', 2)
xlabel('Time','fontsize',12,'fontweight', 'bold')
ylabel('Artist Densities','fontsize',12,'fontweight', 'bold')
title('Artist Densities vs Time','fontsize',12,'fontweight', 'bold')
legend(Legend,'fontsize', 12,'fontweight', 'bold')
grid on; grid minor
hold off


%selecting segment for analysis:
t_range = 200; 
t_start = t(end) - t_range; 
[~,t_start_idx] = min(abs(t-t_start)); 
t_seg = t(t_start_idx:length(t));
a_mat_seg = a_mat(t_start_idx:length(t),:); 

[peak_data,ind_indicators,grp_indicators, fixed_pt, base_osc, grp_osc, sep_grp_osc, chaos] = InterestingBehavior(n,t,a_mat_seg,params);

%looking at covariance matrix: 
%{
cor_mat = zeros(n,n);
for i = 1:n
    for j = 1:n
        temp_cor_mat = corrcoef(a_mat_seg(:,i),a_mat_seg(:,j));
        cor_mat(i,j) = temp_cor_mat(1,2);
    end
end
%}

%frob norm/normalized matrix: 
%{
cor_coeff = norm(cor_mat,'fro')/n;

cor_mat
cor_coeff
%}

%plotting segment
figure
hold on
for i=1:n
plot(t_seg,a_mat_seg(:,i), 'linewidth', 1.5)
Legend{i} = ['Region ', num2str(i)];
end
xlabel('Time','fontsize',12,'fontweight', 'bold')
ylabel('Artist Densities','fontsize',12,'fontweight', 'bold')
title('Artist Densities vs Time','fontsize',12,'fontweight', 'bold')
legend(Legend,'fontsize', 12,'fontweight', 'bold')
grid on; grid minor
hold off

%investigating chaos: 
%  GetChaosInfo(n,t,A0,N0,C0,a_mat,n_mat,c_mat,params, 'Non-spatial')
%}

%% 4th-Order Runge-Kutta Simulation: 
%{
x0 = [A0; N0; C0]; 
[t, xOut] = ode45(@(t,x) RungeKutta3VarSim(t,x,n,params), tspan, x0);

a_mat = xOut(:,1:n); 
n_mat = xOut(:,(n+1):(2*n)); 
c_mat = xOut(:,(2*n+1):(3*n)); 

%plotting:
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


%selecting segment for analysis:
t_range = 300; 
t_start = t(end) - t_range; 
[~,t_start_idx] = min(abs(t-t_start)); 
t_seg = t(t_start_idx:length(t));
a_mat_seg = a_mat(t_start_idx:length(t),:); 

[peak_data, fixed_pt, base_osc, ind_indicators, grp_indicators] = InterestingBehavior(n,t_seg,a_mat_seg,params);

sorted_peak_data = sortrows(peak_data); 

cluster_mat = [(1:length(peak_data))', sorted_peak_data(:,1)];

figure
hold on
plot(1:length(peak_data),sorted_peak_data(:,1),'o')
hold off

idx = dbscan(cluster_mat,10, 2)
gscatter(cluster_mat(:,1), cluster_mat(:,2),idx)

%plotting segment
figure
hold on
for i=1:n
plot(t_seg,a_mat_seg(:,i), 'linewidth', 1.5)
Legend{i} = ['Region ', num2str(i)];
end
plot(peak_data(:,1),peak_data(:,2),'r*')
xlabel('Time','fontsize',12,'fontweight', 'bold')
ylabel('Artist Densities','fontsize',12,'fontweight', 'bold')
title('Artist Densities vs Time','fontsize',12,'fontweight', 'bold')
legend(Legend,'fontsize', 12,'fontweight', 'bold')
grid on; grid minor
hold off
%}

toc
