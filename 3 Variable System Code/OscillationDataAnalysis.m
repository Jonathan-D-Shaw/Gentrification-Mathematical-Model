%% Gentrification Model: Oscillatory Solutions Analysis

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

%% Main

load('OscillationData.mat')

%storing tau_n, tau_c data in matrices for quicker plotting: 
fixed_point_mat = zeros(length(fixed_point),1); 
base_oscillation_mat = zeros(length(base_oscillations),1); 
grp_oscillation_mat = zeros(length(grp_oscillations),1); 
chaos_mat = zeros(length(chaos),1); 


    for i = 1:length(fixed_point)
         fixed_point_mat(i,1) = fixed_point(i).parameters(2);
         fixed_point_mat(i,2) = fixed_point(i).parameters(3); 
    end
    
    for i = 1:length(base_oscillations)
         base_oscillations_mat(i,1) = base_oscillations(i).parameters(2);
         base_oscillations_mat(i,2) = base_oscillations(i).parameters(3); 
    end
    
    for i = 1:length(grp_oscillations)
        grp_oscillations_mat(i,1) = grp_oscillations(i).parameters(2);
         grp_oscillations_mat(i,2) = grp_oscillations(i).parameters(3);
    end
    
    for i = 1:length(chaos)
        chaos_mat(i,1) = chaos(i).parameters(2); 
        chaos_mat(i,2) = chaos(i).parameters(3); 
    end

    
%generating linear stability curve: 
n = 5; 
a=1/n; 
z = 0.01; 
tau_n = 1:0.2:50; 
epsilon =0.1;

sig_p = sigma_prime(a,z,epsilon); 

tau_c = 1./(sig_p - 1./tau_n);
num2remove = length(find(tau_c<1)); %number of values of tau_c<1

%removing indices corresponding to tau_c<1:
tau_c = tau_c((num2remove+1):end);
tau_n_temp = tau_n((num2remove+1):end); 
    
%plotting:    
figure('HandleVisibility','off');
hold on
plot(tau_n_temp,tau_c,'linewidth',1.5)
plot(fixed_point_mat(:,1), fixed_point_mat(:,2), 'ro')
plot(base_oscillations_mat(:,1), base_oscillations_mat(:,2), 'bo')
plot(grp_oscillations_mat(:,1), grp_oscillations_mat(:,2), 'go')
plot(chaos_mat(:,1), chaos_mat(:,2), 'ko')
xlabel('\tau_N', 'fontsize', 12, 'fontweight', 'bold')
ylabel('\tau_C', 'fontsize', 12, 'fontweight', 'bold')
title('Behavior of Oscillatory Solutions as a Function of \tau_N and \tau_C', 'fontsize', 12, 'fontweight', 'bold')
legend('Linear Stability Condition', 'Fixed Point', 'Basic Oscillations', 'Group Oscillations', 'Chaos', 'fontsize', 12, 'fontweight', 'bold')
hold off

    