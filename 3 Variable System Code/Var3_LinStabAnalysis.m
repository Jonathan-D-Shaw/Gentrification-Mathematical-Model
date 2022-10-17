%% Gentrification Model: 3-Variable System Linear Stability Analysis

%{
Authors: 
Juan G. Restrepo
Nancy Rodriguez
Jonathan D. Shaw

Date created: 06/15/2022
Date modified: 07/16/2022
%}

%housekeeping: 
clear; close all; clc; 

%% Phase Plot
%
n =4; 
a=1/n; 
z = 0.01; 
tau_n = 0.1:0.1:20; 
epsilon =0.1;


sig_p = sigma_prime(a,z,epsilon); %sigma'(a*). see overleaf doc for details

%
figure
hold on

     tau_c = 1./(sig_p - 1./tau_n);
    
    %plotting: 
%     plot(tau_n,tau_c,'linewidth',1.5); 
%     xline(1/sig_p, 'r-.', 'linewidth', 1.2); yline(1/sig_p, 'r-.', 'linewidth', 1.2)
%     Legend{i} = strcat('n=', num2str(n), ', z=', num2str(z), ', \epsilon=', num2str(epsilon(i)));


xlabel('\tau_n','fontweight','bold','fontsize',12)
ylabel('\tau_c','fontweight','bold','fontsize',12)
title('Stabilty Phase Plot','fontweight','bold','fontsize',12)
% legend(Legend, 'fontweight','bold','fontsize',12)
grid on; grid minor
hold off

%surface to directly show inequality: 
tau_n = 0.1:0.1:20;
[x,y] = meshgrid(tau_n); 
unstable_cond = 1./x + 1./y > sig_p;
surf(x,y,double(unstable_cond))
view(0,90)

