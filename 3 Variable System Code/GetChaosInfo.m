function GetChaosInfo(n, t, A0, N0, C0, a_mat, n_mat, c_mat, params, ModelType)

%generating close initial conditions: 
delta_A = 1e-12; 
delta_N = 1e-12; 
delta_C = 1e-12; 

A0_new = A0; 
N0_new = N0;
C0_new = C0; 

if mod(n,2) %testing for odd n
    for i = 1:n-1
        A0_new(i) = A0(i) + (-1)^(i-1)*delta_A;
        N0_new(i) = N0(i) + (-1)^(i-1)*delta_N;
        C0_new(i) = C0(i) + (-1)^(i-1)*delta_C;
    end
else
    for i = 1:n
        A0_new(i) = A0(i) + (-1)^i*delta_A;
        N0_new(i) = N0(i) + (-1)^i*delta_N;
        C0_new(i) = C0(i) + (-1)^i*delta_C;
    end
    
end

%simulation (by type): 

if ModelType == 'Non-spatial'
    [a_mat_new, n_mat_new, c_mat_new] = Euler3VarSim(n,A0_new,N0_new,C0_new,params); 
elseif ModelType == 'SpatialRing'
    [a_mat_new, n_mat_new, c_mat_new] = Euler3VarSimSpatialRing(n,A0_new,N0_new,C0_new,params); 
end

%testing for chaos: 
[norm_dif] = GetDifNorm(a_mat, n_mat, c_mat, a_mat_new, n_mat_new, c_mat_new); 

%plotting:
figure
hold on
for i=1:n
plot(t,a_mat(:,i),'r', 'linewidth', 1.3)
Legend{i} = ['Region ', num2str(i)];
end
xlabel('Time','fontsize',12,'fontweight', 'bold')
ylabel('Artist Densities','fontsize',12,'fontweight', 'bold')
title('Artist Densities vs Time','fontsize',12,'fontweight', 'bold')
legend(Legend,'fontsize', 12,'fontweight', 'bold')
grid on; grid minor

for i=1:n
plot(t,a_mat_new(:,i),'b', 'linewidth', 1.3)
Legend{i} = ['Region ', num2str(i)];
end
xlabel('Time','fontsize',12,'fontweight', 'bold')
ylabel('Artist Densities','fontsize',12,'fontweight', 'bold')
title('Artist Densities vs Time','fontsize',12,'fontweight', 'bold')
legend(Legend,'fontsize', 12,'fontweight', 'bold')
grid on; grid minor
hold off



%plotting norm of differences vs time: 
figure
hold on
plot(t,log(norm_dif), 'linewidth', 1.5)
xlabel('Time', 'fontsize', 12, 'fontweight', 'bold')
ylabel('ln|difference|', 'fontsize', 12, 'fontweight', 'bold')
title('Norm of Difference Between State Vectors vs Time', 'fontsize', 12, 'fontweight', 'bold')
grid on
grid minor
hold off

end