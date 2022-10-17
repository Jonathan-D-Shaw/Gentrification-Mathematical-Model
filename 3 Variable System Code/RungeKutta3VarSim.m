function dxdt = RungeKutta3VarSim(t,x,n,params)

%extracting variables:
var_mat = zeros(3,n); 
    for i = 1:3
        var_mat(i,:) = x(1:n); 
        x(1:n) = []; 
    end
    
    
A = var_mat(1,:)'; 
N = var_mat(2,:)'; 
C = var_mat(3,:)'; 

difC = C*params.un'-params.un*C'; %matrix of cost differences

%passing to function acting on cost differences: 
f_C = f_x(difC);    

%sums in matrix form: 
Adot = 1/params.tau_a * (f_C'*A - A.*(f_C*params.un)) + 0*rand(n,1);
Ndot = 1/params.tau_n * (sigma(A,params.z,params.epsilon) - N) + 0*randn(n,1);
Cdot = 1/params.tau_c * (N-C); 

dxdt = [Adot; Ndot; Cdot]; 
end