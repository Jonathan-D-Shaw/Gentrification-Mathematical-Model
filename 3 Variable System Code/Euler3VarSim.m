%3-variable system simulation (Euler's method)
function [a_mat, n_mat, c_mat] = Euler3VarSim(n,A,N,C,params)

%initializing zeros matrices for computational efficiency:
a_mat = zeros(length(params.t),n); c_mat = a_mat; n_mat = a_mat;

%populating first rows with IC:
a_mat(1,:) = A; 
n_mat(1,:) = N; 
c_mat(1,:) = C; 

    for i = 1:(length(params.t)-1)

        difC = C*params.un'-params.un*C'; %matrix of cost differences

        %passing to function acting on cost differences: 
        f_C = f_x(difC);    

        %derivative calcs: 
        Adot = 1/params.tau_a * (f_C'*A - A.*(f_C*params.un)) + 0*rand(n,1);
        Ndot = 1/params.tau_n * (sigma(A,params.z,params.epsilon) - N) + 0*randn(n,1);
        Cdot = 1/params.tau_c * (N-C); 

        A = A + Adot*params.dt;
        N = N + Ndot*params.dt; 
        C = C + Cdot*params.dt;

        a_mat(i+1,:) = A'; 
        n_mat(i+1,:) = N'; 
        c_mat(i+1,:) = C';  
    end
end