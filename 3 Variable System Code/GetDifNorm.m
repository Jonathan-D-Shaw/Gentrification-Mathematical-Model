%function to test for chaos: 
function [norm_dif] = GetDifNorm(a_mat1, n_mat1, c_mat1, a_mat2, n_mat2, c_mat2)

%conacenating everything: 

state_mat1 = [a_mat1, n_mat1, c_mat1]; 
state_mat2 = [a_mat2, n_mat2, c_mat2]; 
norm_dif = zeros(length(state_mat1),1); 
    for i = 1:length(state_mat1)
        dif_vec = state_mat1(i,:) - state_mat2(i,:); 
        norm_dif(i) = norm(dif_vec); 
    end
end