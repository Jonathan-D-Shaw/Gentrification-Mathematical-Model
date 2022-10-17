%function approximating h(x): 
function fx = f_x(del_mat)

fx = del_mat.*(del_mat>0); 

%{
m = 709; 
fx = log(1+exp(m*del_mat))/m; 
%}
end