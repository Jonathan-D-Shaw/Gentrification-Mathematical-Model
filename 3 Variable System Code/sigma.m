%sigma function. see overleaf doc for details as needed. 
function sig = sigma(a,z,epsilon)
sig = 0.5 * (1 + tanh((a-z)/epsilon)); 
end
