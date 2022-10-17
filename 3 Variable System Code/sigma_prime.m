%d/da(sigma(a,z,epsilon)):
function sig_p = sigma_prime(a,z,epsilon)
sig_p = (sech((a-z)./epsilon).^2)./(2*epsilon); 
end