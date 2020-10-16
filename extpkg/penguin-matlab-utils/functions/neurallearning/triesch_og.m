function [a, b] = triesch_og(input, output, a, b, eta_ip, mu_output) 
da = eta_ip*(1./a + input - (2+1./mu_output).*input.*output + input.*output.*output./mu_output);
db = eta_ip*(1 - (2+1./mu_output).*output + output.*output./mu_output);
a = a + da; 
b = b + db; 
end
