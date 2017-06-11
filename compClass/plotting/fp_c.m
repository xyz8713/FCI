function y = fp_c(coefs, shift, x)
n = length(x);
x = x(:);
y = zeros(n,1);
%%-------------------- loop over poles
for ii =1:length(shift)
    z = 1./  (x - shift(ii)) ;
    y = y+coefs(ii)*z;
end
%%-------------------- this assumes both parts of 
%%                     real axis. otherwise double it
y = real(y);
end
