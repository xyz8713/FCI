function y = fh_c(coefs, shift, x)
n = length(x);
x = x(:);
y = zeros(n,1);
%%-------------------- loop over poles
for ii =1:length(shift)
    z = 1./  (shift(ii) -  x) ;
    y = y+coefs(ii)*z;
end
y = real(y)/2;
end
