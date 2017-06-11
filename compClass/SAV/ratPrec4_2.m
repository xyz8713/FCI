 function [y] = ratPrec4_2(prec, rhs)
%% function y = ratPrec(PRE, rhs)
%% main precconditioning function 
%% PRE = struct for preconditioner
%%--------------------
InnerIt = 6;
PRE1 = prec.PRE1;
PRE2 = prec.PRE2;
nC1 = prec.nC1;
nC2 = prec.nC2;
rlcase = isreal(prec.B);
%%rlcase = 0;
if (~rlcase) 
    nC1 = nC1+nC1;
    nC2 = nC2+nC2;
end
coefs1 = prec.coefs1;
coefs2 = prec.coefs2;
n = length(rhs);
y = zeros(n,1);
B = prec.B;
%%----

for it =1: InnerIt 
    r = rhs - B*y;
    d = zeros(n,1);
for k =1:nC1
     Lk = PRE1(k).Lk;
     Uk = PRE1(k).Uk;
     cp = coefs1(k);
     yr = Uk\( Lk \ r);
     d = d + cp*yr;
 end
 if (rlcase)
     d = 2*real(d);
 end
 y = y+d;
 r = rhs - B*y;
 d = zeros(n,1);
 for k =1:nC2
     Lk = PRE2(k).Lk;
     Uk = PRE2(k).Uk;
     cp = coefs2(k);
     yr = Uk\( Lk \ r);
     d = d + cp*yr;
 end
 if (rlcase)
     d = 2*real(d);
 end
 y = y +d;
end
end


