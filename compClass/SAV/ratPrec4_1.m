 function [y] = ratPrec4_1(prec, rhs)
%% function y = ratPrec(PRE, rhs)
%% main precconditioning function 
%% PRE = struct for preconditioner
%%--------------------
InnerIt = 8;
PRE = prec.PRE;
nC = prec.nC;
coefs1 = prec.coefs1;
coefs2 = prec.coefs2;
%%shift = prec.shift;
n = length(rhs);
y = zeros(n,1);
B = prec.B;
rlcase = isreal(B);
%%rlcase = 0;
if (~rlcase) 
    nC= nC+nC;
end

for it =1: InnerIt 
 r = rhs - B*y;
 for k =1:nC
     Lk = PRE(k).Lk;
     Uk = PRE(k).Uk;
     cp = coefs1(k);
     yr = Uk\( Lk \ r);
     y = y + cp*yr;
 end
    if (rlcase)
        y = y+conj(y);
    end
 r  = rhs - B*y;
 for k =1:nC
     Lk = PRE(k).Lk;
     Uk = PRE(k).Uk;
     cp = coefs2(k); 
     yr = Uk\( Lk \ r);
     y = y + cp*yr;
 end
    if (rlcase)
        y = y+conj(y);
    end
end
end


