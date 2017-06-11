 function [y] = ratPrec4_0(prec, rhs)
%% function y = ratPrec(PRE, rhs)
%% main precconditioning function 
%% PRE = struct for preconditioner
%%--------------------
 %% method = 1 ; standard combination rational + solve with ratPrecD3_c   
 %% method = 2 ; PAinv only --  with ratPrecD3_c  -- not (I-P)Ainv part
 %% method = 3 ; (I-P)Ainv contour integral only -- ignore inside of circle.
 
opts = prec.ITopts;
method = opts.outer;

PRE = prec.PRE;
nC = 2*prec.nC;

y  = zeros(size(rhs));
y1 = zeros(size(rhs));
y2 = zeros(size(rhs));

n = length(rhs);
B = prec.B;
n = size(B,1);

if (method ==2) 
    y2 = rhs;
else
    coefs1 = prec.coefs1;
    coefs2 = prec.coefs2;
    %%%s2 = prec.shift;
    %%---- form P*rhs
    for k =1:nC
        Lk  = PRE(k).Lk;
        Uk  = PRE(k).Uk;
        cp  = coefs1(k);
        cp2 = coefs2(k); 
        yr  = Uk\( Lk \ rhs);
        y1  = y1 + cp*yr;
     
        if (method ~= 3) 
            y2  = y2 + cp2*yr;
        end
    end
end
if (method == 3) 
    y = y1;
    return
end
%%-------------------- USE P as the preconditioner to solve Ax = y2
    %%                          %% Nevec = 0 means no deflation
    [y3,res2,its2] = fgmrez_dr2(B,prec,'ratPrecD3_c',y2,y2,opts);
%------y3 = B\y2;
    y =  y1+y3 ;
end


