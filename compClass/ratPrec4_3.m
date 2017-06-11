  function [y] = ratPrec4_3(prec, rhs)
%% function y = ratPrec(PRE, rhs)
%% main precconditioning function 
%% PRE = struct for preconditioner
%%--------------------
%% outer-method = 1; standard combination rational + solve with ratPrecD3_c  
%% outer-method = 2; only with ratPrecD3_c   -- not countour integral
%% outer-method = 0; ILU preconditioning only -- FOR
%% TESTING/COMPARISONS.  

%%--------------------
opts   = prec.ITopts;
method = opts.outer;

%%====================>> this is to test ILU only *BY ITSELF* 
if (method == 0)
    Lk = prec.Lk;
    Uk = prec.Uk;
    y =  Uk\ (Lk\rhs);
    return
end

n = length(rhs);
B = prec.B;
n = size(B,1);
y1 = zeros(n,1);
if (method == 2) 
    y2 = rhs;
else
    opts.Type = 1;
    opts.lev = 1;
    [y1, y2] = ffom_mp(prec,rhs,opts);
end
%%-------------------- USE P as the preconditioner to solve Ax = y2
%% for compatibility only:
%%    opts.shift = 0.0;
    [y3,res2,its2] = fgmrez_dr2(B,prec,'ratPrecD3_d',y2,y2,opts);
%------y3 = B\y2;
%------y3 = B\y2;
    y =  y1+y3;
 end


