 function [y] = ratPrec4_5(prec, rhs)
%% function y = ratPrec(PRE, rhs)
%% main precconditioning function 
%% PRE = struct for preconditioner
%%--------------------
 %% method = 1 ; standard combination rational + solve with ratPrecD3_c   
 %% method = 2 ; solve A P x = Pb by a two-shift method... [test]
%%--------------------
opts   = prec.ITopts;
method = opts.outer;
Nlev  = prec.Nlev;

n = length(rhs);
y = zeros(n,1);
B = prec.B;
rhsi = rhs;
%%-------------------- iteration circles 
for lev = 1:Nlev
    opts.Type = 2;
    opts.lev  = lev;
    [y1, unused] = ffom_mp(prec,rhsi,opts);
%%-------------------- USE P as the preconditioner to solve Ax = y2
    y =  y + y1;
    %%    rhsi = rhs-B*y;
end



