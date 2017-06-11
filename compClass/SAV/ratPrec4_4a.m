 function [y] = ratPrec4_4a(prec, rhs)
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
B = prec.B;
n = size(B,1);
y = zeros(n,1);
rhsi = rhs;
%%-------------------- iteration circles 
for lev = 1:Nlev
    opts.Type = 3;
    opts.lev  = lev;
    [y1, y2] = ffom_mp(prec,rhsi,opts);
    y = y+y2;
    rhsi = rhs-B*y;
end
 


