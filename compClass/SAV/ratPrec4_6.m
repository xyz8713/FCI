 function [y] = ratPrec4_6(prec, rhs)
%% function y = ratPrec_6(PRE, rhs)
%% main precconditioning function 
%% PREc = struct for preconditioner
%%--------------------
%%--------------------
opts   = prec.ITopts;
method = opts.outer;

n = length(rhs);
y = zeros(n,1);
B = prec.B;
%%-------------------- 
    opts.Type = 1;
    opts.lev  = 1;
    
    rhsi = rhs;
    
    [y1, y2] = ffom_mp(prec,rhsi,opts);
        
    y3 = fgmrez_dr2(B,prec,'subsSol',y2,y2,opts);
    y = y1+y3;

    end
 


