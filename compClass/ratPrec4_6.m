 function [y] = ratPrec4_6(prec, rhs)
%% function y = ratPrec_6(PRE, rhs)
 %% one circle only  --  subspace preconditioner for
 %% gmres when solving second part. 
 %% method = 0 corresponds to inner-outer gmres
 %% with same shift for system with B. 
%%--------------------
opts   = prec.ITopts;
method = opts.outer;

n = length(rhs);
y = zeros(n,1);
B = prec.B;
%%-------------------- 
    opts.Type = 1;
    opts.lev  = 1;
    %%    
    rhsi = rhs;
    if (method == 0)
        rhsi = rhs;
        for iter = 1:1
            y3 = fgmrez_dr2(B,prec,'nopre',rhsi,rhsi,opts);
            y = y+y3;
            rhsi = rhs-B*y;
            opts.shift = opts.shift*0.25;
        end
        return;
    end
    %%-------------------- methods 1, 2, and 3 are the same. 
    [y1, y2] = ffom_mp(prec,rhsi,opts);

    %%-------------------- Note: when prec.imSubs==0 subSol== nopre.
    y3 = fgmrez_dr2(B,prec,'subsSol',y2,y2,opts);
    y = y1+y3;

 


