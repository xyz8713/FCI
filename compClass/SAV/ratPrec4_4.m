 function [y] = ratPrec4_4(prec, rhs)
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
    opts.Type = 1;
    opts.lev  = lev;
    [y1, y2] = ffom_mp(prec,rhsi,opts);
%%-------------------- Solve Ax = P y2
    if (method <=2)
        if(method == 2) 
            opts.shift = opts.shift*0.25;
        end
        [y3,res2,its2] = fgmrez_dr2(B,prec,'ratPrecD3_d',y2,y2,opts);
        %%        [y3,res2,its2] = fgmrez_dr2(B,prec,'nopre',y2,y2,opts);
    elseif(method == 3) 
        opts   = prec.ITopts;
        it = opts.maxits*lev;
        tol = opts.tolIts;
        sh = opts.shift/lev;
        M = B - sh*speye(n);
        [y3,res2] = CGNR (M,zeros(n,1), y2, it, tol) ;
    end
    y =  y + (y1 + y3);
    %%# debug 
    %% - ex = B\ones(n,1);
    %% - er1 = real(ex- y1);
    %% - er2 = real(ex -B\y2);
    %% - er = real(ex -y2);
    %% - yy = real(rhs - B*y);
    %% - plot(yy(1:50))
    %% - keyboard
    %% - error(' stop here') 
    %% - %y =  y + y3;
    rhsi = rhs-B*y;
end
 


