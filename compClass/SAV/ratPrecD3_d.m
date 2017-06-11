 function y = ratPrecD3_d(prec, rhs)%% function y = ratPrec(PRE, rhs)
%% PRE = struct for preconditioner
%% USE P as the Preconditioner
%%--------------------
%% method = 0;  %% no preconditiong [Default!]
%% method = 1;  %% 1 ILU preconditioning only
%% method = 2;  %% 2 Both projector and ILU 
%% method = 3;  %% 3 projector only     
    opts = prec.ITopts;
    method = opts.inner;

    if (method == 0)
        y = rhs;
        return
    end
    %%
    rlcase = isreal(prec.B);
%%---- form P*rhs
   if (method >= 2)      
       opts.Type    =  1;
       opts.lev    =  1;
       [y,ydummy] = ffom_mp(prec,rhs,opts);
   else
       y = rhs;
   end
   %%   
   if (method <=2) 
        Lk = prec.Lk;
        Uk = prec.Uk;
        y = Uk \ (Lk \ y);
        %% in case I use a complex L/U
    end
    if (rlcase) 
        y = real(y);
    end
end

