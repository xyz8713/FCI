function y = ratPrecD3_c(prec, rhs)
%% function y = ratPrec(PRE, rhs)
%% PRE = struct for preconditioner
%% USE P as the Preconditioner
%%--------------------
    method = 2;  %% 2 Both projector and ILU 
    method = 3;  %% 3 projector only 
    method = 1;  %% 1 ILU only
   
    method = 2;  %% 2 Both projector and ILU 

    PRE = prec.PRE;
    %%    B = prec.B;
    rlcase = isreal(prec.B);
    coefs2 = prec.coefs2;
    nC = prec.nC;
    if (~rlcase)
        nC = nC+nC;
    end
    y = zeros(size(rhs));
%%---- form P*rhs
   if (method >= 2)      
       parfor k =1:nC
           Lk = PRE(k).Lk;
           Uk = PRE(k).Uk;
           cp2 = coefs2(k); 
           y = y + cp2*Uk\( Lk \ rhs);
       end
   else
       y = rhs;
   end
   
   if (method <=2) 
        k = nC+1;
        if (rlcase)
            k = 2*nC+1;
        end
        Lk = PRE(k).Lk;
        Uk = PRE(k).Uk;
        y = Uk \ (Lk \ y);
        %% in case I use a complex L/U
    end
    if (rlcase) 
        y = real(y);
    end
end

