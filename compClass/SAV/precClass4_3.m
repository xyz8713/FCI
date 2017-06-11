classdef precClass4_3 < handle 
%% original one -- as in paper 
%% 
 properties
     nC
     r 
     coefs1
     coefs2
     ITopts
     shift
     Lk       %% precon for 2nd part 
     Uk       %% optional ..  
     B  
 end

methods
  function prec = precClass4_3(B,nC,r,ITopts)
  %%-------------------- CONSTRUCTOR --------------------
  %%-------- nC = number of poles on the upper half plane we need 2nC for
  %%  non-Hermitian case               
  %% construct preconditoner class. 
  %%  r = 30;
     bDo_sILU = 0;
     if ((ITopts.outer == 0) | (ITopts.inner)==1 | (ITopts.inner)==2 )
         bDo_sILU = 1;
     end
     %%

     n = size(B,1);
     %%-------- intial center
     c = 0.0;
     %% Part 4: Compute the weights and nodes for rational approximation
     %% to h                            
     %%-------- Get nodes (z) and weights (om) on a unit circle 
     [z, om] = contQuad(nC,2);
     %%-------- Get rational approximation to 1/z outside circle 
     %% 1/z \approx \sum omega(k)/(sigma(k)+c\sigma(k)^2)
     %% *(1/(1/sigma(k)+c-z))                          
     z = z(:); om = om(:);
     sigma  = z/r;
     omega  = om/r;
     c      = -max(real(z))*r-1;
     coefs1 = -omega./(sigma + c*sigma.*sigma);
     shift  = 1./sigma +c;
     coefs1 = [coefs1; conj(coefs1)];
     shift  = [shift ; conj(shift)];
     %%--------------------=
     %%     shift2 = [shift2, conj(shift2)];
     coefs2 = conj(om*r);
     coefs2 = [coefs2;  conj(coefs2)];
     %%-------------------- ILU preconditoner case.
     if (0)
         opts.type = 'crout';
         opts.milu = 'off';
         opts.droptol = 0.01;
     else
         opts.type = 'ilutp';
         opts.droptol = 0.01;
     end
     tic;
    %%----------------    find a good complex shift for ILU. 
    prec.nC = nC;
    prec.coefs1 = coefs1;
    prec.coefs2 = coefs2;
    prec.shift= shift;
    prec.B = B;
    prec.ITopts = ITopts;
    %%-------------------- get Lk Uk factors in certain cases.
    opts.droptol = 0.01;
    if (bDo_sILU) 
        if (1)
            t = find_shift(B);
            [Lk, Uk] = ilu(B+t*speye(n),opts);
        else
            d = find_shiftd(B);
            [Lk, Uk] = ilu(B+spdiags(d,0,n,n),opts);
        end
        (nnz(Lk)+nnz(Uk)-n)/nnz(B)
        t2 = toc;
        fprintf(1,'factorize time is %f\n',t2);
        %%-------------------- initial U for subspace preconditioning
        prec.Lk = Lk;
        prec.Uk = Uk;
    end
end
end
end
