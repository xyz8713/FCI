classdef precClass4_0 < handle 
%% original one -- as in paper 
%% 
 properties
     nC
     r 
     PRE
     coefs1
     coefs2
     ITopts
     U 
     B
 end

methods
  function prec = precClass4_0(B,nC,r,ITopts)
  %%-------------------- CONSTRUCTOR --------------------
  %%-------- nC = number of poles on the upper half plane we need 2nC for
  %%  non-Hermitian case               
  %% construct preconditoner class. 
  %%  r = 30;
     bDo_sILU = 1;
     if (ITopts.inner == 0) 
         bDo_sILU = 0;
     end
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

     sigma  = z/r;
     omega  = om/r;
     c      = -max(real(z))*r-1;
     coefs1 = -omega./(sigma + c*sigma.*sigma);
     shift  = 1./sigma +c;
     coefs1 = [coefs1 conj(coefs1)];
     shift  = [shift conj(shift)];
     %%--------------------=
     %%     shift2 = [shift2, conj(shift2)];
     omega2 = conj(om*r);
     coefs2 = omega2;
     coefs2 = [coefs2 conj(coefs2)];
     
     %%     shift2 = c -sigma2;
     %%     shift2 = shift;
     %%     shift2 = ones(size(sigma2));
     if (0)
         opts.type = 'crout';
         opts.milu = 'off';
         opts.droptol = 0.01;
     else
         opts.type = 'ilutp';
         opts.droptol = 0.01;
     end
     tic;
     %Lk = cell(2*nC,1);
     %Uk = cell(2*nC,1);
    for k = 1:2*nC  
        if (isreal(B) & k>nC)
            k1 = k-nC;
            Lk = conj(PRE(k1).Lk);  
            Uk = conj(PRE(k1).Uk);
        else
            fprintf(1,'factorize the %dth shifted matrix\n',k);          
            [Lk, Uk] = ilu(B-speye(n)*shift(k),opts);
        end
        PRE(k) = struct('Lk',Lk,'Uk',Uk) ;
    end 
    %%----------------    find a good complex shift for ILU. 
    opts.droptol = 0.01;   
    if (bDo_sILU) 
         if (1)
             t = find_shift(B);
             opts.type = 'nofill';
             [Lk, Uk] = ilu(B+t*speye(n),opts);
        else
            d = find_shiftd(B);
            [Lk, Uk] = ilu(B+spdiags(d,0,n,n),opts);
        end
        (nnz(Lk)+nnz(Uk)-n)/nnz(B)
         PRE(2*nC+1) = struct('Lk',Lk,'Uk',Uk) ;
    end
    t2 = toc;
    fprintf(1,'factorize time is %f\n',t2);
    %%-------------------- initial U for subspace preconditioning
    
    prec.U = [];
    prec.nC = nC;
    prec.PRE = PRE;
    prec.coefs1 = coefs1;
    prec.coefs2 = coefs2;
    prec.B = B;
    prec.ITopts = ITopts;
end
end
end
