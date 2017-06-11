classdef precClass4_2 < handle 
%% 2 circles.
 properties
     nC1
     nC2
     r1
     r2
     PRE1
     PRE2
     coefs1
     coefs2
     B 
end

methods

    function prec = precClass4_2(B, nC1, r1, nC2, r2) 
%%  function prec = precClass4(B, nC1, r1, nC2, r2) 
%%-------------------- CONSTRUCTOR --------------------
%% construct preconditoner class. 
    n = size(B);
    if (nargin < 4)
        nC2 = nC1;
        r2 = r1;
    end
    bDo_Plot = 0;
    %%non-Hermitian case               
    %%-------- Get rational approximation to 1/z outside circle                 
    %% 1/z \approx \sum omega(k)/(sigma(k)+c\sigma(k)^2)
    %% *(1/(1/sigma(k)+c-z))    
    %%-------------------- circle
    [z1, om1] = contQuad(nC1,2);
    sigma1 = z1/r1;
    omega1 = om1/r1;
    c1 = -max(real(z1))*r1-1;
    %%-------------------- first projector
    coefs1 = -omega1 ./(sigma1 + c1*sigma1.*sigma1);
    shift1 = 1./sigma1 + c1;
    %%x2
    coefs1 = [coefs1 conj(coefs1)];
    shift1 = [shift1, conj(shift1)];
    %%-------------------- first projector
    [z2, om2] = contQuad(nC2,2);
    %% PRE2:
    %%% r2 = (r1-c1)/2;
    if (1) 
        sigma2 = z2/r2;
        omega2 = om2/r2;
        c2 = -min(real(z2))*r2+1;
        %%-------------------- second projector
    coefs2=omega2 ./(sigma2+c2*sigma2.*sigma2);
    shift2 = 1./sigma2 + c2;
    %%x2 
    else
        c2 = -r2*1.05;
        shift2 = c2+z2*r2;   %% refactor.
        coefs2 = om2*r2 ; %%./ shift2;
    end
    coefs2 = [coefs2 conj(coefs2)];
    shift2 = [shift2, conj(shift2)];
    %%x2 
      %%
    if (bDo_Plot)
        figure(1)
        compp_c(coefs2,shift2,r2,c2);
        pause
        close
    end
    %%-------------------- tye of precon.
    if (1)
        opts.type = 'ilutp';
        opts.milu = 'off';
        opts.droptol = 0.01;
    else
        opts.type = 'ilutp';
        opts.droptol = 0.01;
    end
    tic;
%%-------------------- factors for 1st projector 
    %%parfor k = 1:2*(nC+nC1)
    for k = 1:2*nC1  
        if (isreal(B) & k>nC1)
            k1 = k-nC1;
            Lk = conj(PRE1(k1).Lk);  
            Uk = conj(PRE1(k1).Uk);
        else
            fprintf(1,'factorize the %dth shifted matrix\n',k);
            [Lk, Uk] = ilu(B-speye(n)*shift1(k),opts);
        end
        PRE1(k) = struct('Lk',Lk,'Uk',Uk) ;
    end 
    %%-------------------- factorizations for 2nd projector 
    for k = 1:2*nC2  
        if (isreal(B) & k>nC2)
            k1 = k-nC2;
            Lk = conj(PRE2(k1).Lk);  
            Uk = conj(PRE2(k1).Uk);
        else
            fprintf(1,'factorize the %dth shifted matrix\n',k);
            [Lk, Uk] = ilu(B-speye(n)*shift2(k),opts);
        end
        PRE2(k) = struct('Lk',Lk,'Uk',Uk) ;
    end 

     t2 = toc;
     fprintf(1,'factorize time is %f\n',t2);
     %%-------------------- data for subspace preconditioning
     prec.nC1 = nC1;
     prec.nC2 = nC2;
     prec.PRE1 = PRE1;
     prec.PRE2 = PRE2;
     prec.coefs1 = coefs1;
     prec.coefs2 = coefs2;
     prec.B = B;
  end
 end
end
