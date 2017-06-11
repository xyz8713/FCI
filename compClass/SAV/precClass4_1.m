classdef precClass4_1 < handle 

 properties
     nC
     r
     PRE
     coefs1
     coefs2
     %%     shift
     B 
end

methods

    function prec = precClass4_1(B, nC, r) 
%%  function prec = precClass4_1(B, nC, r) 
%%-------------------- CONSTRUCTOR --------------------
%% construct preconditoner class. 
    n = size(B);
    %%non-Hermitian case               
    %%-------- Get rational approximation to 1/z outside circle                 
    %% 1/z \approx \sum omega(k)/(sigma(k)+c\sigma(k)^2)
    %% *(1/(1/sigma(k)+c-z))    
    %%-------------------- circle
    [z1, om1] = contQuad(nC,2);
    sigma1 = z1/r;
    omega1 = om1/r;
    c1 = -max(real(z1))*r-1;
    %%-------------------- first projector
    coefs1 = -omega1 ./(sigma1 + c1*sigma1.*sigma1);
    shift  = 1./sigma1 + c1;
    %%##--= Very strange that coefs*4 works better -- 
    coefs1 = [coefs1 conj(coefs1)];
    shift  = [shift, conj(shift)];
    %%-------------------- second projector
    %%---> : recycle factorization
    coefs2 = om1 ./ z1; %%% shift(1:nC);
    %%coefs2 = conj(r*om1);
    coefs2 = [coefs2 conj(coefs2)];
    %%-------------------- type of precon.
    if (1)
        opts.type = 'ilutp';
        opts.milu = 'off';
        opts.droptol = 0.01;
    else
        opts.type = 'ilutp';
        opts.droptol = 0.01;
    end
    tic;
    %%-------------------- factorizations for 1st projector 
    %%parfor k = 1:2*(nC+nC)
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
    %%-------------------- factorizations for 2nd projector 
    t2 = toc;
    fprintf(1,'factorize time is %f\n',t2);
    %%-------------------- data for subspace preconditioning
    prec.nC = nC;
    prec.PRE = PRE;
    prec.coefs1 = coefs1;
    prec.coefs2 = coefs2;
    prec.B = B;
    %%prec.shift = shift;
 end
end
end
