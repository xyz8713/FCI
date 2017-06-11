classdef precClass4_4 < handle 
%% original one -- matvec-only preconditioner 
%% 
 properties
     Nlev
     nC
     r         %% now a vector of length Nlev 
     coefs1    %% 2D array coefs1(:,k) = coefs for circle k
     coefs2    %% similar coefs for second part
     ITopts
     shift    %% 2D array shift(:,k) = shifts for circle k
     B  
     Lk
     Uk
 end

methods
  function prec = precClass4_4(B,nC,r,ITopts)
  %%-------------------- CONSTRUCTOR --------------------
  %% nC     = number of countour points [upper] same for each circle
  %% r      = r(1:nlev) = radii of the circles
  %% ITopts = parameters for solver 
  %%-------- nC = number of poles on the upper half plane we need 2nC for
  %%  non-Hermitian case               
  %% construct preconditoner class. 
     n = size(B,1);
     %%-------- intial center
     c = 0.0;
     %% Part 4: Compute the weights and nodes for rational approximation
     %% to h                            
     %%-------- Get nodes (z) and weights (om) on a unit circle 
     [z, om] = contQuad(nC,2);
     z = z(:); om = om(:);
     %%-------- Get rational approximation to 1/z outside circle 
     %% 1/z \approx \sum omega(k)/(sigma(k)+c\sigma(k)^2)
     %% *(1/(1/sigma(k)+c-z))                          
     %% Nlev   = number of levels
     Nlev   = length(r);
     sgn = -1;
     for lev = 1:Nlev
         %%-------------------- for each circle do:
         rj     = r(lev);    %% radius of jth circle 
         sigma  = z/rj;
         omega  = om/rj;
         c      = sgn*(max(real(z))*rj)-sgn;
         sgn    = -sgn;
         co1    = -omega./(sigma + c*sigma.*sigma);
         sh     = 1./sigma +c;
         %%     shift2 = [shift2, conj(shift2)];
         co2    = conj(om*rj);
         coefs1(:,lev) = [co1; conj(co1)];
         coefs2(:,lev) = [co2; conj(co2)];
         shift(:,lev)  = [sh; conj(sh)];
     end
    %%----------------    find a good complex shift for solving
     d = abs(B)*ones(n,1)-abs(diag(B)); 
     d = diag(B)-d;
     t = min(d);
     %% will make matrix B diagonally dominant. 
     ITopts.shift = 0.5*t;
     %%  ITopts.shift = 0.0;
    if (ITopts.inner==1 | ITopts.inner==2)
        tic;
        %%-------------------- get Lk Uk factors in certain cases.
        opts.droptol = 0.01;
        if (1)
            % t = find_shift(B);
            [Lk, Uk] = ilu(B-t*speye(n),opts);
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
    %%\
    prec.nC     = nC;
    prec.Nlev   = Nlev;
    prec.coefs1 = coefs1;
    prec.coefs2 = coefs2;
    prec.shift  = shift;
    prec.B      = B;
    prec.ITopts = ITopts;
     %%-------------------- get Lk Uk factors in certain cases.
end
end
end
