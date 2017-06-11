classdef precClass4_6 < handle 
%% original one -- matvec-only preconditioner 
%% 
 properties
     nC
     r         %% Only one circle
     coefs1    %% coeffs for (I-P)
     coefs2    %% coeffs for P
     ITopts    %% various iterations parameters 
     shift     %% shifts == poles on circle. 
     B  
     V         %% invariant subspace basis
     H         %% representation.
     U         %% for deflation;
 end

methods
  function prec = precClass4_6(B,nC,r,ITopts)
  %%-------------------- CONSTRUCTOR --------------------
  %% nC     = number of countour points [upper] same for each circle
  %% r      = radius of the circle
  %% ITopts = parameters for solver 
  %%-------- nC = number of poles on the upper half plane we need 2nC for
  %%  non-Hermitian case               
  %% construct preconditoner class. 
     n = size(B,1);
     r = r(1:1);  %% <<<< safety -- use only one radius
     %%-------- intial center
     %% Part 4: Compute the weights and nodes for rational approximation
     %% to h                            
     %%-------- Get nodes (z) and weights (om) on a unit circle 
     [z, om] = contQuad(nC,2);
     z = z(:); om = om(:);
     %%-------- Get rational approximation to 1/z outside circle 
     %% 1/z \approx \sum omega(k)/(sigma(k)+c\sigma(k)^2)
     %% *(1/(1/sigma(k)+c-z))                          
     %%-------------------- for each circle do:
     sigma  = z/r;
     omega  = om/r;
     c      = max(real(z))*r-1.0;
     c = 0.0;
     co1    = -omega./(sigma + c*sigma.*sigma);
     sh     = 1./sigma +c;
     %%     shift2 = [shift2, conj(shift2)];
     co2    = conj(om*r);
     coefs1 = [co1; conj(co1)];
     coefs2 = [co2; conj(co2)];
     shift  = [sh; conj(sh)];
     %%
    %%----------------    find a good complex shift for solving
    d = abs(B)*ones(n,1)-abs(diag(B)); 
    d = diag(B)-d;
    t = 0.5*min(d);
    
    ITopts.shift = t;
   
    %% will make matrix B diagonally dominant. 
    %%---------------- construct basis of near null space 
    im = ITopts.imSubs;
    %%--------------------Arnoldi loop
    if (im == 0) 
        V = [];
        H = [];
    else
        [U S V] = svds(B,im,0);
        H = V'*B*V;
    end
    %%-------------------- 
    %%-------------------- pack into struct.
    prec.nC     = nC;
    prec.coefs1 = coefs1;
    prec.coefs2 = coefs2;
    prec.shift  = shift;
    prec.B      = B;
    prec.ITopts = ITopts;
    prec.V      = V;
    prec.U      = [];
    prec.H      = H;
     %%-------------------- get Lk Uk factors in certain cases.
end
end
end
