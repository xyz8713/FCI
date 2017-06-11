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
     %%Z         %% invariant subspace basis
     H         %% representation.
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
    %%---------------- construct basis of eigenspace. 
    im = ITopts.imSubs;
    %%--------------------Arnoldi loop

    tmpst =struct('shift',shift,'coefs1',coefs1,'coefs2',coefs2,'nC',nC,'B',B);
    ITopts.Type = 3;
    ITopts.lev = 1;
    V = randn(n,1);
    [dummy, V] = ffom_mp(tmpst,V,ITopts);
    V = V/norm(V);

    for i=1:im 
%%-------------------- for first vectors Use U vector
          [dummy, v] = ffom_mp(tmpst,V(:,i),ITopts);
          %%Z(:,i) = z;
          %          v = B*z;
     %%--------------------       modified GS  ;
     for j =1:i
         t  = V(1:n,j)'*v  ; 
         H(j,i) = t  ;
         v  = v - t*V(1:n,j)  ;
    end
    t =norm(v,2)  ;
    
    i1 = i + 1 ;
    H(i1,i) = t  ;
    %% 
    if (t  ~= 0.0)
        t = 1.0 / t  ;
        V(:,i1) = v*t  ;
    end
    %%
end    %% end of arnoldi while (im) loop
    H = V'*(B*V);

    %%----------------    find a good complex shift for solving
     d = abs(B)*ones(n,1)-abs(diag(B)); 
     d = diag(B)-d;
     t = min(d);
     %% will make matrix B diagonally dominant. 
     ITopts.shift = 0.5*t;
     %%-------------------- pack into struct.
    prec.nC     = nC;
    prec.coefs1 = coefs1;
    prec.coefs2 = coefs2;
    prec.shift  = shift;
    prec.B      = B;
    prec.ITopts = ITopts;
    prec.V      = V;
    %%prec.Z      = Z;
    prec.H     = H; 
     %%-------------------- get Lk Uk factors in certain cases.
end
end
end
