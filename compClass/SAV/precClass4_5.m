classdef precClass4_5 < handle 
%% original one -- matvec-only preconditioner 
%% 
 properties
     Nlev
     nC
     r         %% now a vector of length Nlev 
     coefs1     %% 2D array coefs1(:,k) = coefs for circle k
     coefs2     %% NOT USED - FOR COMPATIBILITY ONLY
     ITopts    %% parameters for solvers 
     shift    %% 2D array shift(:,k) = shifts for circle k
     B  
 end

methods
  function prec = precClass4_5(B,nC,r,ITopts)
  %%-------------------- CONSTRUCTOR --------------------
  %% nC     = number of countour points [upper] same for each circle
  %% r      = r(1:nlev) = radii of the circles
  %% ITopts = parameters for solver 
  %%-------- nC = number of poles on the upper half plane we need 2nC for
  %%  non-Hermitian case               
  %% construct preconditoner class. 
     n = size(B,1);
     %%-------------------- Part 4: Compute  weights + nodes 
     %%-------- Get nodes (z) and weights (om) on a unit circle 
     [z, om] = contQuad(nC,2);
     z = z(:); om = om(:);
     %%-------- Get rational approximation to 1/z outside circle 
     %% 1/z \approx \sum omega(k)/(sigma(k)+c\sigma(k)^2)
     %% *(1/(1/sigma(k)+c-z))                          
     %% Nlev   = number of levels
     Nlev   = length(r);
     sgn = 1;
     for lev = 1:Nlev
         %%-------------------- for each circle do:
         rj     = r(lev);    %% radius of jth circle 
         c      = rj*sgn*1.05;
         sgn    = -sgn;
         sh     = c + rj*z; 
         co2    = (rj * om) ./ sh;
         coefs1(:,lev) = [co2; conj(co2)];
         shift(:,lev)  = [sh; conj(sh)];
     end
    %%----------------    find a good complex shift for solving
    prec.nC     = nC;
    prec.Nlev   = Nlev;
    prec.coefs1 = coefs1;
    %% for compabitibility -- not used:
    prec.coefs2 = coefs1;
    prec.shift  = shift;
    prec.B      = B;
    prec.ITopts = ITopts;
     %%-------------------- get Lk Uk factors in certain cases.
end
end
end
