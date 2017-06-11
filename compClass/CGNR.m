  function [sol,res] = CGNR (A, x, rhs, iter, tol) 
%----------------------------------------------------------------------- 
% function [res,sol] = CGNR (A, x, rhs, iter, tol) 
%  cgnr (A, sol, rhs, iter, tol) 
% sol  = initial guess on entry -- app least-squ. sol
% on return 
% rhs  = right-hand-side, 
% iter = max number of iterations..
% tol  = tolerabce 
% res  returns the residuals for each step [for plotting
% purposes.] 
%----------------------------------------------------------------------- 
%%-------------------- initial vector and residual 
 sol = x; 
 r = rhs - A * sol;
 z = A'*r; 
 p = z ; 
 ro1 = z'*z; 
 res(1) = sqrt(ro1);
 tol1 = tol*tol*ro1; 
%%-------------------- main iteration
 j = 0; 
 while (j < iter & ro1 > tol1) 
	 j = j+1; 
	 ro = ro1; 
	 w = A * p; 
	 alp = (z'*z) / (w' * w ) ;
	 sol = sol + alp * p ;
	 r = r - alp * w; 
         z = A'*r;
	 ro1= z'*z;
         res(j+1) = sqrt(ro1);
	 bet = ro1 / ro ;
	 p = z + bet * p; 
 end
         
