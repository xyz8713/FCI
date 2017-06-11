  function [sol,res] = CGNE (A, x, rhs, iter, tol) 
%----------------------------------------------------------------------- 
% function [res,sol] = CGNE (A, x, rhs, iter, tol) 
%  cgne (A, sol, rhs, iter, tol) 
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
 p = A'*r; 

 ro1 = r'*r; 
 res(1) = sqrt(ro1);
 tol1 = tol*tol*ro1; 
%%-------------------- main iteration
 j = 0; 
 while (j < iter & ro1 > tol1) 
	 j = j+1; 
	 ro = ro1; 
	 alp = (r'*r) / (p' * p ) ;
	 sol = sol + alp * p ;
	 r = r - alp *(A*p);
	 ro1= r'*r;
         res(j+1) = sqrt(ro1);
	 bet = ro1 / ro ;
         w = A'*r;
	 p = w + bet * p; 
 end
         
