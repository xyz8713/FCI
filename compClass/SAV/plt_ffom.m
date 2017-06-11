   function [sol1, sol2] = plt_ffom(prec,rhs,ITSopt)
%% function [sol1, sol2] = ffom_mp(prec,rhs,ITSopt) 
%% barebone arnoldi multiple-pole solver
   %%
nC = prec.nC;
A = prec.B;
n = size(A,1);
lev = ITSopt.lev;
%%
coefs1 = prec.coefs1(:,lev)
coefs2 = prec.coefs2(:,lev)
shift = prec.shift(:,lev);
Type = ITSopt.Type;
outputG = ITSopt.outputG;
n = size(A,1)  ;
%  ### DEBUG ONLY
%-------------------- main loop -- no restarts
%
ro =1.0;
%%-------------------- print its/residual info
%-------------------- arnoldi loop
%       now compute solutions. first solve upper triangular system.  ;
%
sol1 = zeros(n,1);
sol2 = zeros(n,1);
%%-------------------- else type==1 or == 2
hh = prec.B;

for ii = 1:2*nC
    H1 = hh - shift(ii)*eye(n,n);
    %%-------------------- type 1 : solve both.. 
    %%                     corresponds to method == 1 in ratPrec4_0/
    %%-------------------- TYPE 1 AND TYPE 2:
    z =  H1 \ rhs;
    s1 = coefs1(ii)*ro;
    sol1 = sol1 + s1*z;
    %%-------------------- TYPE 1: add second solution
    if (Type == 1)
%%-------------------- solve second set of systems first
       s2 = coefs2(ii)*ro;
       sol2 = sol2 + s2*z;
   end
end
