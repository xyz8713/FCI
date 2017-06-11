   function [sol1, sol2] = ffom_mp(prec,rhs,ITSopt)
%% function [sol1, sol2] = ffom_mp(prec,rhs,ITSopt) 
%% barebone arnoldi multiple-pole solver
   %% Type 0 SOlve Ax = b only  --> sol1
   %% Type 1 do (I-P)A\inv b  and  A\inv P b --> sol1 and sol2
   %% Type 2 do (I - P)A \inv b only  --> sol1
   %% Type 3 do P b only  -> sol2
   %%
nC  = prec.nC;
A   = prec.B;
n   = size(A,1);
lev = ITSopt.lev;
im  = ITSopt.imFOM;
%%
coefs1 = prec.coefs1(:,lev);
coefs2 = prec.coefs2(:,lev);
shift = prec.shift(:,lev);
Type = ITSopt.Type;
outputG = ITSopt.outputG;
n = size(A,1)  ;
%
%-------------------- main loop -- no restarts
%
vv(1:n,1) = rhs ;
%%
ro = norm(vv(1:n,1))  ;
t = 1.0/ro;
vv(1:n,1)= vv(1:n,1) * t  ;
%%-------------------- print its/residual info
if (outputG)
   fprintf(1,' its %d  res %e \n',its,ro);
end
%-------------------- arnoldi loop
for i=1:im 
    %%-------------------- for first vectors Use vv vector
    z = A*vv(:,i) ;
    %%--------------------       modified GS  ;
    for j =1:i
        t  = vv(1:n,j)'*z  ; 
        hh(j,i) = t  ;
        z  = z - t*vv(1:n,j)  ;
    end
    t =norm(z,2)  ;
    i1 = i + 1 ;
    hh(i1,i) = t  ;
    if (t  ~= 0.0)
        t = 1.0 / t  ;
        vv(1:n,i1) = z*t  ;
    end 
    %%
end    %% end of arnoldi while (im) loop
       
%       now compute solutions. first solve upper triangular system.  ;
%
sol1 = zeros(n,1);
sol2 = zeros(n,1);
hh = hh(1:im,1:im);
%%-------------------- type == 0 -> solve one system - return it
%%                     corresponds to sol1 ~~ A\inv b
if (Type == 0)
    y  = hh \ eye(im,1) ;
    z = vv(:,1:im)*y;
    sol1 = sol1 + ro*z;
    return
end
%%-------------------- else type==1 or == 2
for ii = 1:2*nC
    H1 = hh - shift(ii)*eye(im,im);
    %%-------------------- type 1 : solve both.. 
    %%                     corresponds to method == 1 in ratPrec4_0/
    %%-------------------- TYPE 1 AND TYPE 2:
    if (Type <= 2)
        y  =  H1 \ eye(im,1) ;
        z  = vv(:,1:im)*y;
        s1 = coefs1(ii)*ro;
        sol1 = sol1 + s1*z;
        %%-------------------- TYPE 1: add second solution
    end
    if (Type ~= 2)
%%-------------------- solve second set of systems first
       s2 = coefs2(ii)*ro;
       sol2 = sol2 + s2*z;
   end
end
