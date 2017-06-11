%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Complex Non-Hermitian rational function preconditioners 3 %%%%%%%%%
%%%%%%%            Apply PA^{-1} with inner deflation             %%%%%%%%%
%%%%%%%                 Yuanzhe Xi, 02/29/2016
%%%%%%%                 %%%%%%%%%
%%%%%%%                 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

%% Part 1: Setup params for Krylov iteration
ITopts.tolIts  = 0.001; %% Outer- tolerance for stopping
ITopts.maxits  = 100;    %% Outer = max its
ITopts.outputG = 1;      %% print info during iteration
ITopts.im     =  30;     %% Outer=krylov subspace dimen.
%% Part 2: Setup the test matrix A and the right hand side b

n = 200;
d = 0.5+[0:199];
A = diag(d);
n = size(A,1);
    %%-------------------- -A = Lapl/
    h = 7.0;
    B = A-h*eye(n);
    %%-------------------- rhs
    rhs = ones(n,1);
%%-------------------- initial guess
sol0 = zeros(n,1);
nC1  =  32; 
%%r1  = [50.0:-30:20] ;

r1  = [60; 20];
%% r1 = 15.0;
%%r1  = [10 4 1];
%%%%r1  = [60; 30; 15; 5];
%%-------------------- set parameters for inner and outer
%%                     iterations in preconditioner.
%% note inner refers to the   A\inv P   part.
PREopts.outer   = 1;   %% (I-P)A\inv part see ratPrec4_x for explanations 
PREopts.inner   = 0;   %% PA\inv part. see ratPrecD3_x for explanations
                       %% default = 0 --> noprec.                      
PREopts.tolIts  =      0.01 ;    %% tolerance for stopping *inner* 
PREopts.maxits  =        20;     %% max its --  for inner
PREopts.outputG =         0;     %% print info during inner iteration
PREopts.im      =         5;     %% krylov subspace dimen. for inner
                                 %% gmres ( P A\inv) -
PREopts.imFOM   =        50;     %% krylov subspace dimen. for
                                 %% projector (I-P) -- FOM is used 
PREopts.Nvec    =         0;     %% number of deflation
                                 %% vectors. always zero 
PREopts.imSubs  =  30;     %% Outer=krylov subspace dimen.

%%-------------------- build preconditioner -
pre = precClass4_6(B, nC1, r1,PREopts); 
%%-------------------- iterate
tic;
[sol1,res2,its2] = fgmrez(B,pre,'ratPrec4_6', rhs, sol0,ITopts) ;
%%-------------------- plot
semilogy([0:its2],res2,'linestyle','-.','marker','v','LineWidth', ...
    2,'color','b')
t3 = toc;
fprintf(1,'iteration time is %f\n',t3);

%%------------------- 
