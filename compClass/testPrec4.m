%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Complex Non-Hermitian rational function preconditioners 3 %%%%%%%%%
%%%%%%%            Apply PA^{-1} with inner deflation             %%%%%%%%%
%%%%%%%                 Yuanzhe Xi, 02/29/2016
%%%%%%%                 %%%%%%%%%
%%%%%%%                 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% three cases: Helmholtz, shifted laplacean, external.
close all;
clear all;
%% three cases: 'H'elmholtz, shifted 'L'aplacean, or 'E'xternal.
Pb = 'H';

%% Part 1: Setup params for Krylov iteration
ITopts.tolIts  = 1.e-07; %% Outer- tolerance for stopping
ITopts.maxits  = 100;    %% Outer = max its
ITopts.outputG = 1;      %% print info during iteration
ITopts.im     =  30;     %% Outer=krylov subspace dimen.
%% Part 2: Setup the test matrix A and the right hand side b
if (Pb == 'L' || Pb == 'l')
    nx = 30;
    ny = 30 ;
    nz = 30;
    A = fd3d(nx,ny,nz,1.0, 0.0, 0.0, 0.0);
    n = size(A,1);
    %%-------------------- -A = Lapl/
    h = 1.05;
    %%h = 0;
    B = A-h*speye(n);
    %%-------------------- rhs
    rhs = B*randn(n,1);
elseif (Pb == 'H')
    addpath ('./Helmholtz')
    nx = 40;
    ny = 40;
    nz = 40;
    f  = (nx/40)*8;
    [B,rhs,label] = mywarp(nx,ny,nz,f); %%Helmholtz equation nx=ny=nz
    %%f is proportional to nx to fix 8 points per wavelength
    %%nx = 80  f = 10hz nx = 40 f = 5hz
elseif (Pb == 'E')
    provideMatrix;
else 
    disp('problem error');
    return;
end
n = size(B,1);
%%-------------------- initial guess
sol0 = randn(n,1);
nC1  = 8 ; 
%%r1  = [ 50; 30];
%%--> r1  = [40; 25; 5];
%% %% for class 4_6 only:
%%
r1  = [10; 10; 5; 5; 2; 2; 1; 1];
r1  = [100; 60; 60; 20; 40; 5];
r1  = [10; 5];
r1  = [60; 40; 20; 5 ];
r1 =  [60];
%%r1  = [100; 2; 60; 2; 20; 2; 5; 2];
%%r1 =  50.0*ones(10,1);
%r1  = [60; 30; 15; 5];
%%-------------------- set parameters for inner and outer
%%                     iterations in preconditioner.
%% note inner refers to the   A\inv P   part.
PREopts.outer   = 2;    %% (I-P)A\inv part see ratPrec4_x for explanations 
PREopts.inner   = 0;    %% PA\inv part. see ratPrecD3_x for explanations
                        %% default = 0 --> noprec.                      
PREopts.tolIts  =  1.e-03;    %% tolerance for stopping *inner* 
PREopts.maxits  =      20;    %% max its --  for inner
PREopts.outputG =      0;    %% print info during inner
                                 %% iteration 
PREopts.im      =      10;     %% krylov subspace dimen. for inner
                               %% gmres ( P A\inv) -
PREopts.imFOM   =      20;   %% krylov subspace dimen. for
                               %% projector (I-P) -- FOM is used 
PREopts.Nvec    =       0;    %% number of deflation
                                %% vectors. always zero -- useful
                                %% for large values of PREopts.im 
PREopts.imSubs   =      0;   %% Dimension of approximate
                               %% eigenspace for method 4_6 only 

%%-------------------- build preconditioner -
pre = precClass4_4(B, nC1,r1,PREopts); 
%%-------------------- iterate
tic;
[sol1,res2,its2] = fgmrez(B,pre,'ratPrec4_4', rhs, sol0,ITopts) ;
%%-------------------- plot
semilogy([0:its2],res2,'linestyle','-.','marker','v','LineWidth', ...
    2,'color','b')
t3 = toc;
fprintf(1,'iteration time is %f\n',t3);

%%------------------- 
