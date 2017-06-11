%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Complex Non-Hermitian rational function preconditioners 3 %%%%%%%%%
%%%%%%%            Apply PA^{-1} with inner deflation             %%%%%%%%%
%%%%%%%                 Yuanzhe Xi, 02/29/2016
%%%%%%%                 %%%%%%%%%
%%%%%%%                 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;
clear all;
addpath ('./Helmholtz')
Pb = 'L';

%% Part 1: Setup params for Krylov iteration
ITopts.tolIts  = 1.e-2; %% tolerance for stopping
ITopts.maxits  = 10;    %% max its
ITopts.outputG = 1;      %% print info during iteration
ITopts.im     =  10;     %% krylov subspace dimen.
%% Part 2: Setup the test matrix A and the right hand side b
if (Pb == 'L' || Pb == 'l')
    nx = 20;
    ny = 20;
    nz = 20;
    A = fd3d(nx,ny,nz, 1, 0.0, 0.0, 0.0);
    n = size(A,1);
    %%-------------------- -A = Lapl/
    h = 0.5;
    B = A-h*speye(n);
    %%-------------------- rhs
    rhs = B*randn(n,1);
elseif (Pb == 'H')
    nx = 80;
    ny = 80;
    nz = 80;
    f = nx/40*5;
    [B,rhs,label] = mywarp(nx,ny,nz,f); %%Helmholtz equation nx=ny=nz
    %%f is proportional to nx to fix 8 points per wavelength
    %%nx = 80  f = 10hz nx = 40 f = 5hz
else 
    disp('problem error');
    exit;
end
n = size(B,1);
%%-------------------- initial guess
sol0 = randn(n,1);
nC  = 4; r  = 30.0;
pre = precClass(B, nC, r);

tic;
%%--------------------
ITopts.tolIts = 1.e-07;
[sol1,res2,its2] = fgmrez(B,pre,'ratPrecD32_c', rhs, sol0,ITopts) ;
semilogy([0:its2],res2,'linestyle','-.','marker','v','LineWidth', ...
    2,'color','b')
t3 = toc;
fprintf(1,'iteration time is %f\n',t3);

