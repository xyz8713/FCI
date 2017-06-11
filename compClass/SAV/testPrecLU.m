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
Pb = 'H';

%% Part 1: Setup params for Krylov iteration
ITopts.tolIts  = 1.e-07; %% tolerance for stopping
ITopts.maxits  = 200;    %% max its
ITopts.outputG = 1;      %% print info during iteration
ITopts.im     =  30;     %% krylov subspace dimen.
%% Part 2: Setup the test matrix A and the right hand side b
if (Pb == 'L' || Pb == 'l')
    nx = 40;
    ny = 40;
    nz = 40;
    A = fd3d(nx,ny,nz, 1, 0.0, 0.0, 0.0);
    n = size(A,1);
    %%-------------------- -A = Lapl/
    h = 0.70;
    B = A-h*speye(n);
    %%-------------------- rhs
    rhs = B*randn(n,1);
elseif (Pb == 'H')
    nx = 80;
    ny = 80;
    nz = 80;
    f  = (nx/40)*8;

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
%% 
d = find_shiftd(B);

%%t = find_shift(B);

opts.type = 'ilutp';
opts.droptol = 0.01;  

tic  
[Lk, Uk,P] = ilu(B+spdiags(d,0,n,n),opts);
%%[Lk, Uk,P] = ilu(B+t*speye(n),opts);
(nnz(Lk)+nnz(Uk)-n)/nnz(B)
tp = toc;
 fprintf(1,'precon time is %f\n',tp);
PRE.L = Lk;
PRE.U = Uk;
PRE.P = P;
tic;
%%--------------------

[sol1,res2,its2] = fgmrez(B,PRE,'precLU', rhs, sol0,ITopts) ;
semilogy([0:its2],res2,'linestyle','-.','marker','v','LineWidth', ...
    2,'color','b')
t3 = toc;
fprintf(1,'iteration time is %f\n',t3);

ff = (nnz(Lk)+nnz(Uk)-n)/nnz(B);
fprintf(1,'fill-factor  is   %f\n',ff);

%%------------------- 
