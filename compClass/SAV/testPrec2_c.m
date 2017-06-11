%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     Complex Non-Hermitian rational function preconditioners  2                                  %%%%%%%%
%%%%%%%           Apply PA^{-1} without FEAST                                                                             %%%%%%%%
%%%%%%%            Yuanzhe Xi, 02/01/2016                                                                                       %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
addpath ('./Helmholtz')
%% Part 1: Setup params for Krylov iteration
ITopts.tolIts  = 1.e-5; %% tolerance for stopping
ITopts.maxits  = 10;    %% max its
ITopts.outputG = 1;      %% print info during iteration
ITopts.im     =  10;     %% krylov subspace dimen.
%% Part 2: Setup the test matrix A and the right hand side b
if (0)
%     nx = 20;
%     ny = 20;
%     nz = 20;
%     A = fd3d(nx,ny,nz, 1, 0.0, 0.0, 0.0);
%     n = size(A,1);
%     %%-------------------- -A = Lapl/
%     h = 0.25;
%    B = A-h*speye(n);
    %%-------------------- rhs 
    load fem_hifreq_circuit;
    B = Problem.A;
    n = size(B,1);
    rhs = B*randn(n,1);
else
    nx = 40;
    ny = 40;
    nz = 40;
    f = nx/40*5;
    [B,rhs,label] = mywarp(nx,ny,nz,f); %%Helmholtz equation
    %%f is proportional to nx to fix 8 points per wavelength
    %%nx = 80  f = 10hz;   nx = 40 f = 5hz
end
n = size(B,1);
%%-------------------- initial guess
sol0 = randn(n,1);
%% Part 3: Setup params for the circle \Gamma
%%-------- radius
 
%% add timer
nC = 4;  r = 30;
pre = precClass(B,nC,r);

tic;
%%--------------------
[sol1,res2,its2] = fgmrez(B, pre,'ratPrecD32_c', rhs, sol0,ITopts) ;
semilogy([0:its2],res2,'linestyle','-.','marker','v','LineWidth', ...
    2,'color','b')
t3 = toc;
fprintf(1,'iteration time is %f\n',t3);
%%visualize the solution
if(0)
    sol1 = reshape(sol1,nx,ny,nz);
    x = sol1(label);
    figure; slice(real(x),round(nx/2),round(ny/2),round(nz/2)); shading interp; axis image; colorbar;
end
% %%visualize the exact solution
% figure(5)
% sol2 = B\rhs;
% if(1)
%     sol1 = reshape(sol2,nx,ny,nz);
%     x = sol2(label);
%     figure; slice(real(x),round(nx/2),round(ny/2),round(nz/2)); shading interp; axis image; colorbar;
% end
