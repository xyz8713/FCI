%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      Mixed grid compact Finite Difference stencil       %%%%%%%%
%%%%%%% 3D Helmholtz matrix after Nested Dissection reordering  %%%%%%%%
%%%%%%%                Shen Wang, 12/16/2009                    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% profile on;
set(0,'defaultaxesfontsize',20);

% paramter list:
FS       = 1;    % 0: PML   1: free surface
Nx       = 40;
Ny       = 40;
Nz       = 40;
Lvl      = 9;
Lswitch  = 0;
tol      = 1.0e-4;
N_delta  = 6;
freq     = 2*pi*5;
sigma0   = 2*pi*70;
sourceX  = round(Nx*0.5);
sourceY  = round(Ny*0.5);   
sourceZ  = round(Nz*0.5);
% model paramters
Vp       = ones(Nx,Ny,Nz);

display('----- all parameters, velocity and density models have been built up -----');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  1. nested dissection preprocessing  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
[ TT, Tlabel, label, NB, hssTT, hssTL, rg ] = NDpre3D( 1, 1, Lvl, Lswitch, Nx, Ny, Nz );
Tpre = toc;
display('----- preprocessing finished: altogether 7 objects have been created: TT, TL, pts, Tlabel, label, mt and NB -----');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  2. matrix A generating  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('number of wavelengths in the domain: %.1f\n',1/min(Vp(:))*(freq/2/pi) );
fprintf('points per wavelengths: %.1f\n', min([Nx,Ny,Nz]) *min(Vp(:))/(freq/2/pi) );
tic;
A = A3d_ISO(freq, FS, Nx, Ny, Nz, N_delta, sigma0,label, Vp);   
TA = toc;

display('----- full matrix A has been generated -----');
% figure; spy(A); title(['after nested dissection: Level = ',num2str(L)]); pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  3. extend-add based factorization: using stack  %%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tic;
% [ Q1, Q2, L, DL, U, B, W, V, FlopsFac, storage ] = factorization( A, tol, Lswitch, TT, Tlabel, NB, hssTT, hssTL, rg );
% Tfac = toc;
%     
%     
% %%%%%%%%%%%%%%%%%  multiple RHS solver  %%%%%%%%%%%%%%%%
% [X,Y,Z] = meshgrid(1:Ny,1:Nx,1:Nz);
%         x  = 1.0e-15 + exp(-((X-sourceX).^2+(Y-sourceY).^2+(Z-sourceZ).^2)/8);
%         x(label) = x;
%         x = reshape(x,Nx*Ny*Nz,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  4. system solution of multiple RHS   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% [ x, FlopsSol ] = mfsolvePar( x, Lswitch, TT, Tlabel, NB, hssTT, hssTL, rg, Q1, Q2, L, DL, U, B, W, V );   % system solution
% Tsol = toc;
% 
% 
%         x = reshape(x,Nx,Ny,Nz);
%         x = x(label);
%         figure; slice(real(x),round(Nx/2),round(Ny/2),round(Nz/2)); shading interp; axis image; colorbar;


%display('----- solving Ax = b with local dense LU factor ends -----');

% display(['relative residule norm(Ax-b)/norm(b) is: ', num2str(norm(A*x-b)/norm(b))]);
% display('----------------------  Conclusions  -------------------------');
% display(['-----  Nx = ',num2str(Nx-2),', Ny = ',num2str(Ny-2),', Nz = ',num2str(Nz-2),', order = ',num2str(size(A,1)),' -----']);
% display(['-----  Time for preprocessing : ',num2str(Tpre),' -----']);
% display(['-----  Time for generating A  : ',num2str(TA),' -----']);
% display(['-----  Time for factorizing A : ',num2str(Tfac),' -----']);
% display(['-----  Time for solution      : ',num2str(Tsol),' -----']);
% fprintf('-----  Factorization flops is: %e -----\n',FlopsFac);
% fprintf('-----  Solution flops is:      %e -----\n',FlopsSol);
% fprintf('-----  Storage is:             %e -----\n',storage);
% display('--------------------------------------------------------------')

% profile viewer;
