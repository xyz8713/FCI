function [A,rhs,label] = mywarp(nx,ny,nz,f)
if (nz >1)
    % paramter list:
    FS       = 1;    % 0: PML   1: free surface
    Nx       = nx;
    Ny       = ny;
    Nz       = nz;
    Lvl      = 7+(nx/40)*3;
    Lswitch  = 0;
    tol      = 1.0e-4;
    N_delta  = 6;
    freq     = 2*pi*f;
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

        
    % %%%%%%%%%%%%%%%%%  multiple RHS solver  %%%%%%%%%%%%%%%%
    [X,Y,Z] = meshgrid(1:Ny,1:Nx,1:Nz);
    rhs  = 1.0e-15 + exp(-((X-sourceX).^2+(Y-sourceY).^2+(Z-sourceZ).^2)/8);
    rhs(label) = rhs;
    rhs = reshape(rhs,Nx*Ny*Nz,1);  
else
    % paramter list:
    FS       = 0;    % 0: PML;  1: free surface
    Nz       = nx;
    Nx       = ny;
    Lvl      = 9+nx/40*3;
    Lswitch  = 0;
    tol      = 1.0e-4;
    N_delta  = 10;
    freq     = 2*pi*f;
    sigma0   = 2*pi*70;
    sourceZ  = round(Nz*0.5);
    sourceX  = round(Nx*0.5);
    % model paramters
    Vp       = ones(Nz,Nx);
    
    display('----- all parameters, velocity and density models have been built up -----');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%  1. nested dissection preprocessing  %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % a quick check of sampling
    fprintf('number of wavelengths in the domain: %.1f\n',1/min(Vp(:))*(freq/2/pi) );
    fprintf('points per wavelengths: %.1f\n', min([Nz,Nx]) *min(Vp(:))/(freq/2/pi) );
    tic;
    [ TT, Tlabel, NB, label, hssTT, hssTL, rg ] = NDpre2D( 1, 1, Lvl, Lswitch, Nz, Nx );
    Tpre = toc;
    display('----- preprocessing finished: altogether 7 objects have been created: TT, TL, pts, Tlabel, label, mt and NB -----');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%  2. matrix A generating  %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    A = A2d_ISO( freq, FS, Nz, Nx, N_delta, sigma0, label, Vp );
    TA = toc;
    
    display('----- full matrix A has been generated -----');
    % %%%%%  generating the RHS  %%%%%
    [X,Z] = meshgrid(1:Nx,1:Nz);
    
    rhs = 1.0e-15 + exp(-((X-sourceX).^2+(Z-sourceZ).^2)/8);
    rhs(label) = rhs;
    rhs = reshape(rhs,Nz*Nx,1);
    
end
end