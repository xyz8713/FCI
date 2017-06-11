%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     plot matrix pattern      %%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultaxesfontsize',20);

% control parameters
Nx       = 20;
Ny       = 20;
Nz       = 20;
Lvl      = 6;   
lw       = 0.5;    % LineWidth 

% model parameters
vel_pz   = 4000*ones(Nx,Ny,Nz);
vel_sz   = 2500*ones(Nx,Ny,Nz);
rho      =  1.0*ones(Nx,Ny,Nz);
epsilon  =  0.3*ones(Nx,Ny,Nz);
delta    =  0.0*ones(Nx,Ny,Nz);
theta    =   30*ones(Nx,Ny,Nz);
phi      =   45*ones(Nx,Ny,Nz);

% unimportant parameters
Lswitch = 0;
N_delta = 0;
freq = 2*pi*10;
sigma0 = 2*pi*60; 
hx = 50;
hy = 50;
hz = 50;

tos = 1;
[ TT, Tlabel1, label, NB, hssTT, hssTL, rg ] = NDpre3D( 1, tos, Lvl, Lswitch, Nx, Ny, Nz );
A1 = A3d_ELP( freq, 1, 0, Nx, Ny, Nz, N_delta, sigma0, hx, hy, hz, label, vel_pz, rho, epsilon, theta, phi );
A1 = real(A1); clear TT label NB hssTT hssTL rg;

% tos = 1;
% [ TT, TlabelE, label, NB, hssTT, hssTL, rg ] = NDpre3D( 3, tos, Lvl, Lswitch, Nx, Ny, Nz );
% E1 = E3d_ISO( freq, 3, tos, 0, Nx, Ny, Nz, N_delta, sigma0, hx, hy, hz, label, vel_pz, vel_sz, rho );
% E1 = real(E1); clear label NB hssTT hssTL rg;

tos = 2;
[ TT, Tlabel2, label, NB, hssTT, hssTL, rg ] = NDpre3D( 1, tos, Lvl, Lswitch, Nx, Ny, Nz );
A2 = A3d_TTI( freq, 2, 0, Nx, Ny, Nz, N_delta, sigma0, hx, hy, hz, label, vel_pz, rho, epsilon, delta, theta, phi );
A2 = real(A2); clear label NB hssTT hssTL rg;

% clear unwanted variables
clear vel_pz vel_sz rho epsilon delta theta phi Lswitch N_delta freq sigma0 hx hy hz;

% plotting: tos = 1
figure; spy(A1); title(['t = 1, level = ',num2str(Lvl)]); hold on;
for Ntmp = 1:length(TT)
    plot( Tlabel1(Ntmp):Tlabel1(Ntmp+1)-1, Tlabel1(Ntmp)*ones(1,Tlabel1(Ntmp+1)-Tlabel1(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
    plot( Tlabel1(Ntmp)*ones(1,Tlabel1(Ntmp+1)-Tlabel1(Ntmp)), Tlabel1(Ntmp):Tlabel1(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
    plot( Tlabel1(Ntmp):Tlabel1(Ntmp+1)-1, (Tlabel1(Ntmp+1)-1)*ones(1,Tlabel1(Ntmp+1)-Tlabel1(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
    plot( (Tlabel1(Ntmp+1)-1)*ones(1,Tlabel1(Ntmp+1)-Tlabel1(Ntmp)), Tlabel1(Ntmp):Tlabel1(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
    kidL = TT{Ntmp}.kids(1);
    if kidL ~= 0
        while TT{kidL}.kids(1) ~= 0
            kidL = TT{kidL}.kids(1);
        end
        plot( Tlabel1(kidL):Tlabel1(Ntmp)-1, Tlabel1(Ntmp)*ones(1,Tlabel1(Ntmp)-Tlabel1(kidL)), 'r--', 'LineWidth', lw ); hold on;
        plot( Tlabel1(Ntmp)*ones(1,Tlabel1(Ntmp)-Tlabel1(kidL)), Tlabel1(kidL):Tlabel1(Ntmp)-1, 'r--', 'LineWidth', lw ); hold on;
        plot( Tlabel1(kidL):Tlabel1(Ntmp)-1, (Tlabel1(Ntmp+1)-1)*ones(1,Tlabel1(Ntmp)-Tlabel1(kidL)), 'r--', 'LineWidth', lw ); hold on;
        plot( (Tlabel1(Ntmp+1)-1)*ones(1,Tlabel1(Ntmp)-Tlabel1(kidL)), Tlabel1(kidL):Tlabel1(Ntmp)-1, 'r--', 'LineWidth', lw ); hold on;
        plot( Tlabel1(kidL)*ones(1,Tlabel1(Ntmp+1)-Tlabel1(Ntmp)), Tlabel1(Ntmp):Tlabel1(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
        plot( Tlabel1(Ntmp):Tlabel1(Ntmp+1)-1, Tlabel1(kidL)*ones(1,Tlabel1(Ntmp+1)-Tlabel1(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
        plot( (Tlabel1(Ntmp)-1)*ones(1,Tlabel1(Ntmp+1)-Tlabel1(Ntmp)), Tlabel1(Ntmp):Tlabel1(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
        plot( Tlabel1(Ntmp):Tlabel1(Ntmp+1)-1, (Tlabel1(Ntmp)-1)*ones(1,Tlabel1(Ntmp+1)-Tlabel1(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
    end
end
hold off;

% % plotting: E1
% figure; spy(E1); title(['elastic matrix: level = ',num2str(Lvl)]); hold on;
% for Ntmp = 1:length(TT)
%     plot( TlabelE(Ntmp):TlabelE(Ntmp+1)-1, TlabelE(Ntmp)*ones(1,TlabelE(Ntmp+1)-TlabelE(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
%     plot( TlabelE(Ntmp)*ones(1,TlabelE(Ntmp+1)-TlabelE(Ntmp)), TlabelE(Ntmp):TlabelE(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
%     plot( TlabelE(Ntmp):TlabelE(Ntmp+1)-1, (TlabelE(Ntmp+1)-1)*ones(1,TlabelE(Ntmp+1)-TlabelE(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
%     plot( (TlabelE(Ntmp+1)-1)*ones(1,TlabelE(Ntmp+1)-TlabelE(Ntmp)), TlabelE(Ntmp):TlabelE(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
%     kidL = TT{Ntmp}.kids(1);
%     if kidL ~= 0
%         while TT{kidL}.kids(1) ~= 0
%             kidL = TT{kidL}.kids(1);
%         end
%         plot( TlabelE(kidL):TlabelE(Ntmp)-1, TlabelE(Ntmp)*ones(1,TlabelE(Ntmp)-TlabelE(kidL)), 'r--', 'LineWidth', lw ); hold on;
%         plot( TlabelE(Ntmp)*ones(1,TlabelE(Ntmp)-TlabelE(kidL)), TlabelE(kidL):TlabelE(Ntmp)-1, 'r--', 'LineWidth', lw ); hold on;
%         plot( TlabelE(kidL):TlabelE(Ntmp)-1, (TlabelE(Ntmp+1)-1)*ones(1,TlabelE(Ntmp)-TlabelE(kidL)), 'r--', 'LineWidth', lw ); hold on;
%         plot( (TlabelE(Ntmp+1)-1)*ones(1,TlabelE(Ntmp)-TlabelE(kidL)), TlabelE(kidL):TlabelE(Ntmp)-1, 'r--', 'LineWidth', lw ); hold on;
%         plot( TlabelE(kidL)*ones(1,TlabelE(Ntmp+1)-TlabelE(Ntmp)), TlabelE(Ntmp):TlabelE(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
%         plot( TlabelE(Ntmp):TlabelE(Ntmp+1)-1, TlabelE(kidL)*ones(1,TlabelE(Ntmp+1)-TlabelE(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
%         plot( (TlabelE(Ntmp)-1)*ones(1,TlabelE(Ntmp+1)-TlabelE(Ntmp)), TlabelE(Ntmp):TlabelE(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
%         plot( TlabelE(Ntmp):TlabelE(Ntmp+1)-1, (TlabelE(Ntmp)-1)*ones(1,TlabelE(Ntmp+1)-TlabelE(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
%     end
% end
% hold off;

% plotting: tos = 2
figure; spy(A2); title(['t = 2, level = ',num2str(Lvl)]); hold on;
for Ntmp = 1:length(TT)
    plot( Tlabel2(Ntmp):Tlabel2(Ntmp+1)-1, Tlabel2(Ntmp)*ones(1,Tlabel2(Ntmp+1)-Tlabel2(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
    plot( Tlabel2(Ntmp)*ones(1,Tlabel2(Ntmp+1)-Tlabel2(Ntmp)), Tlabel2(Ntmp):Tlabel2(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
    plot( Tlabel2(Ntmp):Tlabel2(Ntmp+1)-1, (Tlabel2(Ntmp+1)-1)*ones(1,Tlabel2(Ntmp+1)-Tlabel2(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
    plot( (Tlabel2(Ntmp+1)-1)*ones(1,Tlabel2(Ntmp+1)-Tlabel2(Ntmp)), Tlabel2(Ntmp):Tlabel2(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
    kidL = TT{Ntmp}.kids(1);
    if kidL ~= 0
        while TT{kidL}.kids(1) ~= 0
            kidL = TT{kidL}.kids(1);
        end
        plot( Tlabel2(kidL):Tlabel2(Ntmp)-1, Tlabel2(Ntmp)*ones(1,Tlabel2(Ntmp)-Tlabel2(kidL)), 'r--', 'LineWidth', lw ); hold on;
        plot( Tlabel2(Ntmp)*ones(1,Tlabel2(Ntmp)-Tlabel2(kidL)), Tlabel2(kidL):Tlabel2(Ntmp)-1, 'r--', 'LineWidth', lw ); hold on;
        plot( Tlabel2(kidL):Tlabel2(Ntmp)-1, (Tlabel2(Ntmp+1)-1)*ones(1,Tlabel2(Ntmp)-Tlabel2(kidL)), 'r--', 'LineWidth', lw ); hold on;
        plot( (Tlabel2(Ntmp+1)-1)*ones(1,Tlabel2(Ntmp)-Tlabel2(kidL)), Tlabel2(kidL):Tlabel2(Ntmp)-1, 'r--', 'LineWidth', lw ); hold on;
        plot( Tlabel2(kidL)*ones(1,Tlabel2(Ntmp+1)-Tlabel2(Ntmp)), Tlabel2(Ntmp):Tlabel2(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
        plot( Tlabel2(Ntmp):Tlabel2(Ntmp+1)-1, Tlabel2(kidL)*ones(1,Tlabel2(Ntmp+1)-Tlabel2(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
        plot( (Tlabel2(Ntmp)-1)*ones(1,Tlabel2(Ntmp+1)-Tlabel2(Ntmp)), Tlabel2(Ntmp):Tlabel2(Ntmp+1)-1, 'r--', 'LineWidth', lw ); hold on;
        plot( Tlabel2(Ntmp):Tlabel2(Ntmp+1)-1, (Tlabel2(Ntmp)-1)*ones(1,Tlabel2(Ntmp+1)-Tlabel2(Ntmp)), 'r--', 'LineWidth', lw ); hold on;
    end
end
hold off;