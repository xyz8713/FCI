function  A = A2d_ISO( freq, FS, Nz, Nx, N_delta, sigma0, label, vel)
% 2D acoustic Helmholtz matrix

% simpliefied by Xiao
hx = 1/Nx;
hz = 1/Nz;
rx = hz/hx;
tos = 1;


% build up two PML decaying function
Sz = ones(1,Nz+2*tos);
Sx = ones(1,Nx+2*tos);
for ii = 1:N_delta+tos
    tmp = sigma0/sqrt(-1)/freq*(1+cos(pi*ii/(N_delta+tos)))/2;
    Sx(Nx+2*tos-ii+1) = Sx(Nx+2*tos-ii+1) + tmp;
    Sz(Nz+2*tos-ii+1) = Sz(Nz+2*tos-ii+1) + tmp;
    Sx(ii) = Sx(Nx+2*tos-ii+1);
    if FS == 0
        Sz(ii) = Sz(Nz+2*tos-ii+1);
    end 
end

% A generation
len = 9*Nz*Nx - 6*((Nx-2)+(Nz-2)) - 20;
row = zeros(1,len);
col = zeros(1,len);
val = zeros(1,len);
counter = 0;

% No use to follow the tree sequence
for xx = 1:Nx
    for zz = 1:Nz

        row_tmp = label(zz,xx)*ones(3,3);
        col_tmp = zeros(3,3);
        rho_tmp = ones(3,3);
        
        for jj = -1:1
            for ii = -1:1
                if zz+ii >= 1 && zz+ii <= Nz && xx+jj >= 1 && xx+jj <= Nx
                    col_tmp(ii+2,jj+2) = label(zz+ii,xx+jj);
                end
            end
        end
               
        % value preprocessing
        Czz = -1;
        Cxx = -1 * rx^2;
        Czx =  0 * rx;
        
        % value generating
        val_tmp = zeros(3,3);
        
        % 1. d^2/dz^2
        tmp  = Czz * rho_tmp(2,2) / Sz(zz+tos);
        tmp1 = (1/(rho_tmp(1,2)*Sz(zz-1+tos)) + 1/(rho_tmp(2,2)*Sz(zz+tos)))/2;
        tmp2 = (1/(rho_tmp(3,2)*Sz(zz+1+tos)) + 1/(rho_tmp(2,2)*Sz(zz+tos)))/2;
        val_tmp(1,2) = val_tmp(1,2) + tmp * tmp1;
        val_tmp(2,2) = val_tmp(2,2) - tmp *(tmp1 + tmp2);
        val_tmp(3,2) = val_tmp(3,2) + tmp * tmp2;

        % 2. d^2/dx^2
        tmp  = Cxx * rho_tmp(2,2) / Sx(xx+tos);
        tmp1 = (1/(rho_tmp(2,1)*Sx(xx-1+tos)) + 1/(rho_tmp(2,2)*Sx(xx+tos)))/2;
        tmp2 = (1/(rho_tmp(2,3)*Sx(xx+1+tos)) + 1/(rho_tmp(2,2)*Sx(xx+tos)))/2;
        val_tmp(2,1) = val_tmp(2,1) + tmp * tmp1;
        val_tmp(2,2) = val_tmp(2,2) - tmp *(tmp1 + tmp2);
        val_tmp(2,3) = val_tmp(2,3) + tmp * tmp2;

        % 3. d^2/dzdx
        tmp = Czx * rho_tmp(2,2) / (4*Sz(zz+tos)*Sx(xx+tos));
        val_tmp(1,1) = val_tmp(1,1) + tmp * (rho_tmp(2,1)+rho_tmp(1,2))/(2*rho_tmp(2,1)*rho_tmp(1,2));
        val_tmp(3,1) = val_tmp(3,1) - tmp * (rho_tmp(2,1)+rho_tmp(3,2))/(2*rho_tmp(2,1)*rho_tmp(3,2));
        val_tmp(1,3) = val_tmp(1,3) - tmp * (rho_tmp(1,2)+rho_tmp(2,3))/(2*rho_tmp(1,2)*rho_tmp(2,3));
        val_tmp(3,3) = val_tmp(3,3) + tmp * (rho_tmp(3,2)+rho_tmp(2,3))/(2*rho_tmp(3,2)*rho_tmp(2,3));

        % mass term
        val_tmp(2,2) = val_tmp(2,2) - (hz*freq/vel(zz,xx))^2;

        for jj = 1:3
            for ii = 1:3
                if col_tmp(ii,jj) ~= 0
                    counter = counter+1;
                    row(counter) = row_tmp(ii,jj);
                    col(counter) = col_tmp(ii,jj);
                    val(counter) = val_tmp(ii,jj);
                end
            end
        end 
            
    end
end

A = sparse(row,col,val);

return;
end